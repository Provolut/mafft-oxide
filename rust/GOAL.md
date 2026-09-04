# Goal: prove correctness cheaply, then match the C's speed

**The constraint does not move: output must not change.** Every step below
keeps `tests/matrix.rs`, `verify.sh`, `verify_tbfast.sh` and `verify_e2e.sh`
green. An optimisation that changes one alignment is a bug, however fast.

## Where this starts

Measured 2026-09-04, `disttbfast -F` on `test/sample` (36 seqs):

| metric | C | port | ratio |
|---|---|---|---|
| instructions (callgrind Ir) | 764 M | 1298 M | **1.70x** |
| peak RSS | 25.1 MB | 16.1 MB | **0.64x** |
| wall clock (best of 5) | 108 ms | 131 ms | **1.21x** |

`cargo test`: 4.4 s, 358 bit-identity cases, 86.2% regions / 85.9% lines.

## 1. Testing — a floor to hold, not a milestone to reach

- **>= 85% region coverage from `cargo test` alone.** At 86.2%; do not regress.
- **Suite stays far under 5 min.** At 4.4 s. Keep it seconds: tiny inputs, no
  sleeps, no network. A slow suite stops being run.
- **Bit-identity is what the tests are for.** `tests/matrix.rs` replays cases
  whose expected bytes came from the real C via `gen_golden.py`; a fixture that
  changes is the port changing its answer. Add a feature -> add its argument
  vector to `gen_golden.py` and regenerate. Cases the C refuses stay in, so the
  refusal paths are covered too.
- Laggards worth raising: `addfunctions` 57%, `bin/tbfast` 56%, `pairalign`
  59%, `dalign` 75%, `constants` 80%.
- Correctness is not finished while `-@` is refused. The `Gaplen` layer of
  `Dalignmm.c` is ported and checked bit-for-bit; the ~640-line DP is not.

Coverage: `RUSTFLAGS=-C instrument-coverage cargo test`, then `llvm-profdata
merge` + `llvm-cov report` including both binaries as `--object` (the golden
tests drive them as subprocesses; omitting them under-reports by ~5 points).

## 2. A profiling database

Build and keep `bench/profile.tsv`, refreshed by `bench/profile.sh`:

- One row per (input, argument vector, binary), C and port side by side.
- **Compute: callgrind `Ir`.** Deterministic, so a diff is a real change and
  not machine noise — which matters more here than wall-clock realism.
- **Memory: peak RSS** from `/proc/<pid>/VmHWM`, plus massif or DHAT where a
  row needs explaining.
- Per-function `Ir` for the top ~20 functions, so a regression is attributable
  without re-profiling from scratch.
- Commit it. The database is the evidence for every claim in section 3, and a
  diff in it is reviewable.

`perf`, `hyperfine` and `heaptrack` are not in this environment; valgrind's
callgrind/cachegrind/DHAT/massif are. Keep the port dependency-free — no
criterion.

## 3. Performance

**Definition of done: across the profiling database, the port's compute and
memory equal the C's on average.** Better is welcome and is defined as <= 90%
of the C's compute *or* memory.

Memory is already there (0.64x). Compute is the gap: 1.70x instructions,
1.21x wall clock. Wall clock beats instruction count by a wide margin, which
says the port's extra work is cache-friendly — so read the database before
optimising anything.

Rules while doing it:

- **Do not trade evaluation order for speed.** f64 arithmetic must stay in the
  C's order, including its redundant operations. Reassociating changes results.
- Allocation is the obvious first suspect: the C hoists `static` scratch
  buffers that the port often rebuilds per call. Hoisting is safe — it changes
  no arithmetic.
- Bounds checks are the second. Prefer restructuring to slices and iterators
  over `get_unchecked`; `#![forbid(unsafe_code)]` stays.
- Measure per change, against the database. "Should be faster" is not evidence.

## Method

Transliterate; do not improve. When output diverges, instrument a scratch copy
of `core/`, dump intermediates as f64 bit patterns, dump the same from Rust
behind an env var, and bisect. Every mismatch so far was found that way.
