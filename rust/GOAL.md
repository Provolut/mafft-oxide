# Goal: prove correctness cheaply, then match the C's speed

**The constraint does not move: output must not change.** Every step below
keeps `tests/matrix.rs`, `verify.sh`, `verify_tbfast.sh` and `verify_e2e.sh`
green. An optimisation that changes one alignment is a bug, however fast.

## Where this starts

`bench/profile.tsv`, 8 cases, C and port side by side. A single case is not
the picture: `disttbfast -F` on `test/sample` alone reads 1.70x / 0.63x /
1.21x, which is the flattering end of the range.

| metric | aggregate | worst case |
|---|---|---|
| instructions (callgrind Ir) | **2.15x** work-weighted | memsave 2.46x |
| wall clock | **1.36x** work-weighted | memsave 1.44x |
| peak RSS | **1.25x** mean per run | linsi 2.67x |

`cargo test`: 4.1 s, 553 bit-identity cases, 89.6% regions / 89.1% lines.

## 1. Testing — a floor to hold, not a milestone to reach

- **>= 85% region coverage from `cargo test` alone.** At 89.6%; do not regress.
  Measure with `bench/coverage.sh`, which gets the two easy mistakes right.
- **Suite stays far under 5 min.** At 4.1 s. Keep it seconds: tiny inputs, no
  sleeps, no network. A slow suite stops being run.
- **Bit-identity is what the tests are for.** `tests/matrix.rs` replays cases
  whose expected bytes came from the real C via `gen_golden.py`; a fixture that
  changes is the port changing its answer. Add a feature -> add its argument
  vector to `gen_golden.py` and regenerate. Cases the C refuses stay in, so the
  refusal paths are covered too.
- Laggards worth raising: `dalign` 75% (capped until the `-@` DP lands),
  `constants` 80%, `msalign` 85%.
- Correctness is not finished while `-@` is refused. The `Gaplen` layer of
  `Dalignmm.c` is ported and checked bit-for-bit; the ~640-line DP is not.

`bench/coverage.sh` encodes the two mistakes that under-report badly: omitting
the binaries from `llvm-cov`'s `--object` list, and a `cargo` call missing the
same `RUSTFLAGS`.

## 2. A profiling database

Build and keep `bench/profile.tsv`, refreshed by `bench/profile.sh`:

- One row per (input, argument vector, binary), C and port side by side.
- **Compute: callgrind `Ir`.** Deterministic, so a diff is a real change and
  not machine noise — which matters more here than wall-clock realism.
- **Memory: peak RSS** from `/proc/<pid>/VmHWM`, plus massif or DHAT where a
  row needs explaining.
- `bench/hotspots.sh` gives the per-function `Ir` breakdown for one case.
- A row is recorded only if both binaries exit 0 *and* produce the same bytes.
  Without that check a crashing port reads as a 1000x speedup, which is how a
  profiling database becomes a lie -- it happened once already.

`perf`, `hyperfine` and `heaptrack` are not in this environment; valgrind's
callgrind/cachegrind/DHAT/massif are. Keep the port dependency-free — no
criterion.

## 3. Performance

**Definition of done: across the profiling database, the port's compute and
memory equal the C's on average.** Better is welcome and is defined as <= 90%
of the C's compute *or* memory.

"On average" means the **work-weighted** ratio -- the sum of the port's cost
over the sum of the C's -- not the unweighted mean of per-case ratios. The two
disagree by 60% here (2.15x against 1.33x) because the three smallest inputs
run *faster* than the C, startup dominating them, which drags an unweighted
mean down and hides that the two heaviest cases are the two worst.
`bench/summarize.py` prints both, weighted first. Peak RSS is the exception:
peaks do not sum, so memory is judged on the per-run mean and the worst case.

Nothing is close to done. Instructions are 2.15x and wall clock 1.36x; memory
averages 1.25x with two cases above 2x. Wall clock is consistently better than
instruction count, which says the extra instructions are cache-friendly -- so
read the database before optimising, and do not assume the two will reach
parity together.

Where the cost is concentrated, from the database:

- **`sample.memsave` 2.46x instructions** -- the linear-space aligner, and the
  most recently written code.
- **`sample.linsi` 2.67x peak RSS** and **`sample.nwns2` 2.19x** -- the two
  places memory goes the wrong way, against 0.54-0.63x on the FFT path.

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
of `core/` and bisect on f64 bit patterns. Every mismatch so far was found
that way, never by reading.
