# Goal: prove correctness cheaply, then match the C's speed

**The constraint does not move: output must not change.** Every step below
keeps `tests/matrix.rs`, `verify.sh`, `verify_tbfast.sh` and `verify_e2e.sh`
green. An optimisation that changes one alignment is a bug, however fast.

## Where this stands

`bench/profile.tsv`, 8 cases, C and port side by side. One case is not the
picture: `disttbfast -F` on `test/sample` alone reads the flattering end of
every range.

| metric | aggregate | worst case | vs done |
|---|---|---|---|
| wall clock | **1.05x** work-weighted | nwns2 1.23x | within 5% |
| instructions (callgrind Ir) | **1.40x** work-weighted | nwns2 1.69x | not met |
| peak RSS | **0.85x** mean per run | linsi 1.78x | **stretch met** |

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

**Memory is past the stretch bar**: 0.85x mean per run.

**Wall clock is within 5% of parity**: 1.05x work-weighted, and 0.84x as an
unweighted mean -- most individual cases are already faster than the C.

**Instructions are not**: 1.40x. They have been consistently the worse of the
two metrics, so the port's extra instructions are cache-friendly; do not
assume the two reach parity together, and do not chase Ir at the cost of
wall clock.

Four changes so far, each verified byte-identical before being measured, and
every one of them a transliteration miss rather than a cost of the language:

| | Ir | wall | RSS |
|---|---|---|---|
| start | 2.15x | 1.36x | 1.25x |
| release profile matched to the C's -O3 | 2.03x | 1.27x | 0.92x |
| `L__align11` warp check hoisted out of the row | 1.89x | 1.22x | 0.89x |
| `falign` profile buffers hoisted out of the column loop | **1.40x** | **1.05x** | **0.85x** |

The last was the big one and shows the shape of the remaining work.
`hotspots.sh` had the port at 43% of its instructions inside malloc/free;
DHAT counted 18,810,901 allocations against the C's 406,972, while peak
footprint was already *lower* than the C's -- churn, not size. Two adjacent
sites in `falign::falign` held 18.2M of them: `prf1`/`prf2` allocated per
column, where the C re-zeroes one hoisted buffer. Now 625,467 allocations and
180 MB total against the C's 407,062 and 548 MB.

What is left, in order of what the database says:

- **`sample.nwns2`, 1.69x Ir and 1.23x wall** -- now the worst case on both.
- **`sample.linsi`, 1.78x peak RSS** -- the one case where memory is well
  above the C, against 0.35x on the FFT path.
- Look for more allocation churn first: 625K against the C's 407K says there
  is still some, and it has been the cheapest thing to find every time.

## Method

Transliterate; do not improve. When output diverges, instrument a scratch copy
of `core/` and bisect on f64 bit patterns. Every mismatch so far was found
that way, never by reading.
