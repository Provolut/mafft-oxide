# Goal: prove correctness cheaply, then match the C's speed

**The constraint does not move: output must not change.** Every step below
keeps `tests/matrix.rs`, `verify.sh`, `verify_tbfast.sh` and `verify_e2e.sh`
green. An optimisation that changes one alignment is a bug, however fast.

## Where this stands

`bench/profile.tsv`, 8 cases, C and port side by side. One case is not the
picture: `disttbfast -F` on `test/sample` alone reads the flattering end of
every range.

| metric | work-weighted | mean per case | worst | vs done |
|---|---|---|---|---|
| wall clock | **1.02x** | 0.81x | linsi 1.06x | at parity |
| peak RSS | n/a | **0.86x** | linsi 1.36x | **stretch met** |
| instructions (Ir) | 1.38x | 1.00x | nwns2 1.60x | not met |

`cargo test`: 4.1 s, 553 bit-identity cases, 89.6% regions / 89.1% lines.

## 1. Testing — a floor to hold, not a milestone to reach

The floor is enforced, not just written down: `.githooks/pre-push` refuses a
push whose tests fail or whose coverage is under 85%, and refreshes the
sampling profile on the way through. Install with
`git config core.hooksPath .githooks`; about 70 s, most of it the coverage run.

- **>= 85% region coverage from `cargo test` alone.** At 89.8%; do not regress.
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
- `bench/hotspots.sh` gives the per-function `Ir` breakdown for one case,
  under callgrind -- exact, and about 50x, so it is a tool for deciding
  whether a change helped rather than something to run routinely.
- `bench/sample.sh` is the cheap counterpart: gprofng sampling, measured at
  3-4% here against an 8% budget it enforces per case. Over budget it discards
  the profile rather than recording numbers taken under conditions that no
  longer resemble a normal run. The pre-push hook runs it on every push.
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

**Memory is past the stretch bar** at 0.86x mean per run.

**Wall clock is at parity**: 1.02x work-weighted, 0.81x as a mean -- the
typical case is already faster than the C, and `sample.fftns2` runs at 0.98x.

**Instructions are not**: 1.38x work-weighted, though the mean is exactly
1.00x. They have been the worse of the two compute metrics at every single
step, which says the port's extra instructions are cache-friendly. Do not
chase Ir at the cost of wall clock; wall clock is what a caller experiences.

Five changes, each verified byte-identical before being measured:

| | Ir | wall | RSS |
|---|---|---|---|
| start | 2.15x | 1.36x | 1.25x |
| release profile matched to the C's -O3 | 2.03x | 1.27x | 0.92x |
| `L__align11` warp check hoisted out of the row | 1.89x | 1.22x | 0.89x |
| `falign` profile buffers hoisted out of the column loop | 1.40x | 1.05x | 0.85x |
| `match_calc` loop transposed, row slices hoisted | **1.38x** | **1.02x** | **0.86x** |

Every one was a transliteration miss, not a cost of the language: a build
profile that did not match the C's, a check inside a loop the C splits in two,
buffers allocated per column where the C re-zeroes one, and a matrix walked
down its columns where the C's loop order walks along rows.

Two method notes worth keeping:

- **Measure one change at a time against the database.** Of the two
  `match_calc` fixes, the one that looked obvious (hoisting the row slices)
  was worth 1%, and the one that looked like a rounding detail (transposing
  the loops) was worth 4%.
- **Instrumentation can corrupt the metric.** Putting line tables on the
  release profile added 750 KB to the mapped binary and read as a memory
  regression, 0.85x -> 0.93x, entirely on the small cases. Only the per-case
  breakdown showed it; the aggregate just looked like a loss. Line tables live
  in the `profiling` profile now.

What is left, in the order the database ranks it:

- **`sample.nwns2`, 1.60x Ir and 1.17x wall** -- the worst case on both.
- **`sample.linsi`, 1.36x peak RSS and 1.06x wall** -- worst on memory.
- Allocation churn is still worth checking: 625K against the C's 407K.

## Method

Transliterate; do not improve. When output diverges, instrument a scratch copy
of `core/` and bisect on f64 bit patterns. Every mismatch so far was found
that way, never by reading.
