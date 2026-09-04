# Optimisation attempts, 2026-09-04

The run GOAL.md describes: beat the port as it stood at `77dd5a6` by 10% on
**both** wall clock and peak RSS, or spend 90 minutes and write down what was
tried. This is that write-up. The target was not met; what follows is where
the time goes and which stones are now turned.

Output did not change at any point. Every kept change was measured with
`bench/profile.sh`, which refuses to record a row unless both binaries exit 0
and produce identical bytes.

## Result

Not met. Both sides re-measured with the corrected harness (see *Measurement*),
`77dd5a6` rebuilt from a worktree so the two are measured identically:

| metric | baseline `77dd5a6` | now | target (<= 90%) |
|---|---|---|---|
| wall clock, 8 cases | 2372 ms | **2361 ms** (99.5%) | <= 2135 ms |
| mean peak RSS | 6084 KB | **6144 KB** (101.0%) | <= 5476 KB |

Per case, port against port:

```
                    instructions                       peak RSS        wall ms
sample.fftns2   1107938947 -> 1103611601  -0.4%      8864 ->  8940     112 -> 110
sample.nwns2    1515694127 -> 1515693464  -0.0%      8436 ->  8440     141 -> 138
sample.memsave  8950539338 -> 8283448864  -7.5%      5036 ->  5104    1096 -> 1078
sample.linsi    7606390975 -> 7606562578  +0.0%     13440 -> 13468    1001 -> 1014
aa.fftns2         16584195 ->   16545652  -0.2%      2956 ->  3348       5 ->   5
nuc.fftns2        19029295 ->   19006136  -0.1%      3364 ->  3316       6 ->   5
tandem.fftns2     19390787 ->   19222878  -0.9%      3268 ->  3288       6 ->   6
dup.nwns2         12795777 ->   12794559  -0.0%      3308 ->  3252       5 ->   5
```

One line of that is the whole story: **`sample.memsave` gave up 7.5% of its
instructions and 1.6% of its wall clock.** The work removed was real, and it
was nearly free to begin with. That case is not instruction-bound, and neither
is the total -- the port runs 1.33x the C's instructions at 1.015x its wall
clock.

Two of the eight cases carry the workload: `sample.memsave` (1078 ms) and
`sample.linsi` (1014 ms) are 88% of the total. The other six are 5-138 ms and
dominated by process startup, so nothing done to them can move the number. The
memory figure is `sample.linsi`'s: at 13.4 MB it is more than twice the mean,
and it did not move at all this session.

Against the frozen `bench/baseline-2026-09-04.tsv`, which was recorded with the
old harness, the same numbers read 99.7% and 105.9%. Its RSS column is not
comparable with anything measured since; see *Measurement*.

## What was tried

### Kept

1. **`prf1`/`prf2` hoisted out of the per-position loop** (`falign.rs`). The
   port allocated two `vec![0.0; na]` per position where the C re-zeroes one
   pair of buffers. Replaced by two buffers refilled with `fill(0.0)`.
   Allocations on `sample.memsave`: 625,467 -> 555,921. Instructions on that
   case: 8.95 G -> 8.28 G, a 7.5% cut.

2. **Lag-loop slice vectors hoisted** (`falign.rs`). The `Vec<&[u8]>` for each
   group was collected fresh per lag; now built once and `clear()`/`extend()`ed.
   Part of the same measurement as (1).

3. **Bounds check folded in the 26-element profile scan** (`falign.rs`).
   Slicing to `[..26]` up front lets the check fold. Worth **0%** wall clock,
   even though that scan is 21% of the case's instructions -- the checks were
   predicted perfectly. Kept only because it costs nothing.

4. **`FalignState`'s three `MAXSEG` arrays grow on demand.** They were
   allocated at 100000 entries x 24 B = 2.4 MB each regardless of input.
   First recorded as "no change", then as a 2.4 MB saving on `sample.linsi`.
   Both readings were artefacts of the measurement below. In a clean
   directory `sample.linsi` is unchanged by it to within 0.1% -- which is what
   should have been expected, since these arrays belong to `Falign`, the FFT
   path `disttbfast` uses, and the L-INS-i case never touches them. The
   allocation itself is still worth not making.

### Reverted or rejected

5. **`impmtx` walked as a rectangle** rather than per row. Measurably worse;
   reverted.

6. **`-C target-cpu=native`.** AVX-512 is emitted and the output stays
   byte-identical, and it is worth **0%** wall clock. That is not a surprise
   in hindsight: the DP inner loops carry an argmax with an index, which
   blocks vectorisation. `core/Lalign11.c` says exactly this in a comment
   above its own hand-written AVX-512 path -- and that path is compiled out in
   this build anyway, because `CFLAGS = -O3` never defines `__AVX512F__`. Both
   sides are running the scalar loop.

## Measurement, which was the real finding

Four conclusions in this session were wrong the first time, all from the same
metric: **peak RSS as `profile.sh` measured it was not measuring the program.**

The first symptom was spread. The same binary on the same input varied by 31%
between runs; `sample.linsi` clustered at 11.0, 13.4 and 14.4 MB while the C
sat inside 1%. On single samples that produced a change which only removes
allocations reading as a 7.6% memory *regression*, the same change reading as
"no difference", and a median compared against a single draw reading as 5.5%
worse. Sampling five times and taking the median fixed the spread but not the
cause.

The cause is the harness. Every measured run shared one working directory, and
`tbfast` and `disttbfast` write `hat2`, `hat3`, `order` and `pre` into the cwd
and read them back. So the second run of a case was not doing the work the
first one did, and each case inherited the previous case's files. Run in a
fresh empty directory, eleven repeats of `sample.linsi` land inside 1%:

```
port  14332 14356 14420 14420 14428 14432 14436 14448 14452 14468 14472
C      9160  9164  9176  9176  9176  9196  9216  9224  9236  9248  9340
```

`profile.sh` now gives every run -- callgrind, RSS poll and each timing
iteration -- its own directory. The wall-clock column was measured the same
polluted way, so it is re-measured too, on both sides.

`summarize.py` had a matching hole: a baseline whose rows all failed to
measure divided by zero and printed **TARGET MET**. It fired for real, on a
re-measurement pointed at a tree with no compiled C. It now refuses to compare
when either side has no rows marked `ok`, or when the two sides cover a
different number of cases.

The frozen `bench/baseline-2026-09-04.tsv` predates all of this. **It has been
left in place** rather than rewritten: re-measuring it raises the RSS target,
which is a change in my own favour, and that is not mine to make. The
re-measured file sits beside it and the comparison is reported against both.
Which one is the yardstick is a call for a human.

## One more thing the gates caught

The push was refused by `bench/coverage.sh`, which runs the suite a second
time. `blosum`, `io` and `tree` each had a test helper that swaps the
process-wide working directory under a `static LOCK` of its own -- four locks,
one directory, several test threads. About one run in eight, two helpers
interleaved and a test ran in the other's directory. Fixed with a single
`crate::cwd_lock()` taken by all four; no assertion changed. 20 consecutive
instrumented runs clean. Coverage 89.86% regions, 89.36% lines.

Worth recording here because the flake was invisible to the fast gate for as
long as only one run happened per push.

## Leads not taken

- **Three `commonIP` buffers against the C's one.** `align.rs`, `genalign.rs`
  and `lalign.rs` each keep their own traceback matrix with its own
  `common_alloc1`/`common_alloc2` growth. In the C, `commonIP` is a single TLS
  global shared by every aligner (`dp.h:15`, written from `Lalign11.c`,
  `Salignmm.c`, `Dalignmm.c`, `MSalign11.c`, `Galign11.c` and others). At
  `sample.linsi`'s sizes one such matrix is about 2.4 MB, and the L-INS-i run
  has a pairwise stage and a progressive stage that use different ones.
  The buffer is pure scratch -- every cell a traceback reads was written by
  the same call -- so sharing it should be output-neutral, and it is what the
  C does, so it is more faithful rather than less. **Not measured**: it needs
  the buffer threaded through the three states' call sites, which is more
  surgery than the remaining time allowed.

- **The argmax inside the DP row.** Both the port and the C's scalar tail
  update `maxwm`/`endali`/`endalj` inside the inner loop, and that is what
  stops either from vectorising. The C's AVX-512 path lifts it into a
  per-row reduction; doing the same in safe Rust means either a second pass
  over stored pre-clamp scores (extra stores, extra traffic) or a
  hand-blocked reduction. Worth trying, and untried here.

- **`sample.nwns2` at 1.60x the C's instructions** is the worst compute ratio
  left, but only 141 ms absolute; even halving it moves the total by 3%.

## What the ratios say

Instructions overstate the gap consistently. On `sample.memsave` the port
executes 1.23x the C's instructions and runs at 0.98x its wall clock; on
`sample.linsi`, 1.42x instructions for 1.07x wall clock. The port's extra work
is cache-friendly and pipelines behind the loads the C is also waiting on.
Optimising against `Ir` alone would therefore chase work that is already free.
