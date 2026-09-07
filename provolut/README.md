# Provolut's MAFFT work

Why this exists: MAFFT is ~87% of the wall-clock of a GEMME/ESCOTT score in the
suggestion pipeline. On a 550-residue parent with a 4,927-sequence alignment,
scoring spends ~41s of 47s inside 100 parallel MAFFT realignments, and `--auto`
picks L-INS-i ("probably most accurate, very slow"). Profiling put ~68% of that
inside one function, `L__align11`.

## The constraint

**Output must not change.** Everything here is verified byte-identical against
the unmodified upstream build, never argued to be equivalent. Anything that
changes an alignment -- a different strategy, `--thread`, `-ffast-math` -- is
out of scope, however fast.

`--thread` is worth naming explicitly: it is fast (3.5x) and it is *not*
result-preserving. Even `--thread 1` differs from no flag at all, on every
sequence of every subset tested. It is not a race; the threaded path is a
different computation.

## What is safe to reassociate, and why

The DP runs entirely on exact integers. Measured, not assumed: every entry of
the scoring matrix is integral (probe over the live matrix reported
`nonintegral=0`, range [-374, 1069]), and `penalty` / `penalty_ex` are C `int`s
(-1199 / -59). Sequence lengths here are in the hundreds, so accumulated values
stay around 1e6 -- far inside the 2^53 where double addition is exact. So
reordering additions cannot change a result, which is what makes a hoisted
prefix scan and a vector kernel legitimate rather than approximate.

## Verification

- `provolut/run_upstream_tests.sh <prefix> <srcdir>` runs upstream's own
  `test/script` cases per-case instead of exit-on-first. 9/9 on both builds,
  including `--allowshift`, which turns on the `trywarp` path this work leaves
  on the original scalar code.
- 50 subsets sampled from a real 4,927-sequence alignment, sizes 8-120:
  50/50 byte-identical.
- End-to-end ESCOTT scores: within the pipeline's own run-to-run variance. Note
  that variance exists -- two runs of the *unmodified* binary differ by ~5e-7,
  which is larger than the difference this work introduces.

## Results

| build | 40-subset corpus | vs conda |
|---|---|---|
| conda-forge mafft 7.526 (what production runs) | 59.4s (10-subset) | 1.00x |
| this source, gcc -O2 | 48.6s (10-subset) | 1.22x |
| + AVX-512 kernel in `L__align11` | 104.6s (40-subset, vs 152.5s stock) | **1.45x** |

Real scorer path, 550aa parent, cold cache: 39.0s -> 25.4s.

Things that did *not* help, so nobody repeats them: `-O3 -march=native`
(no change -- the DP is a serial recurrence with nothing to auto-vectorise),
clang 18 (9% slower than gcc), and branchless-scalar (6% slower; branches were
already predicted at 3.5% miss rate). The win came only once the loop was both
straight-line *and* explicitly vectorised.

## A Rust port of the progressive aligner

`rust/` holds `disttbfast-rs`, a transliteration of MAFFT's progressive MSA
algorithm into Rust, written for bit-exactness rather than speed. It is a
drop-in replacement for `libexec/mafft/disttbfast`, covering FASTA I/O, the
BLOSUM62 score matrix, the 6-mer distance matrix, the UPGMA guide tree and
weights, the group-to-group DP (`A__align`) with its profile-history cache,
and the FFT anchoring (`Falign`, `fft`, `blockAlign2`) that MAFFT is named
for. That is the whole of FFT-NS-1/2 and NW-NS-1/2; the `-i` variants come out
identical too, since `dvtditr` consumes this binary's output unchanged.

The same rule as the rest of this directory applies: output must not change,
and it is verified, never argued. `rust/verify/verify.sh` diffs it against the C
binary over a generated corpus of 74 inputs x 15 argument sets, all
byte-identical; `rust/verify/verify_e2e.sh` does the same through the real `mafft`
wrapper; and upstream's own `test/script` passes 9/9 with the Rust binary
installed. Configurations outside the ported path are refused rather than
approximated. See `rust/README.md`.

Not covered, and the reason this is not yet the binary the scorer would use:
G-INS-i and L-INS-i do not go through `disttbfast` at all. They run `tbfast`
plus `pairlocalalign`, which is where `L__align11` -- the 68% hotspot this
directory is about -- lives. Porting those is the next step.
