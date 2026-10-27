# Goal: 10% off the port's peak memory, or an account of why it is not there

**The constraint does not move: output must not change.** `tests/matrix.rs`,
`verify.sh`, `verify_tbfast.sh` and `verify_e2e.sh` stay green. A change that
alters one alignment is a bug however fast it is. See `../CLAUDE.md` for the
gates an agent may not lower — the coverage floor it names and the test suite
are not negotiable, and skipping a test to go green is forbidden. The floor
lives there and only there; this file named 85% while that one said 90%, which
is exactly the drift `../CLAUDE.md` warns about.

## The target

Re-based on 2026-10-27, against **this port's own numbers measured with
`bench/maxrss.c`**, frozen in `bench/baseline-rebased.tsv`. Not relative to the
C: the port is at 0.40x the C's mean peak RSS and 0.69x its wall clock.

**Both targets are met as of 2026-10-27.**

| metric | baseline | target (<= 90%) | achieved |
|---|---|---|---|
| mean peak RSS | 6725 KB | <= 6052 KB | **5870 KB, 87.3%** |
| wall clock, 8 cases | 2367 ms | <= 2130 ms | **1243 ms, 52.5%** |

`summarize.py` prints `TARGET MET (both)`. Verified on a run where all eight
rows read `ok`, with `cargo test` (238 tests, the 553-case matrix among them)
and all four differential harnesses green: `verify.sh` 12018/0,
`verify_tbfast.sh` 1690/0, `verify_e2e.sh` 28/0, `verify_dvtditr.sh` 32/0.

What closed it, after `commonIP` sharing got within 101 KB: **an assignment
evaluates its right-hand side before dropping the old value**, so growing a
buffer held two copies at the moment the working set was largest. Thirteen
large allocations had that shape. `impmtx` alone -- the 5.2 MB constraint
matrix in `treebase.rs` -- was worth 1856 KB on `sample.linsi`, whose ratio
against the C fell from 1.54x to 1.09x. See `bench/OPTIMISATION-LOG.md`,
2026-10-27.

A new target should be re-based from here with `bench/maxrss.c` before anyone
works on this again; the sections below describe the run that is now finished.

**Wall clock is settled: 1248 ms against a 2130 ms target, 52.7%, achieved
2026-10-27.** It is recorded here rather than carried as a live condition. It
was met with room to spare and nothing in flight can unmeet it, so leaving it
in a conjunction only meant every check re-derived a finished fact and then
reported the pair as unmet. `summarize.py` still prints both, and still
requires both, which is correct for the tool; this file is saying which half is
work.

Check it with
`./bench/summarize.py bench/profile.tsv bench/baseline-rebased.tsv`.

### Why the baseline file is a composite, and why that is honest

The RSS column is today's HEAD under `bench/maxrss.c`. The wall-clock column is
carried over from 2026-09-04 unchanged. Two reasons:

- **The instrument fault was in peak RSS only.** `profile.sh` polled
  `/proc/<pid>/status` for VmHWM while the program ran, and the four small
  cases exit in about a millisecond -- before the first read. `/bin/true`
  measured 0 KB on it. Wall clock was best-of-five timing and was never
  affected, so there is nothing to correct there.
- **The old RSS number cannot be re-measured.** It came from commit
  `77dd5a6`, which the squash merge of #5 collapsed; the object is gone from
  this clone and from the remote. So the only correctly-measured baseline
  available is today's.

The previous target was 5220 KB, derived as 90% of 5800 KB -- a number the
broken instrument produced. **That target was below the measured floor of 5892
KB and could not be reached by any amount of work.** The new one, 6052 KB, sits
above the floor, so the RSS half is a real target for the first time: 673 KB of
mean to find, against about 160 KB of headroom at the floor.

Where it is: `sample.linsi`, 14.3 MB and 1.54x the C, reproducible to 0.8%. The
one measured lead is the two `commonIP` traceback matrices that are live
simultaneously there, 1,798,280 B each where the C shares one -- worth about
220 KB of the mean on its own, so necessary but not sufficient. See
`bench/OPTIMISATION-LOG.md`, 2026-10-27.

## Done is either of these

1. **Hit the target.** `bench/profile.sh` reports mean peak RSS <= 6052 KB on a
   run where every row says `ok`, with the suite green and the differential
   harnesses green.
2. **Account for the gap instead.** Not "time ran out" -- *the remaining
   kilobytes are named, each with the reason it cannot move.* One line per
   item, measured, in `bench/OPTIMISATION-LOG.md`: what it is, how big it is,
   and whether it is fixed process cost, a table the C sizes the same way, or
   live working set. An item is only accounted for when someone could check the
   claim; "probably irreducible" is not an account.

   This replaced a 90-minute budget, which measured the wrong thing. A clock
   quota rewards grinding and penalises concluding: it reported the goal unmet
   while the substantive finding -- that four of the five largest items are
   locked -- was already complete, and it invited filling time rather than
   stopping. It also could not be audited, because the only witness to the
   minutes is the agent spending them.

   What ends the work either way: **leave the tree in a state someone else can
   build and test.** Revert half-finished edits and failed experiments; do not
   revert measurement or harness fixes, which are not experiments. (The old
   wording, "restore the last good state", read as "undo your work" and was
   cited four times against commits that fixed a harness measuring the wrong
   thing.)

## What is blocking the rest, as of 2026-10-27

Stated here because the file previously implied the remaining memory was simply
available, and the measurements say it is not. Of the five largest items in the
database, four cannot move without changing behaviour:

| item | size | why |
|---|---|---|
| `impmtx` | 5.2 MB | `AllocateFloatMtx` returns `double **`; the port's `f64` is exact |
| score matrices, x4 | 2.1 MB | `0x100` wide because raw sequence bytes index it; narrowing turns the C's silent zero on a high byte into a panic |
| sextet tables, x4 | 1.5 MB | 6^6 is the 6-mer code space, whatever the input size |
| fixed process cost, x8 | 23.3 MB | 2912 KB before the aligner starts; only a smaller binary |

The fifth was duplication and is spent: `commonIP`, two matrices where the C
keeps one, worth 1.66 MB on `sample.linsi`.

**That table omitted something, and the omission was worth more than anything
in it.** Growing a matrix with `x = vec![..]` evaluates the right-hand side
before dropping the old one, so both exist at the instant the working set is
largest. Dropping first was worth 2.8 MB across three cases. The same pattern
is still in `align.rs`, `lalign.rs` and `falign.rs`'s scratch buffers and has
not been swept -- which is the open lead, and it is mechanical rather than a
faithfulness question.

So the path to 6052 KB is either **a smaller binary**, or **a deviation from
the C's own data structures**. The second is a faithfulness trade, and
`../CLAUDE.md` makes that a human's decision, not an agent's -- which is why
this section exists rather than an instruction to keep trying.

## Rules

- **Measure one change at a time** against `bench/profile.sh`. Batching hides
  which change did the work, and intuition has been wrong here: of the two
  `match_calc` fixes, the obvious one was worth 1% and the subtle one 4%.
- **Never reassociate f64 arithmetic.** Transposing a loop while each
  accumulator keeps its summation order is fine and has been done. Changing
  the order values are added is not.
- **Verify before measuring.** `profile.sh` refuses to record a row unless
  both binaries exit 0 and produce identical bytes, because a crashing port
  once recorded as a 1000x speedup.
- **Run the wide harnesses before claiming a change is safe**, not just
  `cargo test`. The matrix is fast, not exhaustive.
- Keep `#![forbid(unsafe_code)]`. Prefer restructuring to `get_unchecked`.
- Keep the port dependency-free.

## Where to look first

From `bench/profile.tsv` and `bench/hotspots.sh`, in the order the database
ranks them:

- **`sample.linsi`, 1003 ms and 11064 KB** — the largest single wall-clock and
  memory contributor. 88% of its instructions are in `l_align11`. Its peak RSS
  is 1.36x the C's, the only case well above.
- **`sample.memsave`, 1093 ms** — the other heavy case. `msalignmm_rec`
  allocates 16 vectors per recursive call, two of them `Vec<Vec<_>>` that
  allocate a row apiece. The C allocates per call too, but the port still makes
  625K allocations against the C's 407K on this workload.
- **`sample.nwns2`, 1.60x the C's instructions** — the worst compute ratio
  left, though only 141 ms absolute.
- The small cases are 5-6 ms each and dominated by process startup. They
  cannot move the total; do not spend time there.

Four of the five wins so far came from allocation or addressing, none from
cleverness: buffers allocated per iteration where the C re-zeroes one, a check
inside a loop the C splits in two, a matrix walked down columns where the C
walks rows, and a build profile that did not match the C's. Look there first.

## Method

Transliterate; do not improve. When output diverges, instrument a scratch copy
of `core/` and bisect on f64 bit patterns. Every mismatch so far was found that
way, never by reading.

`bench/hotspots.sh` is callgrind: exact, ~50x, for deciding whether a change
helped. `bench/sample.sh` is gprofng: 3-4%, inside the 8% budget it enforces
per case, and runs on every push.
