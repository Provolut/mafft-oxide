# Goal: beat today's port by 10%, or spend 90 minutes proving you cannot

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

| metric | baseline | target (<= 90%) | today |
|---|---|---|---|
| wall clock, 8 cases | 2367 ms | <= 2130 ms | **1249 ms, MET** |
| mean peak RSS | 6725 KB | **<= 6052 KB** | 6725 KB |

Check it with
`./bench/summarize.py bench/profile.tsv bench/baseline-rebased.tsv`, which
exits non-zero until both hold.

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

1. **Hit both targets.** `bench/profile.sh` reports wall clock <= 2130 ms
   *and* mean peak RSS <= 6052 KB, on a run where every row says `ok`, with
   the suite green. Check it with
   `./bench/summarize.py bench/profile.tsv bench/baseline-rebased.tsv`,
   which exits non-zero until both hold. Wall clock is already there at
   1249 ms; the RSS half is the live one.
2. **90 minutes of genuine attempts without hitting it.** Then stop, restore
   the last good state, and write up what was tried and what each attempt
   measured. A negative result that says *where the remaining time goes and
   why it did not move* is the deliverable in that case, and it is worth
   having — it tells the next person which stones are already turned.

Track the time. "90 minutes" means effort spent on optimisation attempts, not
wall time including a long harness run. **Take the two clock readings when you
start and when you stop, and put them in the write-up.** The 2026-09-04 run
recorded "90 minutes" as an estimate because nobody took them, and an untimed
budget cannot be audited — which is the one complaint against that write-up
that could not be answered after the fact.

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
