# Goal: beat today's port by 10%, or spend 90 minutes proving you cannot

**The constraint does not move: output must not change.** `tests/matrix.rs`,
`verify.sh`, `verify_tbfast.sh` and `verify_e2e.sh` stay green. A change that
alters one alignment is a bug however fast it is. See `../CLAUDE.md` for the
gates an agent may not lower — the 85% coverage floor and the test suite are
not negotiable, and skipping a test to go green is forbidden.

## The target

Relative to **this port as of 2026-09-04**, frozen in
`bench/baseline-2026-09-04.tsv`. Not relative to the C: parity with the C is
already reached and is no longer the objective.

| metric | today | target (<= 90%) |
|---|---|---|
| wall clock, 8 cases | 2367 ms | **<= 2130 ms** |
| mean peak RSS | 5800 KB | **<= 5220 KB** |

**Either one is enough.** Both is better. Instructions (1.378x of the C) are
not a target — they have been the pessimistic metric at every step, and wall
clock is what a caller experiences.

## Done is either of these

1. **Hit the target.** `bench/profile.sh` reports wall clock <= 2130 ms or
   mean peak RSS <= 5220 KB, on a run where every row says `ok`, with the
   suite green.
2. **90 minutes of genuine attempts without hitting it.** Then stop, restore
   the last good state, and write up what was tried and what each attempt
   measured. A negative result that says *where the remaining time goes and
   why it did not move* is the deliverable in that case, and it is worth
   having — it tells the next person which stones are already turned.

Track the time. "90 minutes" means effort spent on optimisation attempts, not
wall time including a long harness run. Record start and finish in the write-up.

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
