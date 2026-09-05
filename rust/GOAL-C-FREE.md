# Goal: run `mafft --auto` with no C process at all

**The constraint does not move: output must not change.** `tests/matrix.rs`,
`tests/utils.rs`, `verify.sh`, `verify_tbfast.sh`, `verify_utils.sh` and
`verify_e2e.sh` stay green. A change that alters one alignment is a bug however
much C it removes. `../CLAUDE.md` owns the gates an agent may not lower; the
85% coverage floor and the test suite are not negotiable here either.

This goal sits **beside** `GOAL.md`, it does not replace it. That one is a
performance target against a frozen baseline and its numbers are not this
goal's to move. Where they conflict — and they can, because a port of a hot
binary may land slower than the C before it lands faster — `GOAL.md`'s
targets are advisory for the duration of this work and binding again after.

## What is already C-free, so that the goal is not overstated

The `rust/` crate has no `build.rs` and no dependencies; `tests/matrix.rs` and
`tests/utils.rs` run the Rust binaries against committed fixtures and only
*mention* `core/` in regeneration comments; `bench/coverage.sh` needs no C; and
`.githooks/pre-push` is exactly those two things. **`cargo build && cargo test`
already passes with `core/` deleted.** The library is not the problem.

What is not C-free is the **runtime**: the `mafft` script dispatches to
binaries, and most of them are still the C's.

## The target

`--auto` picks one of seven brackets (`scripts/mafft:1290-1345`), and the
binaries a run executes follow from that choice. Traced by wrapping every
binary in `$MAFFT_BINARIES` with a logger:

| bracket | selects | binaries executed |
|---|---|---|
| `nlen<3000 & nseq<100` | L-INS-i | version, countlen, replaceu, f2cl, **tbfast**, **dvtditr**, restoreu |
| `nlen<1000 & nseq<200` | L-INS-2 | same set |
| `nlen<10000 & nseq<500` | FFT-NS-i | version, countlen, replaceu, f2cl, **disttbfast**, **dndpre**, **dvtditr**, restoreu |
| `nseq<20000` | FFT-NS-2 | as above without `dvtditr` |
| `nseq<100000` | + memsavetree | as above; `disttbfast-rs` already covers `--memsavetree` |
| `nseq<200000` | + memsavetree, one cycle | as above |
| otherwise | parttree | adds **splittbfast** |

Nine distinct binaries. Four are Rust today — `disttbfast`, `tbfast`,
`countlen`, `version`. Five remain: **`replaceu`, `restoreu`, `f2cl`,
`dndpre`, `dvtditr`**, and a sixth, `splittbfast`, only for the parttree
bracket.

**The target is the first six brackets: every input with fewer than 200000
sequences runs `--auto` with no C process.** That is five binaries, and it
covers every workload this repository actually serves, including the
GEMME/ESCOTT trace in provolut-ml, which never draws a sub-alignment above 51
sequences.

The parttree bracket is deliberately out. Say so in the README rather than
letting it look finished: a 200000-sequence input keeps executing the C
`splittbfast`, and it must keep *working*, not start refusing. Silently
falling back is fine here precisely because it is documented; what is not fine
is claiming the run is C-free when it is not.

## Done is all four of these

1. **Coverage.** A `$MAFFT_BINARIES` directory containing only Rust binaries
   runs `mafft --auto` to completion on an input from each of the first six
   brackets. Prove it by *removing* the C ones from that directory, not by
   trusting a symlink: a missing binary must fail the run, which is what makes
   a passing run mean something.
2. **Identity.** Byte-identical stdout against the all-C configuration on: the
   840-payload ProteinGym corpus, upstream's own test suite (9/9), and the
   committed golden matrix. Plus the end-to-end proof that actually matters —
   provolut-ml's ESCOTT scorer producing SHA-256-identical scores and
   prediction matrices on `SPG1_STRSG_Olson_2014` and `P53_HUMAN_Kotler_2018`,
   each run with its own fresh trace cache.
3. **Independence.** `mv core core.hidden && cargo build --release && cargo test`
   passes, and the run in (1) still works. The C becomes a development
   dependency of the differential harnesses, nothing more.
4. **No speed regression, measured per input shape.** On the ESCOTT corpus,
   against the all-C configuration:

   | shape | four-Rust baseline (2026-09-05) | target |
   |---|---|---|
   | 51x393 sub-alignment (where ESCOTT spends its time) | 0.983x | >= 0.98x |
   | 11x448 sub-alignment (startup-dominated) | 0.971x | >= 0.96x |

   Two floors rather than one, because a single number hides which shape moved.
   The short-input figure is **not** the utilities' doing: the same corpus with
   only `tbfast`/`disttbfast` swapped measured 0.966x before `countlen-rs` and
   `version-rs` existed, so the gap is `tbfast-rs` on calls whose DP is too
   small to amortise fixed costs, and adding the two utilities did not widen
   it. They do cost something, and it is small: 0.31 ms against 0.25 ms per
   `version` invocation and 0.35 ms against 0.31 ms per `countlen`, from
   binaries of 380-420 KB against 16-125 KB. About 0.1 ms per `mafft` call,
   roughly 4% of the 2.9 ms gap on that shape.

   `dvtditr` is 22% of an L-INS-i call and is the one that can still move
   these; a port that lands the binary but costs 15% is not done.

## Or: a written negative result

If a piece proves genuinely unportable — a path whose reference output is
undefined, a dependency on C library behaviour with no faithful equivalent —
stop and write it up: what was tried, what it measured, and why the remaining
gap is where it is. `bench/OPTIMISATION-LOG.md` is the format. A negative
result that names the obstacle is worth having; a port that quietly disagrees
with the C is not.

One such case already exists and is the template. `getnumlen_nogap_countn`
reads an uninitialised `nN` for an input whose first record has an empty body
(upstream issue #8, valgrind-confirmed), so **there is no bit-exact answer to
reproduce**. `countlen-rs` picks the defined value, zero, which agrees with the
C on every input where the C itself is defined, and `verify_utils.sh` scores
that input `indeterminate` rather than as a pass or a failure. Expect more of
these in `dvtditr`, and handle them the same way: name the input, prove the
indeterminacy, choose deliberately, and do not let it hide inside a pass count.

## Order, and why

1. **`replaceu`, `restoreu`, `f2cl`** — small, and every bracket runs them, so
   they are the cheapest reduction in C-process count. The reader they need
   already exists: `load1SeqWithoutName_realloc_casepreserve` (`io.c:1387`)
   applies `charfilter()`, which `io.rs` has as `SeqFilter::Text`. Note
   `getnumlen_casepreserve` (`io.c:2132`) caps sampling at **1000000**, not
   100000, and calls `countATGC` rather than `countATGCandN` — so issue #8's
   defect does not exist on that path. Do not conflate the two.
2. **`dndpre`** — needed by the FFT-NS-i and FFT-NS-2 brackets, and small.
3. **`partA__align`** (`partSalignmm.c`) — 40.7% of `dvtditr`'s instructions on
   a 51x393 sub-alignment and 61.7% on an 11x448 one, measured with callgrind
   on both shapes. It is a sibling of the already-ported `A__align` in
   `align.rs`; check first whether it carries a `cpmxhist`-style history cache,
   because in `align.rs` that cache changed the summation order and was
   load-bearing for bit-exactness.
### What reading `partSalignmm.c` changed about (3)

Measured against the file rather than its line count, and each of these is a
saving someone would otherwise spend a day rediscovering:

- **`partA__align_variousdist` (746 lines) is never reached.** Refuse it. It is
  the single biggest item in the file and none of it is needed.
- **`st_OpeningGapCount` / `st_FinalGapCount` are live, and already ported.**
  An earlier revision of this list called them dead, because the copies at the
  top of `partSalignmm.c` sit inside `#if 0`. That was wrong: those are dead
  *local duplicates*, and `partA__align`'s calls bind to the real definitions
  in `mltaln9.c:12837`, declared in `functions.h`. Either way there is nothing
  to write — `profile.rs` already has `st_opening_gap_count` and
  `st_final_gap_count`.
- **`match_calc` is the function `align.rs` already has** — same signature, same
  sparse-profile walk, same register accumulation (this repo's own patch). The
  only difference is that `partA__align` reads the global
  `n_dis_consweight_multi` where `A__align` takes `n_dynamicmtx`; pass
  `cn.n_dis_consweight_multi` and it is exact. It is now `pub(crate)` for that
  reason. **Do not write a second copy** — two transliterations of one DP inner
  loop are free to drift apart, and nothing would catch it but a differential
  run.
- **`partA__align` carries no `cpmxhist` cache** (`grep -c` is 0). None of
  `A__align`'s summation-order subtlety applies here.
- **The constraint half is done** (`partalign.rs`): `part_imp_match_init_strict`
  is allocation plus the already-ported `fillimp`, and the two gapmap
  accumulators are five lines each.
- **The two tracebacks are NOT reusable — checked, and the guess was wrong.**
  They open the same way (`partSalignmm`'s `Atracking` starts with
  `atracking_localhom`'s `outgap == 0` scan, not `atracking`'s `tailgp` form),
  which is what made reuse look plausible. Diffing the two C
  `Atracking_localhom` bodies with comments and whitespace stripped gives 148
  differing lines out of ~185 and ~140:

  * `Salignmm.c`'s builds **gap tables** (`gt1bk` / `gt2bk`, `gaptable1/2`) and
    takes `ngap1`, `ngap2`, `reuseprofiles`, `gt1`, `gt2` — the profile-reuse
    machinery that belongs with `cpmxhist`.
  * `partSalignmm.c`'s writes **directly into `mseq1` / `mseq2`** and takes
    `start1/end1/start2/end2` plus the gap maps.

  Port them separately. The consolation is that `partSalignmm`'s pair is the
  *older, simpler* shape — no gap-table bookkeeping at all — so they are easier
  than the ones already in `align.rs`, not harder. About 290 lines for the two.

  Recording the negative result on purpose: "these two look alike" is exactly
  the assumption that would have been discovered halfway through a port, by
  which point the wrong shape is already written.

So what remains is the 696-line DP, the two tracebacks (~290 lines, ported
fresh), and the `tditeration` driver plus `dvtditr`'s arguments — about 1000
lines, not the 2125 the file suggests.

**Start it with a full context budget.** It is one long transliteration where a
single reordered addition changes an alignment, no test says which line did it,
and the method that works here is instrumenting a scratch `core/` and bisecting
on f64 bit patterns. That wants the whole budget for the bisect.

### What `dvtditr` actually costs, measured per function

This estimate has been corrected three times in this document, each time
upward, because it was made from file sizes and skim-reading. Measured per
function it is:

| lines | function | needed on the `-C 0` path? |
|---|---|---|
| 1188 | `TreeDependentIteration` (`tditeration.c`) | yes — the sequential refinement loop is the bulk of it |
| 925 | `athread` (`tditeration.c`) | **no** — the threaded worker; `-C 0` takes the sequential branch |
| 66 | `classifypairs` | probably; check |
| 81 | `freelocalarrays`, 81 `Writeoptions`, 20 `shuffle`, 13 `makescoringmatrices` | small |
| ~500 | `dvtditr.c`'s `main` reached path | `hat3` table reading, `topol`/`len`/`eff`, subalignments, `_focus_` targets |

So the driver is roughly **1400 reached lines**, on top of the ~360-line DP
already written. Earlier notes here said "about 1000 lines" for everything
remaining; that was wrong, and the difference is not a rounding error.

The saving is real but smaller than it looks: `athread` is skipped only because
this port is single-threaded, and a future `-C n` would need it or an
equivalent.

### Status: the DP is written and **not verified**

`partalign::PartAlignState::part_a_align` exists and compiles, and the suite is
green — which says nothing about whether it is correct. Nothing has executed
it: `dvtditr-rs` does not exist, so no differential run has reached it. A DP
that compiles has not been tested in any sense that matters here.

What will decide it, in order:

1. `dvtditr-rs` calling it, so the code runs at all.
2. `verify_tbfast.sh` and `verify_e2e.sh` against a built `core/`.
3. The 840-payload ProteinGym corpus through `mafft --auto`, byte-compared.
4. ESCOTT end to end on the two assays, SHA-256 compared.

Until (2) passes, treat the DP as a draft. The two `>` comparisons below are
the most likely place for it to be wrong, and they are also the least likely to
announce themselves.

### The driver is bigger than an earlier note here claimed

That note said `dvtditr.c`'s 1214 lines "are mostly argv that `params.rs`
already models". That is wrong, and the correction matters for planning: its
`main` also reads the `hat3` localhom table (`readlocalhomtable2_half` /
`_target`), builds `topol` / `len` / `eff`, handles `--subalignment` tables and
`_focus_` targets, and only then calls `tditeration`. Add `tditeration.c`'s own
refinement loop and the remaining work is **larger than the DP**, not a thin
wrapper over it.

Still absent from the Rust side, as far as the reached set goes:
`alignableReagion`, `new_OpeningGapCount`'s tail-gap siblings already noted,
the `tditeration` accept/reject loop, and `dvtditr`'s table reading.

### Read this before writing one line of the DP

**`partA__align` and `A__align` disagree on two tie comparisons, in opposite
directions.** Copying `align.rs`'s inner loop — the obvious move, since it is
the same DP and already ported — gives the wrong alignment on every tie, at
identical score. No unit test and no score check can see it.

| recurrence | `Salignmm.c` / `align.rs` | `partSalignmm.c` |
|---|---|---|
| `*prept + *ogcp2pt * gf1vapre` against `mi` | **`>=`** (the `>` form is commented out) | **`>`** (the `>=` form is commented out, `// 2018/Apr`) |
| `*prept + ogcp1va * *gf2ptpre` against `*mjpt` | **`>=`** (`>` commented out) | **`>`** (`>=` commented out, `// 2018/Apr`) |

Both files carry both forms, one live and one commented, and they are mirror
images: upstream changed these two in `partSalignmm.c` in April 2018 and left
`Salignmm.c` on `>=`. It is deliberate, not a bug, so there is nothing to file
upstream — it is a porting hazard, and the reason this section exists.

The other two comparisons in the same loop (`mi + *fgcp2pt * gf1va` and
`*mjpt + fgcp1va * *gf2pt`, both against `wm`) use `>` in **both** files.

So: write the DP from `partSalignmm.c`, and use `align.rs` only for the shape
of the surrounding Rust — buffer ownership, the `previousw`/`currentw` swap,
how `match_calc` is called. Diff the four comparisons against the C by eye
before running anything.

### `partA__align`'s 696 lines, sliced

Read off its control flow (`core/partSalignmm.c`, line numbers at the time of
writing). Nothing here is ported yet; this is the map, so the next session
starts from structure rather than a wall of code.

| lines | slice | notes |
|---|---|---|
| 753-797 | teardown (`seq1 == NULL`) | the C's free-everything call; an owned struct replaces it, as `PartImp` already did |
| 807-832 | degenerate `lgth1 == 0` / `lgth2 == 0` | fills the other side with `newgapstr`; self-contained, unit-testable |
| 843-862 | `trywarp` setup | **refuse** — this port rejects `trywarp` before reaching here |
| 868-985 | scratch growth (`orlgth1/2`, `commonAlloc1/2`) | becomes owned buffers; the C's `static TLS` is not a contract |
| 991-1029 | gap frequencies from `sgap1`/`egap1`, and `legacygapcost` | two branches, both short; testable on their own |
| 986-1060 | profile and gap-frequency setup | **entirely wiring** — see the table below; nothing here needs writing |
| 1064-1097 | boundary row and column, under `outgap` | |
| 1124-1315 | **the main DP loop** | the real work; the `trywarp` blocks inside (1190-1197, 1249-1290, 1309-1326) are refused |
| 1329-1335 | `!outgap` last row/column fixup | |
| 1340-1378 | traceback and result assembly | tracebacks are **done** (`partalign::atracking*`) |

**The whole setup block is already ported.** Every function
`partA__align` calls between lines 986 and 1060 exists in `profile.rs`:

| C | Rust |
|---|---|
| `cpmx_calc_new` | `profile::cpmx_calc_new` |
| `st_OpeningGapCount` / `st_FinalGapCount` | `profile::st_opening_gap_count` / `st_final_gap_count` |
| `new_OpeningGapCount` / `new_FinalGapCount` | `profile::new_opening_gap_count` / `new_final_gap_count` |
| `outgapcount` | `profile::outgapcount` |
| `gapcountf` | `profile::gapcountf` |
| `match_calc` | `align::match_calc` (now `pub(crate)`) |
| `fillimp` | `localhom::fillimp`, behind `partalign::PartImp` |

So that slice is calls and arithmetic — the `0.5 * (1.0 - ogcp) * fpenalty *
gapfreq` scaling, the `legacygapcost` branch, the `1.0 -` inversions — not new
functions.

Two observations that shape the order:

- **Only the main loop resists isolated testing.** Every other slice can be
  unit-tested against hand-built inputs, which is how the tracebacks were done.
  The loop only becomes verifiable when the whole function runs against the C.
- **So do not ship a partial DP.** Land the testable slices if you like, but the
  loop goes in one piece, and the first thing that should pass afterwards is
  `verify_tbfast.sh` / `verify_e2e.sh` against a built `core/` — not a unit
  test, which cannot see a reordered addition.

The four `if( chudanpt && *chudanpt != chudanref )` checks (1129, 1203) are the
C's cooperative-cancellation hooks; nothing in this repo's call path sets
`chudanpt`, so they are dead here.

4. **`tditeration`'s reached set + the `dvtditr` driver.** `movereg`,
   `movereg_swap`, `fillimp`, `intergroup_score`, `cpmx_calc_new`, `gapcountf`
   and `Atracking_localhom` are **already ported** (`localhom.rs`, `score.rs`,
   `profile.rs`, `align.rs`) — the genuinely missing pieces are
   `alignableReagion`, `new_OpeningGapCount`, `new_FinalGapCount`,
   `commongappick_record` and the refinement driver. `dvtditr.c`'s 1214 lines
   are mostly argv that `params.rs` already models.

`Lalignmm.c` (2563 lines), `Dalignmm.c` (6266) and `rna.c` (527) are never
reached on any `--auto` path. Refuse them loudly rather than port them; that
choice is what keeps this at roughly 2-3k lines instead of 11k.

## Rules

- **Measure the reached set, do not infer it.** Every scope claim above came
  from wrapping the binaries and from callgrind, not from reading Makefiles.
  A new binary gets the same treatment before anyone estimates it.
- **Never reassociate f64 arithmetic.** Unchanged from `GOAL.md`. Loop
  transposition that preserves each accumulator's summation order is fine.
- **Transliterate; do not improve.** The C's quirks are the contract. Two are
  already reproduced deliberately — `version`'s missing trailing newline, and
  `countlen`'s N tally carrying past the sampling cap. Where a quirk is
  *undefined* rather than merely odd, see "a written negative result" above.
- **Add to the harnesses, never subtract.** `verify_utils.sh` is the pattern:
  a wide differential against a live C for anyone extending coverage, and
  committed fixtures so the suite still runs with no C at all. Freeze new
  corpora as fixtures; do not narrow an existing one.
- **Every new binary goes into `bench/coverage.sh`'s `--object` list.** Leaving
  one out reports it at 0.00% while its tests are in fact driving it — the same
  under-report that script's header warns about, pointed the other way. It cost
  a re-run here already.
- **Keep `#![forbid(unsafe_code)]` and keep the port dependency-free.**

## Method

Wrap, trace, then port. `MAFFT_BINARIES` pointed at a directory of logging
wrappers gives the exact binary set and the exact argument vector a real run
uses; replaying that vector against a staged working directory is how a single
binary gets benchmarked and diffed in isolation. When output diverges,
instrument a scratch copy of `core/` and bisect on f64 bit patterns — every
mismatch in this port so far was found that way, never by reading.
