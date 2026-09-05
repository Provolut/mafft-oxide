# Goal: run `mafft --auto` with no C process at all

**The constraint does not move: output must not change.** `tests/matrix.rs`,
`tests/utils.rs`, `tests/flh.rs`, `verify.sh`, `verify_tbfast.sh`,
`verify_utils.sh`, `verify_flh.sh` and `verify_e2e.sh` stay green. A change that alters one alignment is a bug however
much C it removes. `../CLAUDE.md` owns the gates an agent may not lower.

This sits **beside** `GOAL.md`, which is a performance target against a frozen
baseline; its numbers are not this goal's to move.

---

## 1. Status

Nine binaries can run under `mafft --auto`. **Eight are Rust and verified
byte-identical.** The ninth, `dvtditr`, is the whole of what remains.

| binary | state |
|---|---|
| `disttbfast`, `tbfast` | Rust, byte-identical (pre-existing) |
| `countlen`, `version`, `replaceu`, `restoreu`, `f2cl`, `dndpre` | Rust, byte-identical |
| **`dvtditr`** | **not done** — its DP and block splitter are, its driver is not; see §5 |
| `splittbfast` | **out of scope** — parttree bracket only, ≥200000 sequences |

Verified with eight Rust binaries in place: the 840-payload ProteinGym corpus
byte-identical through `mafft --auto`, upstream's own suite, and provolut-ml's
ESCOTT scorer SHA-256-identical on two assays. Speed 0.989× (51×393 windows)
and 0.981× (11×448) against the all-C build.

`f2cl` implements the **FASTA path only**; clustal and phylip are refused
rather than approximated. `--clustalout` / `--phylipout` still need the C.

---

## 2. Start here: the environment

None of this is in the repo, and all of it is needed.

```sh
# Rust toolchain and valgrind (bench/profile.sh and bench/hotspots.sh need it).
mamba create -y -n rustdev -c conda-forge rust
mamba install -y -n rustdev -c conda-forge valgrind
export PATH=$HOME/miniforge3/envs/rustdev/bin:$PATH

# The pre-push gate is NOT wired in a fresh clone. Wire it, or you will push
# ungated and not notice — this happened here for thirteen commits.
git config core.hooksPath .githooks

# Building core/ for the differential harnesses: conda's gcc 15 FAILS.
# `dndblast.c` and `filter.c` are built by make's implicit rule without
# $(MYCFLAGS), so they lose -std=c99 and C23 rejects mltaln.h's old-style
# declarations (upstream issue #7). Use the system gcc 13:
cd core && make -j CC=/usr/bin/gcc-13 CFLAGS="-O3"
```

`make` writes the binaries into `../binaries/`. Point a run at a build with
`MAFFT_BINARIES=/path/to/binaries scripts/mafft ...`; that is how a mixed
C/Rust configuration is assembled — copy `binaries/`, drop the `*-rs` binaries
over the C ones under their C names.

---

## 3. How you know you are right

In increasing order of authority. **A unit test cannot see a reordered
addition; only these can.**

1. `cargo test --release` — golden matrix and the utility fixtures. Needs no C.
2. `./verify_utils.sh ../core [corpus-dir]` — the utilities, wide. Pass a
   directory of FASTA files to widen it.
   `./verify_flh.sh ../core` — `Falign_localhom` and, through it,
   `partA__align`, against the real C on 98 staged calls. Needs only a built
   `core/`, runs in seconds, and `--regenerate` refreshes the fixtures
   `cargo test` replays. **This is what proved the DP**, and it is where the
   next slice of `dvtditr` should be proved too: staging one function against
   the C costs an hour and finds what no unit test can.
3. `./verify.sh`, `./verify_tbfast.sh`, `./verify_e2e.sh` — the aligners,
   against a built `core/`. About an hour. **This is what decides `dvtditr`.**
4. `./bench/coverage.sh` — 85% floor.
5. The ProteinGym corpus and ESCOTT, below.

### Rebuilding the ProteinGym corpus

The 840-payload corpus is the realistic workload: the exact FASTA payloads
provolut-ml's GEMME/ESCOTT trace hands MAFFT. It is not in this repo. Rebuild
it from provolut-ml with `pml.bioseq.gemme.pipeline`'s own sampler —
`_aligned_records` → `filter_homologs` → `sample_sub_alignments` at the iGEMME
defaults (`n_iter=7`, `oversample=2`, `max_homologs=40000`), writing each
subset through `_degap` as `>name\nseq`. MSAs live at
`$PROTEINGYM_MSA_DIR` (default `~/gemme_pg/msa_extract/DMS_msa_files`).
`P53_HUMAN_full_04-29-2022_b09.a2m` and `SPG1_STRSG_full_11-26-2021_b07.a2m`
give 100 and 35 payloads at stride 7 / 4.

### The ESCOTT end-to-end check

The strongest signal, and the one that matters to the caller. In provolut-ml:

```sh
GEMME_TRACE_CACHE=$(mktemp -d) PML_MAFFT_BIN=/path/to/mafft-shim \
  python escott_e2e.py SPG1_STRSG_Olson_2014 out
```

Score one assay per build, **each with its own fresh `GEMME_TRACE_CACHE`** —
a warm cache silently skips MAFFT entirely and the run passes without either
aligner being asked anything. Compare the saved `.scores.npy` and
`.matrix.npy` by SHA-256. Sanity-check the harness by pointing
`PML_MAFFT_BIN` at a nonexistent path: the run must fail.

`SPG1_STRSG_Olson_2014` (536962 mutants) and `P53_HUMAN_Kotler_2018` (1048,
5069-sequence MSA) are the two used here.

---

## 4. Hazards — read before writing code

### `partA__align` and `A__align` disagree on two tie comparisons

The single most dangerous detail in the remaining work. Copying `align.rs`'s
inner loop — the obvious move, since it is the same DP and already ported —
gives a different alignment on every tie, at identical score. No unit test and
no score check can see it.

| recurrence | `Salignmm.c` / `align.rs` | `partSalignmm.c` |
|---|---|---|
| `*prept + *ogcp2pt * gf1vapre` vs `mi` | **`>=`** | **`>`** (`// 2018/Apr`) |
| `*prept + ogcp1va * *gf2ptpre` vs `*mjpt` | **`>=`** | **`>`** (`// 2018/Apr`) |

Both files carry both forms, one live and one commented — mirror images.
Upstream changed `partSalignmm.c` in April 2018 and left `Salignmm.c`. It is
deliberate; there is nothing to report upstream. The other two comparisons in
that loop use `>` in both files.

`partalign.rs` already has this right, with the reason at the call site.

### The traceback's "off the front" branch is warp-only

`Atracking` and `Atracking_localhom` both open with a branch on
`ifi == -warpbase && jfi == -warpbase`. With warp refused, `warpbase` is a
sentinel no entry can equal, so **the branch is dead and every border cell goes
down the ordinary path**, whose two gap-run loops handle it and whose tail still
adds that cell's constraint score. Reading the condition as "both negative" —
which is what it looks like, and what `partalign.rs` did — makes it live at
`(0, 0)`, and drops the last `part_imp_match_out_sc` from `*impwmpt`.

That leaves the alignment rows and the DP score **identical** and only
`impmatch` wrong, so nothing except a differential on the constraint score can
see it: it survived the unit tests, and a unit test written afterwards pinned
the wrong value. `verify_flh.sh` caught it on 3 of 98 cases (fixed).
`align.rs` had it right; the same shape appears in every traceback still to be
ported, so compare against `-warpbase` rather than against zero.

### An input with no bit-exact answer

`getnumlen_nogap_countn` reads an uninitialised `nN` when a record has an empty
body — reachable with a four-line input, valgrind-confirmed (upstream issue
#8). There is no bit-exact answer to reproduce. `countlen-rs` picks the defined
value, zero, which agrees with the C wherever the C is defined, and
`verify_utils.sh` scores that input `indeterminate` rather than pass or fail.
Expect more of these in `dvtditr`; handle them the same way — name the input,
prove the indeterminacy, choose deliberately, never hide it in a pass count.

### `f2cl`'s stderr is not reproducible

It prints heap pointers (upstream issue #9), different on every run of the C
against itself. `verify_utils.sh` compares `f2cl` on stdout and exit status
only.

### The misspelled macro in `partSalignmm.c`

`FASTMACHCALC` (no `T`) guards two functions, is defined nowhere, and the
`#else` branch it forces differs semantically from the branch it hides —
`+=` against `=` (upstream issue #11). The port implements the live `#else`.
Fixing the spelling upstream without also fixing the `=` would change
alignments.

---

## 5. What remains: `dvtditr`

~2650 lines. Sizes are measured per function, and the reuse column is what
actually decides the cost — see §7.

| piece | cost | state |
|---|---|---|
| `partA__align` + both tracebacks + constraint half | ~360 Rust | **bit-verified** through `verify_flh.sh` |
| `Falign_localhom` | 223 Rust (875 C) | **done, bit-verified**, 98/98 |
| `TreeDependentIteration` sequential branch | 1188 C | not started |
| — its missing dependencies | **666 C** | not started |
| `dvtditr.c` `main` reached path | ~500 C | not started |

### `partA__align` — verified

`partalign::PartAlignState::part_a_align` is now driven in anger:
`gen_flhharness.c` stages a call in the real C, `examples/flhref.rs` stages the
same one in the port, and `verify_flh.sh` compares the aligned rows, the score
and the constraint score as raw f64 bit patterns over 98 cases. All 98 agree,
and `tests/golden/flh/` replays them in `cargo test` with no C present.

Verifying it found one defect the tests could not — the traceback branch in §4 —
and nothing else.

`partSalignmm.c`'s other entry point, `partA__align_variousdist` (746 lines), is
still never reached: it needs `scoringmatrices`, which `tditeration.c` allocates
only when `specificityconsideration != 0`, and the `mafft` script passes
`dvtditr -s 0.0`. Refuse it.

### `Falign_localhom` — done

`falign::falign_localhom`, 223 lines beside `falign` and `falign_givenanchors`
in the same file, because the helpers it reuses are module-private.

**Only one configuration is implemented, and it is the only one reachable.**
A traced `mafft --auto` run (L-INS-i on `test/sample`) enters `Falign_localhom`
271 times, and every call has `kobetsubunkatsu == 1`, `alg == 'A'`,
`constraint == 2` and `scoringmatrices == NULL`. `-T`, the flag that would clear
`kobetsubunkatsu`, appears on the `disttbfast` and `pairlocalalign` lines of
`scripts/mafft` and never on the `dvtditr` one.

With `kobetsubunkatsu == 1` the FFT half does not run at all: one candidate lag
of zero, so no signal vectors, no cross-correlation, no `getKouho`, no
`blockAlign2` and no `crossscore`. What remains is `alignable_reagion` on the
groups as they stand, `mymergesort` over the anchors, the cut arrays, and per
block `getkyokaigap` + `commongappick_record` + `partA__align`.

**Do not port the FFT half by copying `falign`'s.** They are not the same code:
`Falign_localhom.c` keeps the older `n20or4or2` of 4 / 2 / 20 with one real
vector per property and one FFT each, where `Falign.c` packs polarity and volume
into a single complex vector (`n20or4or2 == 1`). Reusing `falign`'s signal stage
would change the arithmetic; writing the old one out would add unexercised
floating-point for a branch nothing runs. `dvtditr-rs` refuses `-T` instead.

Note the trace also settles the question the previous revision of this document
left open: `tbfast` never reaches `Falign_localhom` because it never calls it —
`tditeration.c` holds the only two call sites, both guarded by
`constraint == 2 && use_fft` with `alg != 'd'`.

### `TreeDependentIteration` — the expensive one

The opposite case. It calls 28 functions and **16 are absent from the port**,
mostly a tree/branch-weighting family `tbfast` never needed:

| lines | function | file |
|---|---|---|
| 157 | `treeCnv` | `treeOperation.c` |
| 68 | `SSPscore` | `mltaln9.c` |
| 65 | `classifypairs` | `tditeration.c` |
| 59 | `nodeFromABranch` | `tddis.c` |
| 54 | `OneClusterAndTheOther_fast` | `tddis.c` |
| 48 | `distFromABranch` | `treeOperation.c` |
| 45 | `countnode_int` | `mltaln9.c` |
| 44 | `calcBranchWeight` | `treeOperation.c` |
| 39 | `weightFromABranch` | `treeOperation.c` |
| 26 | `msshrinklocalhom_fast_half` | `tddis.c` |
| 22 | `substitution_hosei` | `mltaln9.c` |
| 15 | `plainscore` | `mltaln9.c` |
| 13 | `node_eff` | `tddis.c` |
| 11 | `score_calc0` | `mltaln9.c` |
| ? | `conjuctionfortbfast`, `weightimportance4` | not located by the scan |

`athread` (925 lines) is skipped **only** because this port is single-threaded
and the script passes `-C 0`. A future `-C n` needs it or an equivalent.

### `dvtditr.c`'s `main`

`fftWinSize` and `fftThreshold` are **not** set by `dvtditr`'s option loop
beyond `-w` / `-z`: `constants()` applies the defaults
(`constants.c`, `if( fftWinSize == NOTSPECIFIED )` -> `FFT_WINSIZE_D` / `_P`),
which is exactly what `Constants::new_seq` does. A traced run confirms the
sequence — `NOTSPECIFIED` after `arguments()`, 20 after `constants()`. Nothing
to reproduce beyond calling the port's `Constants` in the same place.

Not "mostly argv" — it reads the `hat3` localhom table
(`readlocalhomtable2_half` / `_target`), builds `topol` / `len` / `eff`,
handles `--subalignment` tables and `_focus_` targets, then calls
`TreeDependentIteration`.

### The argument vectors, captured from a real run

```
dvtditr  -W 0.00001 -E 0.0 -s 0.0 -C 0 -t 0 -F -l 2.7 -z 50 -b 62
         -f -1.53 -Q 100.0 -h 0 -I 16 -X 0.1 -p BAATARI2 -K 0
```

`-C 0` sequential, `-p BAATARI2`, `-F` fft on, `-I 16` iterations. Capture
others by wrapping a binary in `$MAFFT_BINARIES` with a script that records
`"$@"` and stdin, then `exec`s the real one — that is how every argv and
staged input in this document was obtained, and how a single binary is
replayed in isolation for benchmarking or diffing.

---

## 6. Order, and why

1. ~~**`Falign_localhom`**~~ — done, and it did what it was ordered first for:
   the DP is verified.
2. **The 666 lines of small dependencies** — independent, unit-testable, and
   they gate everything after.
3. **`TreeDependentIteration`'s sequential branch.**
4. **`dvtditr.c`'s `main`**, then swap in and run §3's ladder.

Do not ship a partial DP or a partial driver. Every slice above except the DP
loop and the refinement loop can be unit-tested; those two land whole.

---

## 7. What is already ported — check here before writing anything

This list exists because four estimates in this document were wrong, every one
of them from assuming instead of checking.

| C | Rust |
|---|---|
| `cpmx_calc_new`, `gapcountf`, `outgapcount` | `profile.rs` |
| `st_OpeningGapCount` / `st_FinalGapCount` | `profile::st_opening_gap_count` / `st_final_gap_count` |
| `new_OpeningGapCount` / `new_FinalGapCount` | `profile::new_opening_gap_count` / `new_final_gap_count` |
| `match_calc` | `align::match_calc` (`pub(crate)`) — **do not write a second copy** |
| `fillimp`, `movereg`, `movereg_swap`, `shrinklocalhom`, `calcimportance_half` | `localhom.rs` |
| `intergroup_score` | `score.rs` |
| `Falign`, `getKouho`, `zurasu2`, `alignableReagion`, `blockAlign2`, `mymergesort`, `getkyokaigap` | `falign.rs` |
| `commongappick`, `commongappick_record` | `addfunctions.rs` |
| `naivepairscore11`, `naivepairscorefast`, `makeskiptable` | `naivescore.rs` |
| `counteff_simple` | `compacttree.rs` |
| `WriteHat2` | `io::write_hat2` |
| MAFFT's `arguments()` option loop | `cargs.rs` |
| `partA__align`, `Atracking`, `Atracking_localhom`, `part_imp_*` | `partalign.rs` |
| `Falign_localhom` (`kobetsubunkatsu == 1` only) | `falign::falign_localhom` |

Two traps in that table, both found the hard way:

- **`io::write_hat2` indexes `mtx[i][j - i]`**, an offset layout, where the C's
  `WriteHat2_pointer` indexes a square matrix. Fill it the way the writer
  reads it or the last distance comes out zero.
- **`Falign_localhom`'s FFT stage is NOT `Falign`'s**, so `falign.rs`'s signal
  code is not reuse for it — the two files disagree on `n20or4or2` (§5). The
  helpers around it (`alignable_reagion`, `mymergesort`, `getkyokaigap`) are.
- **`partSalignmm`'s tracebacks are NOT `align.rs`'s.** Diffing the two C
  `Atracking_localhom` bodies gives 148 differing lines of ~185: `Salignmm`'s
  builds gap tables for profile reuse, `partSalignmm`'s writes straight into
  `mseq`. `partalign.rs` has its own.

---

## 8. Rules

- **Ask what a function calls, not how long it is.** `Falign_localhom` looks
  like 713 lines and is ~300 of Rust; `TreeDependentIteration` looks like 1188
  and is more. Neither was knowable from the file.
- **Never reassociate f64 arithmetic.** Loop transposition preserving each
  accumulator's summation order is fine; changing addition order is not.
- **Transliterate; do not improve.** The C's quirks are the contract —
  `version`'s missing trailing newline, `countlen`'s N tally carrying past the
  sampling cap. Where a quirk is *undefined* rather than odd, see §4.
- **Add to the harnesses, never subtract.** A wide differential against a live
  C for extending coverage, committed fixtures so the suite runs with no C.
- **Every new binary goes in `bench/coverage.sh`'s `--object` list.** Leaving
  one out reports it at 0.00% while its tests drive it — the same under-report
  that script's header warns about, pointed the other way.
- **Keep `#![forbid(unsafe_code)]` and keep the port dependency-free.**

## 9. Method

**`grep` lies about five files.** `constants.c`, `io.c`, `Lalign11.c`,
`Lalignmm.c` and `MSalignmm.c` are ISO-8859, not UTF-8, so GNU grep in a UTF-8
locale treats them as binary and reports *nothing* — not a match, not a count,
not a "binary file matches" line. A "this function is called nowhere" or "this
macro is used nowhere" conclusion drawn without `grep -a` (or `LC_ALL=C`) is
therefore worthless, and this cost an hour here: `FFT_WINSIZE_P` looked unused,
which made `constants()`'s ordinary default assignment look like memory
corruption until a watchpoint pointed straight at `constants.c`. Use `grep -a`
on this tree.

Wrap, trace, then port. `MAFFT_BINARIES` pointed at logging wrappers gives the
exact binary set and argument vector a real run uses; replaying that vector
against a staged working directory benchmarks or diffs one binary in
isolation. When output diverges, instrument a scratch copy of `core/` and
bisect on f64 bit patterns — every mismatch in this port was found that way,
never by reading.

## 10. Upstream defects found while doing this

All in unmodified MAFFT 7.526, all filed with corrected code:
[#7](https://github.com/exazyme/mafft-provolut/issues/7) Makefile drops
`$(MYCFLAGS)` (breaks GCC ≥ 14, mixes TLS storage class),
[#8](https://github.com/exazyme/mafft-provolut/issues/8) uninitialised `nN`,
[#9](https://github.com/exazyme/mafft-provolut/issues/9) `f2cl` prints heap
pointers, [#10](https://github.com/exazyme/mafft-provolut/issues/10)
`restoreu`'s record-count check is unreachable and a mismatch segfaults,
[#11](https://github.com/exazyme/mafft-provolut/issues/11) misspelled
`FASTMACHCALC` hiding an `=`/`+=` difference.
