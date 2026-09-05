# Goal: run `mafft --auto` with no C process at all

**The constraint does not move: output must not change.** `tests/matrix.rs`,
`tests/utils.rs`, `tests/flh.rs`, `verify.sh`, `verify_tbfast.sh`,
`verify_utils.sh`, `verify_flh.sh` and `verify_e2e.sh` stay green. A change that alters one alignment is a bug however
much C it removes. `../CLAUDE.md` owns the gates an agent may not lower.

This sits **beside** `GOAL.md`, which is a performance target against a frozen
baseline; its numbers are not this goal's to move.

---

## 1. Status

**Reached.** Nine binaries can run under `mafft --auto`, and all nine are Rust
and verified byte-identical. A traced `mafft --auto` on `test/sample` executes
`countlen`, `version`, `f2cl`, `tbfast` and `dvtditr` and nothing else, all of
them the Rust builds, and its output matches the all-C run byte for byte.

One strategy outside `--auto` still needs the C: `--allowshift`, which passes
`dvtditr -s 0.8` and so selects the specificity scoring path
(`scoringmatrices`, `classifypairs`, `smalldistmtx`,
`partA__align_variousdist`). `dvtditr-rs` refuses it by name. See §5.

| binary | state |
|---|---|
| `disttbfast`, `tbfast` | Rust, byte-identical (pre-existing) |
| `countlen`, `version`, `replaceu`, `restoreu`, `f2cl`, `dndpre` | Rust, byte-identical |
| `dvtditr` | Rust, byte-identical — except `--allowshift`, refused (§5) |
| `splittbfast` | **out of scope** — parttree bracket only, ≥200000 sequences |

Verified with eight Rust binaries in place: the 840-payload ProteinGym corpus
byte-identical through `mafft --auto`, upstream's own suite, and provolut-ml's
ESCOTT scorer SHA-256-identical on two assays. Speed 0.989× (51×393 windows)
and 0.981× (11×448) against the all-C build.

With `dvtditr-rs` added, the whole ladder is green:

| rung | result |
|---|---|
| `verify_dvtditr.sh` | 32/32 staged runs identical (`pre`, `order`, `infile.tree`) |
| `verify_e2e.sh` | 56 strategies identical, 2 refused (`--allowshift`) |
| upstream `test/script` | 8 of the 9 runnable cases; the ninth is `--allowshift` |
| ProteinGym payloads, rebuilt per §3 | **120/120 identical** through `mafft --auto` |
| provolut-ml's `test_mafft_builds_agree.py` | both tests pass with `PML_MAFFT_BIN` on an all-Rust build — aligner bytes **and** the GEMME/iGEMME/ESCOTT matrices, bit for bit |

The ESCOTT rung was run through provolut-ml's own test rather than by hand, and
sanity-checked the way §3 says to: pointing `PML_MAFFT_BIN` at a path that does
not exist makes it fail, so the pass is the Rust build's and not a skip.

The corpus is 120 payloads — every sub-alignment `sample_sub_alignments` draws
from `P53_HUMAN_full_04-29-2022_b09.a2m` and
`SPG1_STRSG_full_11-26-2021_b07.a2m` at the iGEMME defaults. Not the 840 of the
earlier run, which came from more assays. **It earned its keep**: the first pass
over it failed 105 of 120 and found the `pre`-written-once defect in §4, which
every other gate had passed over.

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
   `partA__align`, against the real C on 98 staged calls **plus a 60-call
   sequence run through one process on each side**. The sequence is not a
   bigger corpus: single calls cannot see what a call leaves in the scratch
   buffers, and a real run makes hundreds of calls that each land on the
   previous one's leftovers. It holds `njob` and the width constant, as
   `tditeration` does — a sequence that varies either stages something the C
   never sees, and both harnesses refuse it rather than run it.
   `./verify_td.sh ../core` — the branch-weighting family on 96 staged guide
   trees.
   `./verify_tdi.sh ../core` — the whole refinement loop on 126 staged runs.
   That one compares a driver, not a function: hundreds of accept/reject
   decisions feed each other, so one differing bit anywhere lands as a
   different alignment.
   Both need only a built `core/`, run in seconds, and `--regenerate`
   refreshes the fixtures `cargo test` replays. **This is what proved the DP**,
   and it is where the next slice should be proved too: staging one function
   against the C costs an hour and finds what no unit test can.
3. `./verify_dvtditr.sh ../binaries target/release/dvtditr-rs` — the whole
   binary, staged out of a real `mafft` run by the §9 wrapper method and
   replayed against both. Minutes, not an hour, and it is the first thing to
   run after touching anything the refinement path reaches.
   `./verify.sh`, `./verify_tbfast.sh`, `./verify_e2e.sh` — the aligners,
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

### `LENGTH OVER` is behaviour, and only one aligner has it

`partA__align` ends by checking `alloclen < resultlen || resultlen > N` and
calling `ErrorExit`. It is an error path, but a reachable one — the refinement
loop can grow an alignment past the buffer `dvtditr` sized from the longest
input (`alloclen = nlenmax * 9`) — and it is a *behavioural* difference of the
worst kind: where the C stops, a port without the check carries on and produces
output the C never would. `part_a_align` now carries it, and `alloclen` is
threaded from `IterParams`.

**`A__align` and `Falign` have the same check and the port still does not.**
Those predate this stack and are reached through `treebase`, whose callers would
all have to thread `alloclen` — a change with its own verification, against
`tbfast-rs` and `disttbfast-rs`. It is left open deliberately: on real input the
limit is nine times the longest sequence and nothing approaches it, and the
synthetic cases that did were random noise at 25% gaps. Close it when
`dvtditr-rs` lands, where `alloclen` is already in hand.

### `pre` is written once, and the flag that says so is easy to miss

`writePre` takes a `force` flag. With `USE_XCED` 0 -- which is how this builds --
its entire body is `if( force ) { rewind( prep_g ); writeData_pointer( ... ); }`.
Every call from inside the refinement loop passes **0**, so `dvtditr` writes
`pre` exactly once, from its own final call.

A port that reads the loop's `writePre( ..., 0 )` as "save progress" and writes
on every improvement gets away with it only while the final alignment is the
longest one: the C *rewinds without truncating*, so the tail of a longer
intermediate would survive. Writing on every improvement and then truncating
produced a `pre` 324 bytes shorter than the C's on `test/sample` -- same
alignment, different file. The staged binary differential caught it; nothing
else would have.

### With two sequences, `dvtditr` does not stop early

The `njob == 2` block that writes the input back and returns is inside the
`#if 0` half of `dvtditr.c`. The live branch is `weight = 0; niter = 1;` -- so a
two-sequence run still builds the tree, still writes `infile.tree`, and still
refines once, with equal weights instead of branch weights. A port that took the
disabled branch produces the right alignment and the wrong `infile.tree`, which
is what `mafft --treeout` hands the user.

### A global the harness has to set, or it stages a configuration nobody runs

`Falign` compacts each block before aligning it when `kobetsubunkatsu &&
fftkeika`. `dvtditr` sets **both** — `fftkeika = 1` is on line 51 of
`dvtditr.c`, and every other binary leaves it 0 — so the two conditions are one
condition in practice. A harness that sets `kobetsubunkatsu` and forgets
`fftkeika` stages a `Falign` no caller has: 52 of 126 refinement runs then
"failed" against a port that was right, and the fix was one line in the C
harness, not in the port. When a differential turns red after a change that
made the port *more* faithful, suspect the staging first.

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
| the branch-weighting family (`treeop.rs`) | 470 Rust | **done, bit-verified**, 96/96 |
| `TreeDependentIteration` sequential branch | 350 Rust (1188 C) | **done, bit-verified**, 126/126 |
| `dvtditr.c` `main` reached path | 300 Rust (~500 C) | **done**, 16/16 staged runs identical |

### `partA__align` — verified

`partalign::PartAlignState::part_a_align` is now driven in anger:
`gen_flhharness.c` stages a call in the real C, `examples/flhref.rs` stages the
same one in the port, and `verify_flh.sh` compares the aligned rows, the score
and the constraint score as raw f64 bit patterns over 98 cases. All 98 agree,
and `tests/golden/flh/` replays them in `cargo test` with no C present.

Verifying it found one defect the tests could not — the traceback branch in §4 —
and nothing else. The state a call leaves behind is covered too: `verify_flh.sh`
replays a 60-call sequence through one process on each side, and they agree.

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

**Its dependency list was wrong, and the correction is large enough to change
the plan.** The previous revision listed sixteen absent functions from a source
scan. A callgrind profile of the real `dvtditr`, replayed from a
wrapper-captured argument vector and stdin across L-INS-i, FFT-NS-i and E-INS-i,
says which of them a run actually enters:

| runs | ported here | note |
|---|---|---|
| `treeCnv`, `calcBranchWeight`, `distFromABranch`(`_rec`), `weightFromABranch`(`_rec`), `OneClusterAndTheOther_fast` | yes — `treeop.rs` | plus the support cast the scan missed: `isLeaf`, `searchParent`, `syntheticLength`, `negativeMember2`, `checkMinusLength`, `calcW`, `stopolInit`, `IntExistsInVec` |
| `msshrinklocalhom_fast_half` | already was | it is `fastshrinklocalhom_half` under another name — same body, one iterating by index and one by pointer. `localhom::fastshrinklocalhom_half` covers it |
| `fastconjuction_noname` | already was | `treebase::fastconjuction_noname` |

Never entered, on any of the three strategies: `SSPscore`, `classifypairs`,
`nodeFromABranch`, `countnode_int`, `node_eff`, `plainscore`,
`substitution_hosei`, `score_calc0`, `msshrinklocalhom_fast_target`,
`fastconjuction_noname_kozo`. They are the `score_check`, `--target` and kozo
paths. **Refuse them**, the way `-T` and `partA__align_variousdist` are refused;
do not port on the strength of a source scan.

Two were dead outright: `weightimportance4` has no definition anywhere in the
tree (an `extern` in `functions.h` and commented-out call sites), and
`conjuctionfortbfast`'s live counterpart in the sequential branch is
`fastconjuction_noname`, which the port already had.

`athread` (925 lines) stays out of scope, and the profile settles why rather
than assuming it: at `-C 0` the callgrind output contains
`TreeDependentIteration` and no `athread` at all. A future `-C n` needs it.

Its body is ported — `tditeration::tree_dependent_iteration`, 350 lines — for
the configuration a traced `dvtditr` reports on every strategy:
`utree=1 weight=4 score_check=1 nkozo=0 rnakozo=0 alg='A' use_fft=1 nthread=0`,
with `constraint` 2 or 0. That single line of trace output is what makes the
port small: `utree=0` alone would have brought a second tree builder
(`substitution_hosei` and four clustering methods) that nothing runs.

Everything outside it is refused with an assertion naming the setting, not
approximated: dynamic tree rebuilding, other weightings, `score_check=2`, kozo,
RNA, `--target`, `specificityconsideration`, `alg != 'A'`, `-C n`, BESTFIRST.

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

### `dvtditr-rs`

`src/bin/dvtditr.rs`, 300 lines: `dvtditr`'s own `arguments()`, the `pre` /
`hat2` / `hat3` reads, the guide tree, `order`, the segment loop and the write
back to `pre`.

Two things it does not do the way the C's `main` does, both measured rather
than assumed:

- **The tree builder is not ported twice.** `dvtditr` calls
  `fixed_musclesupg_double_treeout`, a full-matrix copy of the memsave UPGMA
  `tree::upgma_treeout` already ports. Running both C functions on 40 staged
  distance matrices — 3 to 36 sequences, with ties, zeros, near-equal and
  sub-micro distances — gives **bit-identical branch lengths in every case**,
  and the full member lists are exactly the `Treedep` reconstruction
  `treeop::expand_topol` performs. So the port keeps one builder.
- **`searchAnchors` is ported** (`treeop::search_anchors`): FFT-NS-i runs
  unconstrained with segmentation on, and then the refinement really does work
  segment by segment. Only `-l` — which every constrained call passes — turns
  it off.

**Refused, by name, at startup:** `-s` non-zero (the specificity path;
`--allowshift` is the strategy that reaches it), `-C n`, `-K n` (`--add`, which
needs the `skipthisbranch` computation), a `hat3` carrying `korh == 'k'`
structural constraints, an unconstrained input over 50000 residues (where the C
silently switches to the memsave aligner, which is a *different alignment*), and
any option letter `dvtditr` does not define.

`verify_dvtditr.sh` generates two inputs the repo's test data does not have: a
wide alignment with three conserved blocks, which is what makes `searchAnchors`
cut it into four segments -- otherwise the segment loop only ever runs once and
the concatenation across segments is never exercised -- and a two-sequence
input for the `weight = 0` path above. It compares `pre`, `order` and
`infile.tree`.

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
2. ~~**The small dependencies**~~ — done. Six functions ran, not sixteen; the
   rest are refused. See the table above.
3. ~~**`TreeDependentIteration`'s sequential branch.**~~ — done.
4. ~~**`dvtditr.c`'s `main`**, then swap in and run §3's ladder.~~ — done,
   except the ProteinGym and ESCOTT rungs (§1).

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
| `treeCnv`, `calcBranchWeight`, `distFromABranch`, `weightFromABranch`, `OneClusterAndTheOther_fast`, and their helpers | `treeop.rs` |
| `msshrinklocalhom_fast_half` | `localhom::fastshrinklocalhom_half` — **same function, different name** |
| `TreeDependentIteration` (sequential) | `tditeration::tree_dependent_iteration` |
| `counteff_simple` (the full-`topol` one) | `treeop::counteff_simple` — **not** `tree.rs`'s, which is the memsave form and clamps negative lengths |
| `fixed_musclesupg_double_treeout` | `tree::upgma_treeout` + `treeop::expand_topol` — the two C builders agree bit for bit (§5) |
| `searchAnchors` | `treeop::search_anchors` |
| `readhat2_pointer`, `input_new`, `readlocalhomtable2_half` | `io::read_hat2`, `io::read_localhomtable2_half` |
| `dvtditr.c` `main` + its `arguments()` | `src/bin/dvtditr.rs` |

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

**`grep` used to lie about eight files, and will again after an upstream
import.** Upstream ships `constants.c`, `io.c`, `mltaln9.c`, `Lalign11.c`,
`Lalignmm.c` and `MSalignmm.c` (plus one file each under `MPI/` and
`extensions/`) as EUC-JP, and GNU grep in a UTF-8 locale treats such a file as
binary and reports *nothing* for it — not a match, not a count, not a "binary
file matches" line. A "this function is called nowhere" conclusion drawn over
those files was therefore worthless, and it cost an hour here: `FFT_WINSIZE_P`
looked unused, which made `constants()`'s ordinary default assignment look like
memory corruption until a watchpoint pointed straight at `constants.c`.

They are converted now and `pre-push` refuses a source that is not valid UTF-8,
so ordinary `grep` is trustworthy again. **After importing a new upstream
release, run `scripts/check-encoding.sh --fix`** — the release brings the EUC-JP
files back, and until it is run every "used nowhere" claim over them is
unfounded.

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
