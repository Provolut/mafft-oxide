# MAFFT's alignment algorithms, in Rust

A transliteration of MAFFT 7.526's multiple-sequence-alignment algorithms from
C into Rust, written to produce **bit-identical output**, not to be fast. Where
the C and a cleaner formulation disagree, the C wins.

Four binaries, all drop-in replacements for their `libexec/mafft/` originals —
same argument vector, same input, same output files:

| binary | replaces | covers |
|---|---|---|
| `disttbfast-rs` | `disttbfast` | FFT-NS-1/2, NW-NS-1/2, `--add`, `--addprofile`, `--keeplength`, `--allowshift`, `--aamatrix`, `--text`, `--oneiteration`, `--merge`, `--treein` (plain `mafft`, `fftns`, `nwns`) |
| `tbfast-rs` | `tbfast` | L-INS-1/i, G-INS-1/i, E-INS-1/i, `--add`, `--addprofile`, `--allowshift` (`linsi`, `ginsi`, `einsi`) |
| `countlen-rs` | `countlen` | the record count, longest and shortest ungapped lengths, and sequence type the `mafft` script reads to pick a strategy |
| `version-rs` | `version` | the version string, with the missing trailing newline the script relies on |
| `replaceu-rs` | `replaceu` | folds the input to a known alphabet and stamps the `_os_<i>_oe_<name>` markers |
| `restoreu-rs` | `restoreu` | writes the original characters back onto the finished alignment's gap pattern |
| `f2cl-rs` | `f2cl` | the output writer, **FASTA path only** -- clustal and phylip are refused, not approximated |
| `dndpre-rs` | `dndpre` | re-derives the `hat2` distance matrix from a finished alignment, for the second progressive cycle |

The two aligners handle protein and nucleotide input, which the reader tells
apart the same way the C does.

**All nine binaries `mafft --auto` can execute are Rust** (`splittbfast`, the
parttree bracket, is out of scope and never runs under `--auto`), so a
`mafft --auto` run executes no C process at all: traced, it runs `countlen`,
`version`, `f2cl`, `tbfast` and `dvtditr` and nothing else. L-INS-i, FFT-NS-i,
E-INS-i and G-INS-i all come out byte-identical against an all-C run, as do 120
ProteinGym payloads and provolut-ml's GEMME/ESCOTT matrices. The one strategy
that still needs the C is `--allowshift`; see "What is *not* ported".

The utilities' differential harness is
`verify_utils.sh` (needs a built `core/`), and `tests/utils.rs` pins the same
behaviour against committed fixtures, so the suite still runs on a machine
with no C toolchain.

Three defects in the C were found by building these and are filed upstream:
issue #8 (`countlen` reads uninitialised memory), #9 (`f2cl` prints heap
pointers to stderr, making its stderr non-reproducible) and #10 (`restoreu`
segfaults on a record-count mismatch, because its own check is unreachable).
The port diverges from the C deliberately in two of those, in both cases where
the C's behaviour is undefined; each divergence is commented at the site and
scored `indeterminate` rather than as a pass by `verify_utils.sh`.

Dropped into a MAFFT installation, upstream's own test suite passes 9/9 and
`mafft --auto` is byte-identical. The `-i` (iterative refinement) variants come
out identical too: `dvtditr` is a separate binary that consumes these ones'
output, and it gets byte-identical `pre`, `hat2` and `hat3` to work from.

```
cd rust && cargo build --release
./target/release/disttbfast-rs <same args as disttbfast> < in.fa > out.fa
./target/release/tbfast-rs     <same args as tbfast>     < in.fa   # writes `pre`
```

## What is ported

| stage | C source | Rust |
|---|---|---|
| FASTA reading, type detection, writing | `io.c` | `io.rs` |
| BLOSUM tables (30/45/50/62/80/90/100) | `blosum.c` | `blosum.rs` |
| JTT and transmembrane substitution counts | `JTT.c` | `jtt.rs` |
| score matrix (BLOSUM, JTT/TM, nucleotide), gap penalties, `--fmodel`, FFT scales, `makedynamicmtx` | `constants.c`, `mltaln9.c`, `miyata.h`, `DNA.h` | `constants.rs` |
| 6-mer ("ktuples") distance matrix, protein and nucleotide | `disttbfast.c`, `mltaln9.c` | `dist.rs` |
| UPGMA guide tree, sequence weights | `mltaln9.c` | `tree.rs` |
| profile statistics | `tddis.c`, `mltaln9.c` | `profile.rs` |
| **group-to-group DP** | `Salignmm.c` (`A__align`, `Atracking`, `Atracking_localhom`) | `align.rs` |
| 1-vs-1 global DP | `Galign11.c` (`G__align11`) | `align.rs` |
| **1-vs-1 local DP** | `Lalign11.c` (`L__align11`) | `lalign.rs` |
| **1-vs-1 generalised affine-gap local DP** | `genalign11.c` (`genL__align11`) | `genalign.rs` |
| **FFT** | `fft.c` | `fft.rs` |
| **FFT anchoring** | `Falign.c`, `fftFunctions.c` | `falign.rs` |
| **all-pairs pairwise stage** | `pairlocalalign.c` (`pairalign`) | `pairalign.rs` |
| **position-specific constraints** | `io.c`, `mltaln9.c`, `tddis.c` | `localhom.rs` |
| shift ("warp") DP | the `trywarp` blocks of `Salignmm.c` / `Galign11.c` / `Lalign11.c` | `warp.rs` |
| progressive driver | `disttbfast.c` + `tbfast.c` (`treebase`) | `treebase.rs` |
| MSA-derived distances for the second cycle | `mltaln9.c` | `score.rs` |
| argument parsing | `disttbfast.c`, `tbfast.c`, `pairlocalalign.c` | `params.rs`, `bin/tbfast.rs` |

Four pieces carry most of the difficulty.

**`A__align`** is the affine-gap DP over two *profiles* with position-specific
gap-open and gap-extend costs. Three callers take three different branches of
the C: `treebase` supplies a profile-history slot (`cpmxhist`), which lets a
parent node reuse a child's composition matrix rather than recomputing it;
`Falign` supplies boundary-gap patterns describing what lies just outside the
block being aligned; and with constraints it also takes a per-cell match bonus
and a different traceback. The cache is not an optimisation detail — using it
changes the summation order, so a port that skipped it would not be bit-exact.

**`Falign`** is the step MAFFT is named for (Katoh, Misawa, Kuma & Miyata,
*NAR* 30:3059–3066, 2002). Each group becomes a complex signal whose real and
imaginary parts are the weighted, normalised polarity and volume of each
column; the cross-correlation of the two signals, computed by FFT, has peaks at
the relative shifts where homologous regions line up. At each of the 20 best
shifts a sliding window scores the two groups column by column and keeps runs
that beat a threshold; those are the candidate anchors. `blockAlign2` then
picks a consistent, monotone subset of them by dynamic programming, and the
groups are cut at the surviving anchors and aligned block by block.

**`L__align11`** is the local (Smith-Waterman style) pairwise DP that L-INS-i
runs on every pair of input sequences — the function this repository's
profiling put 68% of the runtime in. Its result feeds both the guide tree and
the constraints.

**Nucleotides** take a different branch almost everywhere. The score matrix is
not tabulated but generated: Kimura's two-parameter model at a
transition/transversion ratio of 2, raised to the 200th power, taken through
`log10`, recentred, rescaled so the mean diagonal is 600, offset and rounded --
after which `'u'` is given `'t'`'s row, the upper-case half of the alphabet is
made a copy of the lower-case half, and each IUPAC ambiguity code gets the mean
of the bases it stands for. The k-tuple distance switches to base 4 with its own
length correction, the reader lower-cases instead of upper-casing, and the FFT
signal changes shape (see the note below). `tests/golden/constants_dna.txt` pins
the resulting matrix against the C, because none of that chain is obvious enough
to trust otherwise.

**`genL__align11`** is what E-INS-i runs instead, and it earns its own module
because of a third recurrence. Besides "gap in seq1" and "gap in seq2" it keeps
a *pair* of unaligned runs, which can be entered from either side for one flat
`penalty_OP` however long the two runs grow. That is what lets E-INS-i leave
large, mutually unalignable insertions out of the alignment rather than paying
extension across them, and it is why its traceback needs two index matrices
instead of one: a cell has to record both the row and the column it came from,
because that transition moves in both axes at once.

**The constraints** ("localhom") are what makes L-INS-i and G-INS-i
consistency-based. Each pairwise alignment is cut at its gaps into ungapped
blocks; `calcimportance_half` reweights every block by how much of the rest of
the data set supports the same region; and during the progressive alignment
`fillimp` paints those weights onto a matrix that `A__align` adds to its match
scores.

## What is *not* ported

Refused loudly rather than approximated, so an unsupported run fails instead of
quietly returning a different alignment:

- external anchors (`--anchors`) whose end positions run past the sequence's
  residue count: `readexternalanchors` range-checks only the two *start*
  positions, and even those with `>` rather than `>=`, so an over-long end
  reads `recountpositions`' map above the entries it filled and the anchor
  comes back ending before it starts. `Falign_givenanchors` writes that
  inversion into its cut list and copies the block with a negative `strncpy`
  length; the C segfaults. Anchors are otherwise ported. (An earlier draft of
  this list called the whole feature unreproducible because `anchcomp` and
  `anchscorecomp` break ties on the *addresses* of the elements being
  permuted. Measured rather than assumed: glibc's `qsort` is a mergesort at
  these element sizes, so those tie-breaks reduce to stable and
  reverse-stable, which a Rust stable sort gives.)
- `--treein auto` / `test` where a picked sequence is gap-free and of maximal
  length: `varpairscore` gives each `pickseq` row exactly `nlenmax` bytes and
  `gappick0` writes its terminator one past that, into memory the allocator
  then reuses, after which the C aligns whatever follows the row. Measured: it
  reports 280 residues before the loop and aligns 282 inside it. No defined
  answer to match.
- `--oneiteration` past the first node where the C would reach `A__align`:
  there the C segfaults (it hands `dooneiteration` a freed, NULL `cpmxhist`
  and then writes through `cpmxhist + l`), so there is no output to match and
  the port refuses instead. 78 of the 88 corpus inputs never reach it.
- RNA secondary-structure strategies (`--qinsi`, `--xinsi`), which call out to
  external programs
- `-@` (`--sp`), the generalized-affine DP in `D__align`.
  The memsave aligner (`-M`, `--memsave`, and the switch the C makes on its
  own past 30000 columns) is ported, both its DP and its FFT anchoring, and so
  is `-&` (`Aalign`).
- `--parttree` / `--dpparttree` (a separate binary, `splittbfast`)
- **threads.** The port is single-threaded and refuses `-C > 0`. For
  `disttbfast` and `tbfast` that is a "not implemented", not a "cannot be":
  measured on 39 corpus inputs x 3 repeats, the C's own binaries give
  byte-identical `pre`, `order`, `hat2` and `hat3` at `-C 4` as at `-C 0`,
  including on an 80-sequence input where four threads are a real 3.3x speedup.
  For `dvtditr` it is closer to "cannot": the C's threaded refinement is
  **nondeterministic** — three `--thread 4` runs of one input give three
  different alignments — so there is no fixed reference to be byte-identical
  to. A callgrind profile confirms `athread` never runs at `-C 0`. Run with
  `-C 0` / `-C 0-0`, which is what the wrapper passes by default.

The six configurations that still need the C, and what porting each would take,
are in [`WHAT-STILL-NEEDS-THE-C.md`](WHAT-STILL-NEEDS-THE-C.md). One of them —
`--thread n` — is not a porting task at all: the C's own threaded `dvtditr` is
nondeterministic, so there is no fixed answer to match.

`dvtditr` -- iterative refinement -- **is** ported, for the configuration
`scripts/mafft` produces: `utree=1 weight=4 score_check=1 nkozo=0 alg='A'
use_fft=1`, with a `hat3` table or without. It refuses, by name at startup:

- `-s` non-zero, the specificity scoring path (`scoringmatrices`,
  `classifypairs`, `smalldistmtx`, `partA__align_variousdist`). `--allowshift`
  is the strategy that reaches it, and the only `mafft` strategy that still
  needs the C `dvtditr`.
- `-T` (`kobetsubunkatsu = 0`), which turns the FFT anchoring back on inside
  `Falign_localhom`. `scripts/mafft` passes it to `disttbfast` and
  `pairlocalalign`, never to `dvtditr`.
- `-C n`, `-K n` (`--add`, which needs the `skipthisbranch` computation), a
  `hat3` carrying `korh == 'k'` structural constraints, `score_check=2`, kozo
  and RNA weighting, `--target`, and `utree=0` (rebuilding the tree from the
  current alignment each round -- a second tree builder nothing runs).
- an unconstrained input over 50000 residues, where the C silently switches to
  the memsave aligner. That is a *different alignment*, not a cheaper route to
  the same one, so the port stops rather than quietly disagree.

## Verification

Bit-exactness is measured, not argued. Three differential harnesses plus a
suite of unit tests.

```
cd core && make                       # build the C reference

rust/verify.sh         core/disttbfast rust/target/release/disttbfast-rs
rust/verify_tbfast.sh  core/tbfast     rust/target/release/tbfast-rs
rust/verify_flh.sh     core                # Falign_localhom, function level
rust/verify_td.sh      core                # branch weighting, function level
rust/verify_tdi.sh     core                # the refinement loop, whole driver
rust/verify_dvtditr.sh binaries rust/target/release/dvtditr-rs   # the binary
rust/verify_lalign_kernel.sh                       # AVX-512 kernel vs scalar
TBFAST_RS=rust/target/release/tbfast-rs MAFFT=/path/to/bin/mafft \
  rust/verify_e2e.sh /path/to/libexec/mafft rust/target/release/disttbfast-rs \
  test/sample test/samplerna
cargo test --manifest-path rust/Cargo.toml
```

`verify_flh.sh` works one level down, because the binary that would exercise
`Falign_localhom` — `dvtditr` — is not ported yet: `gen_flhharness.c` links the
real `core/` objects and stages one call directly, `examples/flhref.rs` stages
the same one in the port, and the two are compared on the aligned rows and on
the score and constraint score as f64 bit patterns. `gen_flhcases.py` builds the
corpus: the group shapes `tditeration` produces (one member against the rest,
even splits) over three alignment widths, both terminal-gap settings, legacy gap
cost, anchoring parameters on and off their traced defaults, both constraint
orientations, and the corners — no constraints at all, a single column, no gaps,
identical rows. `--regenerate` writes the C's answers into `tests/golden/flh/`,
which `tests/flh.rs` replays with no C present.

`verify_td.sh` does the same for `treeop.rs`: it stages a guide tree — a merge
sequence with its branch lengths — and compares every branch weight, every
distance-from-a-branch and path-weight vector, and the group split each
refinement step starts from. Its corpus covers random, ladder and balanced
trees from 3 to 36 sequences, with ordinary, zero, negative and sub-`MINLEN`
branch lengths: the zeros and negatives are the point, since they drive
`checkMinusLength`'s clamp and the two early returns in `calcW` that are the
only non-arithmetic paths in the family.

`verify.sh` and `verify_tbfast.sh` run the C binary and the Rust one on the
same inputs with the same arguments and require every output file to be
byte-identical — stdout and `order` for `disttbfast`, and `pre`, `order`,
`hat2` and `hat3` for `tbfast`. Both generate their corpus with
`gencorpus.py`: `test/sample` and `test/samplerna`, random subsets, sequences
truncated to 1–40 residues (which trip the "too short for a 6-mer" sentinel),
pairs, exact duplicates (ties everywhere in the tree), ragged lengths, ambiguity
codes, sequences containing `X` (which the score-only variants strip),
pre-existing gaps, random sequences with no phylogenetic signal, over-long
names, lower case, `*` characters, and tandem repeats, shuffled domains and
inversions. That last group matters twice over: it is the only kind of input on
which the FFT scan proposes mutually inconsistent anchors and `blockAlign2` has
to prune, and it is what `--allowshift` exists for. Half the corpus is
nucleotide -- families at four divergences, IUPAC codes, `u` for `t`, upper
case, repeats, and sequences below the tuple length -- because every stage
branches on the sequence type. Each input is crossed with every argument vector
the `mafft` wrapper is known to emit for the ported strategies. `verify.sh`
then does a second pass for `--add` / `--addprofile`: the C is run first to
produce a real alignment, `gen_add.py` splits it at several points, and both
binaries are asked to add the tail back.

Both harnesses resolve their two binary arguments to absolute paths and refuse
to run a non-executable one. That is not a nicety: they run each case from a
scratch directory, so relative paths made *both* binaries exit 127, no files
were produced, and every case scored a pass. Every count these harnesses
reported before that guard was added is worthless.

`verify_e2e.sh` does the same one level up: it runs the real `mafft` wrapper
against a stock libexec directory and against a copy with the Rust binaries
swapped in, across `--nofft` and FFT variants, `--localpair` and `--globalpair`
at several `--maxiterate` settings, `--allowshift`, `--auto`, `--reorder`,
`--clustalout` and `--inputorder`. It counts a case as failed when the
*reference* run fails, so a missing wrapper cannot make the comparison pass
vacuously.

`cargo test` needs no C build. Beyond four golden alignments and two bit-exact
tables dumped from the C (the derived score matrix and the FFT sine table), it
checks properties that must hold on *every* input, which is what catches the
cases a corpus misses:

- alignments preserve their inputs — ungapping any output row returns the
  input row, no column is gap-gap, both groups come out the same width;
- the two routines that compute the same thing agree — `L__align11` against
  `L__align11_noalign`, `naivepairscore11` against `naivepairscorefast`, and
  each incremental profile update (`cpmx_calc_add`, `st_OpeningGapAdd`, …)
  against recomputing from scratch;
- the algebra holds — profile columns sum to the total weight, sequence
  weights sum to 1, distances land in range, identical sequences get identical
  weights and distance 0, local alignment is symmetric in its arguments and
  bounded by the self scores;
- the FFT round-trips, and the correlation of a signal with a shifted copy of
  itself peaks at exactly that shift — the property the anchoring depends on;
- constraint blocks are ungapped, in range, disjoint and ordered;
- the reader's documented quirks hold — name truncation at `B-2`, `*` becoming
  a gap, non-alphabetic characters dropped, sequence-type detection;
- unsupported configurations are refused.

Current status, on this source tree:

| check | result |
|---|---|
| `verify.sh`, generated corpus (88 inputs × 36 argument sets, plus `--treein`) | 3520/3520 identical |
| `verify_tbfast.sh`, generated corpus (40 inputs × 19, plus `--treein`) | 920/920 identical |
| `verify_e2e.sh`, 2 inputs × 22 wrapper strategies, all three aligners in Rust | 56/58 identical, 2 refused (`--allowshift`) |
| `verify_flh.sh`, 98 staged `Falign_localhom` calls + a 60-call sequence | 98/98 identical, sequence identical |
| `verify_td.sh`, 96 staged guide trees through the branch-weighting family | 96/96 identical |
| `verify_tdi.sh`, 126 staged refinement runs | 126/126 identical |
| `verify_dvtditr.sh`, 32 staged `dvtditr` runs | 32/32 identical |
| `verify_lalign_kernel.sh`, the AVX-512 kernel against the scalar path | 15/15 identical |
| 120 ProteinGym payloads through `mafft --auto`, all binaries Rust | 120/120 identical |
| upstream `test/script`, the 9 cases runnable here | 8/9; `--allowshift` is refused, not wrong |
| `cargo test` | 221/221 |

(`test/script`'s other two cases, `mafft-qinsi` and `mafft-xinsi`, need external
RNA folding programs and fail identically with the stock C binaries.)

## Speed

The port was written for bit-exactness rather than speed, and for most of its
life it ran at roughly the C's pace. It is now the fastest of the three builds
on the workload that matters here — ESCOTT's `mafft --auto --preservecase`, 12
ProteinGym payloads:

| build | time | vs the C at `-O3` |
|---|---|---|
| C, `gcc -O3` | 12.60 s | 1.00x |
| C, `-march=native` (its AVX-512 kernel compiled in) | 7.34 s | 1.72x |
| **this port** | **6.86 s** | **1.84x** |

Two things got it there, and `bench/OPTIMISATION-LOG.md` has the measurements
for both. The smaller one removed bounds checks and a sentinel test from
`dvtditr`'s hot loops (18% off that binary). The larger one is
`lalign_avx512` — a transliteration of the C's own AVX-512 kernel for
`L__align11`, which is ~69% of the aligner's cost. On one staged input it takes
`tbfast` from 0.770 s to 0.292 s, 2.6x, and lands 17% ahead of the C's
`-march=native` build of the same function.

**The port's kernel dispatches at runtime; the C's is compiled in.** A C build
carrying it (`-march=native`) crashes on a machine without AVX-512, and a
portable C build never gets it — which is why the conda-forge 7.526 that
production runs is the slowest build in the table. One `tbfast-rs` is fast where
the hardware allows and correct everywhere else, and
`MAFFT_RS_NO_AVX512=1` forces the scalar path for bisecting.

## Notes for anyone extending this

Thirteen things about the C are load-bearing and easy to "clean up" by accident.
Most were found the hard way, by a corpus case rather than by reading.

1. **Scratch buffers are never cleared.** `A__align` and friends keep `static`
   arrays that grow but are not zeroed, and several loops index one element
   past the current sequence length — `ogcp2[lgth2]`, `ogcp1[1]` when
   `lgth1 == 1`. The Rust keeps the same buffers with the same growth rule so
   the same stale values are read back.
2. **Off-by-one reads of the gap tables are intentional.**
   `creategapfreqresult` and `createfgresult` read `gaptable[alen]`, the NUL
   terminator, and rely on it not being `'-'`. `G__align11` and `L__align11`
   read `seq1[lgth1]`, the terminator again, on their last row; that row of the
   score matrix is all zeros.
3. **`blockAlign2` compares against the wrong variable, and it matters.** The
   loop that fills `maxi` tests candidates against `maxj` — a `static` holding
   the *previous* cell's value, already penalised by `+= penalty`. Writing that
   as a local (the obvious reading) changes which anchors survive, and shows up
   only on inputs with repeats. Its `result1[]` traceback buffer is also read
   one element below where the current call wrote.
4. **`fillimp` reads a constraint backwards half the time.** A constraint is
   stored once, under the lower-numbered sequence; when the group-1 member has
   the higher original index, `movereg_swap` reads its two halves the other way
   round. Missing that shifts the bonus by a few columns — it first showed up
   at tree step 13 of 35 on `test/sample`.
5. **`new_FinalGapCount` cannot use its `egappat` argument.** The step after
   its loop has already read the NUL terminator into `gb`, so the test is
   always false; the compiled function is identical to `st_FinalGapCount`.
6. **Globals leak between the two stages of `tbfast`.** The pairwise stage's
   argument parser sets `distout = 'h'` unconditionally and `tbfast`'s never
   resets it, so `hat2` is written on every run — including `-+ 0`, where
   nothing asked for it.
7. **The nucleotide FFT correlates one base, not four.** `Falign` sets
   `n20or4or2` to 1 for `scoremtx == -1`, so `seq_vec_3`'s guard
   `n < n20or4or2` admits only alphabet index 0 -- `a`. The four-vector form
   exists, but in `Fgetlag`, which the ported strategies never call.
8. **The two `_noalign` scorers charge extensions differently.**
   `G__align11_noalign` skips the extension penalty on its last row and its
   last column (a 2018 change); `L__align11_noalign`, which looks otherwise
   identical, charges it everywhere. Copying one into the other silently
   shifts every distance that goes through the `X`-stripped path.
9. **Tie-breaking depends on traversal order.** The UPGMA search skips the last
   live element of its linked list and keeps `im` across iterations; the DPs
   use `>=` in some comparisons and `>` in others; `Atracking` picks the first
   equal-scoring end cell and `Atracking_localhom` the last. Every one of those
   decides ties, and ties are common — the duplicate-sequence corpus exists to
   catch exactly this.

10. **`disttbfast` and `tbfast` are not two copies of one program.** Their
   `treebase()` and their `--add` setup look alike and are not: `tbfast`
   divides the added sequences' weights by 100 and skips the
   per-added-sequence `samemembern` pass that `disttbfast` runs (in `tbfast`
   both variants of it sit inside one `#if 0 ... #elseif 0 ... #endif`, and
   `#elseif` is not a directive, so it does not reopen the block). Read the
   copy you are porting, not the one you ported last week. The same goes for
   the two definitions of `samemember()` — the live one compares as sets, the
   `#if 0` one position by position.

11. **The anchor code writes one field twice and reuses one buffer.**
   `checkanchors_strongestfirst` invalidates an inconsistent anchor with
   `a[q].starti = a[q].startj = a[q].startj = a[q].endj = -1` — `startj`
   assigned twice, `endi` never — so a dead anchor keeps its end position, and
   which fields are -1 is what the callers test. `recountpositions` allocates
   its residue-to-column `map` once, outside the loop over the group's
   members, and never clears it, so a position above a short member's residue
   count reads whatever a longer member left there. Both are reproduced;
   neither is a rounding detail, since they decide which anchors survive.

12. **The two group aligners are not the same recurrence.** `Salignmm.c`
   (`A__align`) compiles with `USE_PENALTY_EX 1`, `MSalignmm.c` (`MSalignmm`)
   with 0 — the memsave aligner never charges a gap extension on its running
   scores. And `MSalignmm.c:2459` hands `createogresult` its *local*
   `ogcp1o`/`fgcp1o` where `Salignmm.c:2119` hands over `ogcp1opt`/`fgcp1opt`,
   so a node that reused a child's cached profile passes gap rows of zeros up
   to its parent. Together these make the memsave path a different alignment,
   not a cheaper route to the same one: forcing it with `-M -B` moves 17 of
   the 88 corpus inputs, and the second defect alone moves 5.

13. **The tracebacks' "off the front" branch is warp-only.** Both
   `Atracking` variants open on `ifi == -warpbase && jfi == -warpbase`, and
   with warp refused no entry can equal that sentinel: the branch is dead, and
   the border cells go down the ordinary path, whose gap-run loops handle them
   and whose tail still adds the cell's constraint score. Reading the test as
   "both negative" makes it live at `(0, 0)` and silently drops that last
   `part_imp_match_out_sc` — same rows, same DP score, wrong `impmatch`. Only a
   differential on the constraint score sees it (`verify_flh.sh`, 3 of 98).

One portability note: `sin(PI/n)` in `make_sintbl` is the only libm
transcendental the port depends on. Its argument is exact (`PI` scaled by a
power of two) and Rust's `f64::sin` is the platform libm, so it agrees here;
`tests/golden/sintbl.txt` pins every bit of the resulting table so a platform
where it does not agree fails loudly rather than silently realigning.

## Hazards in the C, and why the port looks the way it does

Written while porting `dvtditr`; every one of them was a defect first.

### `partA__align` and `A__align` disagree on two tie comparisons

The single most dangerous detail in this port. Copying `align.rs`'s
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

## Rules

- **Ask what a function calls, not how long it is.** `Falign_localhom` looks
  like 713 lines and is ~300 of Rust; `TreeDependentIteration` looks like 1188
  and is more. Neither was knowable from the file.
- **Never reassociate f64 arithmetic.** Loop transposition preserving each
  accumulator's summation order is fine; changing addition order is not.
- **Transliterate; do not improve.** The C's quirks are the contract —
  `version`'s missing trailing newline, `countlen`'s N tally carrying past the
  sampling cap. Where a quirk is *undefined* rather than odd, see the hazards
  above.
- **Add to the harnesses, never subtract.** A wide differential against a live
  C for extending coverage, committed fixtures so the suite runs with no C.
- **Every new binary goes in `bench/coverage.sh`'s `--object` list.** Leaving
  one out reports it at 0.00% while its tests drive it — the same under-report
  that script's header warns about, pointed the other way.
- **The crate is `#![deny(unsafe_code)]`, and exactly one module lifts it.**
  `lalign_avx512` carries `#![allow(unsafe_code)]` because explicit
  vectorisation on stable Rust means `core::arch` intrinsics, which are unsafe
  by definition (calling one on a CPU without the feature is UB); portable SIMD
  is still nightly-only. `deny` rather than `forbid` is what makes that single
  exemption expressible — and it is the only one. Anything else that reaches
  for `unsafe` should be measured first: the rest of this port is as fast as the
  C without it. **Keep the port dependency-free.**

## Method

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

## Upstream defects found while doing this

All in unmodified MAFFT 7.526, all filed with corrected code:
[#7](https://github.com/exazyme/mafft-provolut/issues/7) Makefile drops
`$(MYCFLAGS)` (breaks GCC ≥ 14, mixes TLS storage class),
[#8](https://github.com/exazyme/mafft-provolut/issues/8) uninitialised `nN`,
[#9](https://github.com/exazyme/mafft-provolut/issues/9) `f2cl` prints heap
pointers, [#10](https://github.com/exazyme/mafft-provolut/issues/10)
`restoreu`'s record-count check is unreachable and a mismatch segfaults,
[#11](https://github.com/exazyme/mafft-provolut/issues/11) misspelled
`FASTMACHCALC` hiding an `=`/`+=` difference.
