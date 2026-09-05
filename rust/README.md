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

**Eight of the nine binaries `mafft --auto` can execute are now Rust**; only
`dvtditr` is still the C's, plus `splittbfast` for the parttree bracket, which
is out of scope. Both the L-INS-i and FFT-NS-i brackets come out byte-identical
against an all-C run. See `GOAL-C-FREE.md` for the whole
target and what is left. The utilities' differential harness is
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
- **threads.** The port is single-threaded and refuses `-C > 0`. Note this is
  a "not implemented", not a "cannot be implemented": measured on 39 corpus
  inputs x 3 repeats, the C's own `disttbfast` and `tbfast` give byte-identical
  `pre`, `order`, `hat2` and `hat3` at `-C 4` as at `-C 0`, including on an
  80-sequence input where four threads are a real 3.3x speedup. The end-to-end
  `mafft --thread` difference is confined to `--maxiterate > 0`, i.e. to
  `dvtditr`, which this port does not touch. Run with `-C 0` / `-C 0-0`, which
  is what the wrapper passes by default.

Iterative refinement (`--maxiterate`) is not ported either, but it still works
end-to-end, because `dvtditr` is a separate binary fed by these ones' output.

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
| `verify_dvtditr.sh`, 16 staged `dvtditr` runs | 16/16 identical |
| upstream `test/script`, the 9 cases runnable here | 8/9; `--allowshift` is refused, not wrong |
| `cargo test` | 221/221 |

(`test/script`'s other two cases, `mafft-qinsi` and `mafft-xinsi`, need external
RNA folding programs and fail identically with the stock C binaries.)

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
