# Extending mafft-oxide

Read this before changing anything under `src/`.

This is a transliteration, and the C's quirks are the contract. Fourteen of
them are load-bearing and easy to "clean up" by accident; most were found by a
corpus case rather than by reading the code, and tidying any of them changes
alignments. The refusals below are each backed by a measurement too — several
are places the C has no defined answer to reproduce — so removing one is a
porting job, not a flag change.

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

The configurations that still need the C are tracked as issues, one per gap,
indexed by the repository's "what is still not done" issue — each carries what
porting it would take and what verifying it would need. One of them —
`--thread n` — is not a porting task at all: the C's own threaded `dvtditr` is
nondeterministic, so there is no fixed answer to match, and the target has to be
this port's own `-C 0` rather than the C's threaded output.

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

## Notes for anyone extending this

Fourteen things about the C are load-bearing and easy to "clean up" by accident.
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
14. **`setmark_clustal`'s nucleotide table declares one group it never uses.**
   `core/f2cl.c:35-38` assigns `weaker[0] = "AG"`, `weaker[1] = "CT"`,
   `weaker[2] = "CU"` and then sets `nweaker = 2`, so the loop stops before
   `"CU"` and a column of `C` and `U` gets no mark. Counting the assignments
   instead of reading `nweaker` — which is what the code looks like it means —
   marks that column `.` and changes the conservation row on nucleotide input.
   It was 4 of 688 cases on the first differential run of the clustal writer,
   all on `d_rna.fa`, all in the mark row only. The protein branch has no such
   gap: its `nstrong`/`nweaker` match its nine and eleven groups.

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
body — reachable with a four-line input, valgrind-confirmed (this repo's issue
#8). There is no bit-exact answer to reproduce. `countlen-rs` picks the defined
value, zero, which agrees with the C wherever the C is defined, and
`verify_utils.sh` scores that input `indeterminate` rather than pass or fail.
Expect more of these in `dvtditr`; handle them the same way — name the input,
prove the indeterminacy, choose deliberately, never hide it in a pass count.

### `f2cl`'s stderr is not reproducible

It prints heap pointers (this repo's issue #9), different on every run of the C
against itself. `verify_utils.sh` compares `f2cl` on stdout and exit status
only.

### The misspelled macro in `partSalignmm.c`

`FASTMACHCALC` (no `T`) guards two functions, is defined nowhere, and the
`#else` branch it forces differs semantically from the branch it hides —
`+=` against `=` (this repo's issue #11). The port implements the live `#else`.
Fixing the spelling upstream without also fixing the `=` would change
alignments.


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
  exemption expressible — and it is the only one. The module's single entry
  point, `row_checked`, is a safe function and therefore *behaves* like one: it
  asks the CPU for the feature and asserts the buffer lengths the kernel will
  read, rather than documenting them as the caller's problem. A safe function
  that ordinary arguments can drive into undefined behaviour would leave the
  crate's `deny` meaning nothing. Anything else that reaches for `unsafe`
  should be measured first: the rest of this port is as fast as the C without
  it. **Keep the port dependency-free.**


## Upstream defects found while doing this

All in unmodified MAFFT 7.526. Each is written up in this repository's issue
tracker with a reproduction and corrected code. None has been reported to the
MAFFT authors yet; each issue carries that as a to-do:
[#4](https://github.com/Provolut/mafft-oxide/issues/4) three memory-safety
bugs (`--oneiteration`, `--treein auto`/`test`, `--anchors`),
[#7](https://github.com/Provolut/mafft-oxide/issues/7) Makefile drops
`$(MYCFLAGS)` (breaks GCC ≥ 14, mixes TLS storage class),
[#8](https://github.com/Provolut/mafft-oxide/issues/8) uninitialised `nN`,
[#9](https://github.com/Provolut/mafft-oxide/issues/9) `f2cl` prints heap
pointers, [#10](https://github.com/Provolut/mafft-oxide/issues/10)
`restoreu`'s record-count check is unreachable and a mismatch segfaults,
[#11](https://github.com/Provolut/mafft-oxide/issues/11) misspelled
`FASTMACHCALC` hiding an `=`/`+=` difference.

