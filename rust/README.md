# MAFFT's alignment algorithms, in Rust

A transliteration of MAFFT 7.526's multiple-sequence-alignment algorithms from
C into Rust, written to produce **bit-identical output**, not to be fast. Where
the C and a cleaner formulation disagree, the C wins.

Two binaries, both drop-in replacements for their `libexec/mafft/` originals —
same argument vector, same input, same output files:

| binary | replaces | strategies |
|---|---|---|
| `disttbfast-rs` | `disttbfast` | FFT-NS-1/2, NW-NS-1/2 (plain `mafft`, `fftns`, `nwns`) |
| `tbfast-rs` | `tbfast` | L-INS-1/i, G-INS-1/i, `--allowshift` (`linsi`, `ginsi`) |

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
| BLOSUM62 tables | `blosum.c` | `blosum.rs` |
| score matrix, gap penalties, FFT scales, `makedynamicmtx` | `constants.c`, `mltaln9.c`, `miyata.h` | `constants.rs` |
| 6-mer ("ktuples") distance matrix | `disttbfast.c`, `mltaln9.c` | `dist.rs` |
| UPGMA guide tree, sequence weights | `mltaln9.c` | `tree.rs` |
| profile statistics | `tddis.c`, `mltaln9.c` | `profile.rs` |
| **group-to-group DP** | `Salignmm.c` (`A__align`, `Atracking`, `Atracking_localhom`) | `align.rs` |
| 1-vs-1 global DP | `Galign11.c` (`G__align11`) | `align.rs` |
| **1-vs-1 local DP** | `Lalign11.c` (`L__align11`) | `lalign.rs` |
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

**The constraints** ("localhom") are what makes L-INS-i and G-INS-i
consistency-based. Each pairwise alignment is cut at its gaps into ungapped
blocks; `calcimportance_half` reweights every block by how much of the rest of
the data set supports the same region; and during the progressive alignment
`fillimp` paints those weights onto a matrix that `A__align` adds to its match
scores.

## What is *not* ported

Refused loudly rather than approximated, so an unsupported run fails instead of
quietly returning a different alignment:

- nucleotide input, matrices other than BLOSUM62, `--fmodel`
- `--add` / `--addprofile`, sub-alignments, `--treein`, external anchors
- `--anysymbol` for nucleotides (it is a no-op for protein and is accepted)
- `--legacygapcost`, `--keeplength`
- the memory-saving aligner (`-M`) and the other `alg` values
- `--parttree` / `--dpparttree` (a separate binary, `splittbfast`)
- **threads.** MAFFT's threaded path is a different computation, not a race:
  it does not agree with the single-thread result even at `--thread 1`. Run
  with `-C 0` / `-C 0-0`, which is what the wrapper passes by default.

Iterative refinement (`--maxiterate`) is not ported either, but it still works
end-to-end, because `dvtditr` is a separate binary fed by these ones' output.

## Verification

Bit-exactness is measured, not argued. Three differential harnesses plus a
suite of unit tests.

```
cd core && make                       # build the C reference

rust/verify.sh         core/disttbfast core/../rust/target/release/disttbfast-rs
rust/verify_tbfast.sh  core/tbfast     rust/target/release/tbfast-rs
TBFAST_RS=rust/target/release/tbfast-rs \
  rust/verify_e2e.sh /path/to/libexec/mafft rust/target/release/disttbfast-rs test/sample
cargo test --manifest-path rust/Cargo.toml
```

`verify.sh` and `verify_tbfast.sh` run the C binary and the Rust one on the
same inputs with the same arguments and require every output file to be
byte-identical — stdout and `order` for `disttbfast`, and `pre`, `order`,
`hat2` and `hat3` for `tbfast`. Both generate their own corpus: `test/sample`,
random subsets of it, sequences truncated to 1–40 residues (which trip the "too
short for a 6-mer" sentinel), pairs, exact duplicates (ties everywhere in the
tree), ragged lengths, ambiguity codes, sequences containing `X` (which the
score-only variants strip), pre-existing gaps, random sequences with no
phylogenetic signal, over-long names, lower case, `*` characters, and tandem
repeats, shuffled domains and inversions. That last group matters twice over:
it is the only kind of input on which the FFT scan proposes mutually
inconsistent anchors and `blockAlign2` has to prune, and it is what
`--allowshift` exists for. Each input is crossed with every argument vector the
`mafft` wrapper is known to emit for the ported strategies.

`verify_e2e.sh` does the same one level up: it runs the real `mafft` wrapper
against a stock libexec directory and against a copy with the Rust binaries
swapped in, across `--nofft` and FFT variants, `--localpair` and `--globalpair`
at several `--maxiterate` settings, `--allowshift`, `--auto`, `--reorder`,
`--clustalout` and `--inputorder`.

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
| `verify.sh`, generated corpus (70 inputs × 15 argument sets) | 1050/1050 identical |
| `verify.sh`, 7 larger inputs up to 300 sequences × 15 | 105/105 identical |
| `verify_tbfast.sh`, generated corpus × 9 argument sets | 216/216 identical |
| `verify_e2e.sh`, 2 inputs × 19 wrapper strategies | 38/38 identical |
| upstream `test/script` with both Rust binaries | 9/9 |
| `cargo test` | 50/50 |

## Notes for anyone extending this

Eight things about the C are load-bearing and easy to "clean up" by accident.
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
7. **The two `_noalign` scorers charge extensions differently.**
   `G__align11_noalign` skips the extension penalty on its last row and its
   last column (a 2018 change); `L__align11_noalign`, which looks otherwise
   identical, charges it everywhere. Copying one into the other silently
   shifts every distance that goes through the `X`-stripped path.
8. **Tie-breaking depends on traversal order.** The UPGMA search skips the last
   live element of its linked list and keeps `im` across iterations; the DPs
   use `>=` in some comparisons and `>` in others; `Atracking` picks the first
   equal-scoring end cell and `Atracking_localhom` the last. Every one of those
   decides ties, and ties are common — the duplicate-sequence corpus exists to
   catch exactly this.

One portability note: `sin(PI/n)` in `make_sintbl` is the only libm
transcendental the port depends on. Its argument is exact (`PI` scaled by a
power of two) and Rust's `f64::sin` is the platform libm, so it agrees here;
`tests/golden/sintbl.txt` pins every bit of the resulting table so a platform
where it does not agree fails loudly rather than silently realigning.
