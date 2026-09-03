# `disttbfast-rs` — MAFFT's progressive MSA algorithm, in Rust

A transliteration of MAFFT 7.526's progressive multiple-sequence-alignment
algorithm from C into Rust, written to produce **bit-identical output**, not to
be fast. Where the C and a cleaner formulation disagree, the C wins.

The binary is a drop-in replacement for `libexec/mafft/disttbfast`: same
argument vector, FASTA on stdin, FASTA on stdout, and the same `order` file
written into the working directory. Dropped into a MAFFT installation it makes
plain `mafft`, `fftns`, `fftnsi`, `nwns` and `nwnsi` produce byte-identical
results — including their refinement stages, which are separate binaries that
consume this one's output.

```
cd rust && cargo build --release
./target/release/disttbfast-rs <same args as disttbfast> < in.fa > out.fa
```

## What is ported

The complete FFT-NS-1/2 and NW-NS-1/2 pipeline, end to end:

| stage | C source | Rust |
|---|---|---|
| FASTA reading, type detection, writing | `io.c` (`getnumlen`, `readData_pointer`, `writeData_pointer`) | `io.rs` |
| BLOSUM62 tables | `blosum.c` (`BLOSUMmtx`) | `blosum.rs` |
| score matrix, gap penalties, FFT scales | `constants.c` (`constants`, protein/BLOSUM branch), `miyata.h` | `constants.rs` |
| 6-mer ("ktuples") distance matrix | `disttbfast.c` (`seq_grp`, `makepointtable`, …), `mltaln9.c` (`commonsextet_p`) | `dist.rs` |
| UPGMA guide tree, sequence weights | `mltaln9.c` (`fixed_musclesupg_double_realloc_nobk_halfmtx_memsave`, `counteff_simple_double_nostatic_memsave`, `topolorderz`) | `tree.rs` |
| profile statistics | `tddis.c` (`cpmx_calc_new`), `mltaln9.c` (`st_*GapCount`, `new_*GapCount`, `gapcountf`, `outgapcount`, `gapireru`) | `profile.rs` |
| **group-to-group DP** | `Salignmm.c` (`A__align`, `match_calc`, `Atracking`, `create*result`) | `align.rs` |
| 1-vs-1 DP | `Galign11.c` (`G__align11`, `Atracking`) | `align.rs` |
| **FFT** | `fft.c` (`fft`, `make_sintbl`, `make_bitrev`) | `fft.rs` |
| **FFT anchoring** | `Falign.c` (`Falign`, `seq_vec_5`), `fftFunctions.c` (`getKouho`, `zurasu2`, `alignableReagion`, `blockAlign2`, `calcNaiseki`) | `falign.rs` |
| progressive driver | `disttbfast.c` (`treebase`) | `treebase.rs` |
| MSA-derived distances for the second cycle | `mltaln9.c` (`naivepairscore11`, `naivepairscorefast`, `makeskiptable`) | `score.rs` |
| argument parsing | `disttbfast.c` (`arguments`) | `params.rs` |

Two pieces carry most of the difficulty.

**`A__align`** is the affine-gap DP over two *profiles* with position-specific
gap-open and gap-extend costs. It has two callers taking two different branches
of the C: `treebase` supplies a profile-history slot (`cpmxhist`), which lets a
parent node reuse a child's composition matrix rather than recomputing it;
`Falign` supplies boundary-gap patterns describing what lies just outside the
block being aligned. The cache is not an optimisation detail — using it changes
the summation order, so a port that skipped it would not be bit-exact.

**`Falign`** is the step MAFFT is named for (Katoh, Misawa, Kuma & Miyata,
*NAR* 30:3059–3066, 2002). Each group becomes a complex signal whose real and
imaginary parts are the weighted, normalised polarity and volume of each
column; the cross-correlation of the two signals, computed by FFT, has peaks at
the relative shifts where homologous regions line up. At each of the 20 best
shifts a sliding window scores the two groups column by column and keeps runs
that beat a threshold; those are the candidate anchors. `blockAlign2` then
picks a consistent, monotone subset of them by dynamic programming, and the
groups are cut at the surviving anchors and aligned block by block with the
profile DP above.

## What is *not* ported

`Args::check_supported` refuses these rather than approximating them, so an
unsupported run fails loudly instead of quietly returning a different alignment:

- nucleotide input, matrices other than BLOSUM62, `--fmodel`
- `--add` / `--addprofile`, sub-alignments, `--treein`, external anchors
- `--allowshift` / `--unalignlevel > 0`, and the shift ("warp") DP behind them
- `--legacygapcost`, `--anysymbol`, `--keeplength`
- the memory-saving aligner (`-M`) and the other `alg` values
- **threads.** MAFFT's threaded path is a different computation, not a race:
  it does not agree with the single-thread result even at `--thread 1`. Run
  with `-C 0-0`.

Iterative refinement (`--maxiterate`) is not ported either, but it still works
end-to-end: `dvtditr` is a separate binary that consumes this one's output, so
FFT-NS-i and NW-NS-i come out identical as long as the progressive stage does.
The same is true of the `--parttree` strategies, which run `splittbfast`.

The G-INS-i / L-INS-i strategies do not go through `disttbfast` at all — they
use `tbfast` plus `pairlocalalign` — so this binary does not affect them.

## Verification

Bit-exactness is measured, not argued.

`verify.sh` runs the C `disttbfast` and `disttbfast-rs` on the same inputs with
the same arguments and requires byte-identical stdout *and* `order` file:

```
cd core && make                       # build the C reference
rust/verify.sh core/disttbfast rust/target/release/disttbfast-rs
```

It generates its own corpus: `test/sample`, 40 random subsets of it, sequences
truncated to 1–40 residues (which trip the "too short for a 6-mer" sentinel),
pairs, exact duplicates (ties everywhere in the tree), ragged lengths,
ambiguity codes, pre-existing gaps, random sequences with no phylogenetic
signal, over-long names, lower case, `*` characters, and — for the FFT path —
tandem repeats with varying copy number, shuffled domains and inversions. That
last group matters: it is the only kind of input on which the FFT scan proposes
mutually inconsistent anchors and `blockAlign2` has to prune. Each input is
crossed with 15 argument vectors covering `-E 1/2/3`, terminal-gap on and off,
non-default penalties and offsets, `-X 0.0/0.5`, `-u`, `-n`, `-F`, forced FFT
(`-G`) and a non-default FFT window (`-w 40`). Pass a directory as a third
argument to use your own corpus instead.

`verify_e2e.sh` does the same one level up: it runs the real `mafft` wrapper
against a stock libexec directory and against a copy with `disttbfast`
replaced, across `--nofft` and FFT variants of `--retree 1/2`, `--maxiterate 2`,
`--reorder`, `--clustalout` and `--inputorder`.

```
rust/verify_e2e.sh /path/to/libexec/mafft rust/target/release/disttbfast-rs test/sample
```

Upstream's own `test/script` (via `provolut/run_upstream_tests.sh`) passes 9/9
with the Rust binary in place — the FFT-NS cases through this code, the
G/L-INS-i and parttree cases through the untouched C binaries.

`cargo test` additionally checks, without needing a C build: the derived score
matrix against values dumped from the C `constants()`; the FFT sine table
against a bit-for-bit dump from `core/fft.c`; the bit-reversal permutation; and
three golden alignments.

Current status, on this source tree:

| check | result |
|---|---|
| `verify.sh`, generated corpus (70 inputs × 15 argument sets) | 1050/1050 identical |
| `verify.sh`, 7 larger inputs up to 300 sequences × 15 | 105/105 identical |
| `verify_e2e.sh`, 4 inputs × 11 wrapper strategies | 44/44 identical |
| upstream `test/script` with the Rust binary | 9/9 |
| `cargo test` | 7/7 |

Coverage of the paths that are easy to miss was checked by instrumentation. On
a 60-sequence input, of 87 `A__align` calls 47 produced a profile-history entry
(`icyc + jcyc > 20`), 45 consumed a cached child profile, and 19 hit the
`reuseprofiles` shortcut. On `test/sample` with FFT enabled, of 70 `Falign`
calls 45 split into more than one block (up to 7) and 26 took the 1-vs-1 path;
on the repeat corpus, 36 of 114 calls had anchors pruned by `blockAlign2`.

## Notes for anyone extending this

Five things about the C are load-bearing and easy to "clean up" by accident.
Four of them were found the hard way.

1. **Scratch buffers are never cleared.** `A__align` keeps `static` arrays that
   grow but are not zeroed, and several loops index one element past the current
   sequence length — `ogcp2[lgth2]`, `ogcp1[1]` when `lgth1 == 1`. The Rust
   keeps the same buffers with the same growth rule so the same stale values are
   read back.
2. **Off-by-one reads of the gap tables are intentional.**
   `creategapfreqresult` and `createfgresult` read `gaptable[alen]`, the NUL
   terminator, and rely on it not being `'-'`. `G__align11` reads `seq1[lgth1]`,
   the terminator again, on its last row when terminal gaps are free. The port
   reproduces these explicitly rather than trimming the loops.
3. **`blockAlign2` compares against the wrong variable, and it matters.** The
   loop that fills `maxi` tests candidates against `maxj` — a `static` holding
   the *previous* cell's value, already penalised by `+= penalty`. Writing that
   as a local (the obvious reading) changes which anchors survive, and shows up
   only on inputs with repeats. Its `result1[]` traceback buffer is also read
   one element below where the current call wrote, picking up the previous
   call's value.
4. **`new_FinalGapCount` cannot use its `egappat` argument.** The step after its
   loop has already read the NUL terminator into `gb`, so the test is always
   false; the compiled function is identical to `st_FinalGapCount`. Implementing
   what the argument *looks* like it does would be wrong.
5. **Tie-breaking depends on traversal order.** The UPGMA search skips the last
   live element of its linked list, keeps `im` across iterations, and uses `>=`
   rather than `>` in two of the DP's four comparisons. Every one of those
   decides ties, and ties are common — the duplicate-sequence corpus exists to
   catch exactly this.

One portability note: `sin(PI/n)` in `make_sintbl` is the only libm
transcendental the port depends on. Its argument is exact (`PI` scaled by a
power of two) and Rust's `f64::sin` is the platform libm, so it agrees here;
`tests/golden/sintbl.txt` pins every bit of the resulting table so a platform
where it does not agree fails loudly rather than silently realigning.

The natural next target is `tbfast.c` plus `Lalign11.c` / `pairlocalalign.c`,
which is the L-INS-i path this repository's pipeline actually runs.
