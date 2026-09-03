# `disttbfast-rs` — MAFFT's progressive MSA algorithm, in Rust

A transliteration of MAFFT 7.526's progressive multiple-sequence-alignment
algorithm from C into Rust, written to produce **bit-identical output**, not to
be fast. Where the C and a cleaner formulation disagree, the C wins.

The binary is a drop-in replacement for `libexec/mafft/disttbfast`: same
argument vector, FASTA on stdin, FASTA on stdout, and the same `order` file
written into the working directory. Dropping it into a MAFFT installation makes
`mafft --nofft` produce byte-identical results, refinement stages included.

It is a replacement for the `--nofft` path only. Swapped in wholesale it makes
`fftns` / `fftnsi` / plain `mafft` *fail*, with a message naming the missing
path — deliberately, because the alternative is silently returning a different
alignment. Upstream's own `test/script` shows exactly that: 7/9 pass with the
Rust binary in place, and the two that fail are the two FFT strategies.

```
cd rust && cargo build --release
./target/release/disttbfast-rs <same args as disttbfast> < in.fa > out.fa
```

## What is ported

The whole NW-NS-1 / NW-NS-2 pipeline, end to end:

| stage | C source | Rust |
|---|---|---|
| FASTA reading, type detection, writing | `io.c` (`getnumlen`, `readData_pointer`, `writeData_pointer`) | `io.rs` |
| BLOSUM62 tables | `blosum.c` (`BLOSUMmtx`) | `blosum.rs` |
| score matrix and gap penalties | `constants.c` (`constants`, protein/BLOSUM branch) | `constants.rs` |
| 6-mer ("ktuples") distance matrix | `disttbfast.c` (`seq_grp`, `makepointtable`, …), `mltaln9.c` (`commonsextet_p`) | `dist.rs` |
| UPGMA guide tree, sequence weights | `mltaln9.c` (`fixed_musclesupg_double_realloc_nobk_halfmtx_memsave`, `counteff_simple_double_nostatic_memsave`, `topolorderz`) | `tree.rs` |
| profile statistics | `tddis.c` (`cpmx_calc_new`), `mltaln9.c` (`st_OpeningGapCount`, `st_FinalGapCount`, `gapcountf`, `gapireru`) | `profile.rs` |
| **group-to-group DP** | `Salignmm.c` (`A__align`, `match_calc`, `Atracking`, `create*result`) | `align.rs` |
| 1-vs-1 DP | `Galign11.c` (`G__align11`, `Atracking`) | `align.rs` |
| progressive driver | `disttbfast.c` (`treebase`) | `treebase.rs` |
| MSA-derived distances for the second cycle | `mltaln9.c` (`naivepairscore11`, `naivepairscorefast`, `makeskiptable`) | `score.rs` |
| argument parsing | `disttbfast.c` (`arguments`) | `params.rs` |

`A__align` is the interesting one: it is the affine-gap DP over two *profiles*
with position-specific gap-open and gap-extend costs, plus the profile-history
cache (`cpmxhist`) that lets a parent node reuse a child's composition matrix
instead of recomputing it. That cache is not an optimisation detail — using it
changes the summation order, so a port that skipped it would not be bit-exact.

## What is *not* ported

`Args::check_supported` refuses these rather than approximating them, so an
unsupported run fails loudly instead of quietly returning a different alignment:

- **FFT anchoring** (`-F`/`-G`, i.e. plain `mafft`, `fftns`, `fftnsi`) — `Falign.c`
- nucleotide input, matrices other than BLOSUM62, `--fmodel`
- `--add` / `--addprofile`, sub-alignments, `--treein`, external anchors
- `--allowshift` / `--unalignlevel > 0`, and the shift ("warp") DP behind them
- `--legacygapcost`, `--anysymbol`, `--keeplength`
- the memory-saving aligner (`-M`) and the other `alg` values
- **threads.** MAFFT's threaded path is a different computation, not a race:
  it does not agree with the single-thread result even at `--thread 1`. Run
  with `-C 0-0`.

Iterative refinement (`--maxiterate`) is *not* ported either, but it still works
end-to-end: `dvtditr` is a separate binary that consumes this one's output, so
NW-NS-i comes out identical as long as the progressive stage does.

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
signal, over-long names, lower case, and `*` characters — crossed with 8
argument vectors covering `-E 1/2/3`, terminal-gap on and off, non-default
penalties and offsets, `-X 0.0/0.5`, `-u` and `-n`. Pass a directory as a third
argument to use your own corpus instead.

`verify_e2e.sh` does the same one level up: it runs the real `mafft` wrapper
against a stock libexec directory and against a copy with `disttbfast` replaced,
across `--retree 1/2`, `--maxiterate 2`, `--reorder`, `--clustalout` and
`--inputorder`.

```
rust/verify_e2e.sh /path/to/libexec/mafft rust/target/release/disttbfast-rs test/sample
```

`cargo test` additionally checks the derived score matrix against values dumped
from the C `constants()`, and two golden alignments, without needing a C build.

Current status, on this source tree:

| check | result |
|---|---|
| `verify.sh`, generated corpus (63 inputs × 8 argument sets) | 504/504 identical |
| `verify.sh`, 7 larger inputs up to 300 sequences × 8 | 56/56 identical |
| `verify_e2e.sh`, 4 inputs × 6 wrapper strategies | 24/24 identical |
| `cargo test` | 4/4 |

Coverage of the paths that are easy to miss was checked by instrumenting a run
on a 60-sequence input: of 87 `A__align` calls, 47 produced a profile-history
entry (`icyc + jcyc > 20`), 40 consumed a cached child-0 profile, 5 a cached
child-1 profile, and 19 hit the `reuseprofiles` shortcut.

## Notes for anyone extending this

Three things about the C are load-bearing and easy to "clean up" by accident:

1. **Scratch buffers are never cleared.** `A__align` keeps `static` arrays that
   grow but are not zeroed, and several loops index one element past the current
   sequence length — `ogcp2[lgth2]`, `ogcp1[1]` when `lgth1 == 1`. The Rust
   keeps the same buffers with the same growth rule so the same stale values are
   read back.
2. **Off-by-one reads of the gap tables are intentional.**
   `creategapfreqresult` and `createfgresult` read `gaptable[alen]`, the NUL
   terminator, and rely on it not being `'-'`. The port reproduces that
   explicitly rather than trimming the loop.
3. **Tie-breaking depends on traversal order.** The UPGMA search skips the last
   live element of its linked list, keeps `im` across iterations, and uses `>=`
   rather than `>` in two of the DP's four comparisons. Every one of those
   decides ties, and ties are common — the duplicate-sequence corpus exists to
   catch exactly this.

The natural next targets, in rough order of effort: `Falign.c` (FFT anchoring,
which is what the default `mafft` and `fftns` use), then `tbfast.c` plus
`Lalign11.c` for the L-INS-i path this repository's pipeline actually runs.
