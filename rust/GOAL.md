# Goal: feature parity with the C, still byte-identical

Done: every `mafft` invocation that works with the stock `libexec` produces
byte-identical output with `disttbfast-rs` and `tbfast-rs` swapped in, and
nothing is refused but `-C > 0` and the RNA strategies.

The constraint does not move: **output must not change.** Anything not ported
is refused loudly. `verify.sh`, `verify_tbfast.sh`, `verify_e2e.sh`,
`cargo test` and `test/script` must stay green, and new argument vectors go
into the harnesses with the feature.

Done: E-INS-i, every built-in BLOSUM table, JTT/TM, `--fmodel`,
`--legacygapcost`, `-y`/`-T`/`-t`/`-^`/`-U`/`-q`, 10-mer nucleotide tuples.
Inventory: issue #2. C quirks to reproduce, not fix: issue #3 and
`README.md`'s "Notes for anyone extending this".

## Work

Everything cheap is done. What is left is bigger than anything above it in the
history, and each item is sized here from reading the C, not guessed.

1. **`--add` / `--addprofile`** (`-I <n>`, `-K`), with `--keeplength` (`-Y`)
   riding along. The most-used missing mode. `mergeoralign` gains `'1'`, `'2'`
   and `'w'`, decided by `topolorderz` over the tree we already build, so no
   new tree representation -- but `treebase` then needs `findcommongaps`,
   `commongappick`, `adjustgapmap`, `restorecommongaps`, `findnewgaps` and
   `insertnewgaps` from `addfunctions.c`, and `insertnewgaps` calls
   `profilealignment`, a further aligner. Budget several hundred lines plus one
   more DP.

2. **Sub-alignments (`-H`).** Needs a *different* guide tree:
   `fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained` keeps each
   sub-alignment monophyletic and stores full member lists, where `tree.rs`
   stores the memsave single-member form. So it is a second tree function and a
   second `topol` representation, plus `readsubalignmentstable` and the
   `mergeoralign == 'n'` skip in `treebase`.

3. **`--text` (`nblosum == -2`)** and `--aamatrix` (`n == -1`, a matrix read
   from a file). The 256-symbol alphabet: `nalphabets` and `charsize` become
   0x100, so `NALPHABETS` has to stop being a compile-time constant and every
   array and loop bound derived from it follows. Mechanical, wide, and every
   site is a chance to break bit-exactness -- but the harnesses would catch it.

4. **`-r` (one iteration)**: `dooneiteration`, a second copy of `treebase` with
   its own bookkeeping.

5. **External anchors (`-l <file>`)**: `readexternalanchors`, `indexanchors`
   and the `Falign_givenanchors` path. Note the C's `case 'l'` falls through
   into `case 'G'`, so `-l` also forces the FFT on.

6. **The generated guide trees `--treein` can ask for** -- `shuffle`,
   `pileup`, `auto`, `compact`, `stepadd`, `youngestlinkage`, `nodepair` --
   each a separate construction (`createchain`, `varpairscore`, `compacttree`).
   Refused by name today.

7. **`-M` (memsave) and the other `alg` values** (`-&`, `-@`): the linear-space
   aligner in `MSalignmm.c` / `partSalignmm.c`. Largest, lowest demand.

8. **Optional, separate binaries.** `dvtditr` (`--maxiterate`), `splittbfast`
   (`--parttree`), `addsingle` (`--addfragments`), `makedirectionlist`
   (`--adjustdirection`). All work today on our byte-identical output; porting
   them is scope growth, not parity.

## Threads

Not required for parity, and not barred by the output rule either: the C's
`disttbfast` and `tbfast` are thread-invariant (39 inputs x 3 repeats, `-C 4`
vs `-C 0`, all four output files byte-identical). The `mafft --thread`
divergence comes from `dvtditr`, which we do not port. So parallelising the
all-pairs pairwise stage -- embarrassingly parallel, and where the 1.2-1.4x gap
lives -- is a legitimate option if the harnesses stay green under repeats.
Accepting `-C > 0` without implementing it is not.

## Explicitly out of scope

- **`--qinsi` / `--xinsi`**, which shell out to McCaskill/CONTRAfold.

## Method

Transliterate; do not improve. When output diverges, instrument a scratch copy
of `core/` to dump intermediates as f64 bit patterns (`%016llx`), dump the same
from Rust behind an env var, and bisect. Every mismatch so far was found that
way, never by reading.

Performance is not a parity goal (~1.2-1.4x the C today). Do not trade
evaluation order for speed.
