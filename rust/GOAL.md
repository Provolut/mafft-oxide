# Goal: feature parity with the C, still byte-identical

Definition of done: every `mafft` invocation that works with the stock
`libexec` produces byte-identical output with `disttbfast-rs` and `tbfast-rs`
swapped in, and nothing is refused but `-C > 0` and the RNA strategies.

The hard constraint does not move: **output must not change.** Anything not
ported is refused loudly. `verify.sh`, `verify_tbfast.sh`, `verify_e2e.sh`,
`cargo test` and `test/script` must stay green, and new argument vectors go
into the harnesses with the feature.

Done so far: E-INS-i, every built-in BLOSUM table, JTT/TM, `--fmodel`,
`--legacygapcost`, `-y`/`-T`.
Inventory: issue #2. C quirks to reproduce, not fix: issue #3 and
`README.md`'s "Notes for anyone extending this".

## Work, roughly in order of value

1. **`--text` (`nblosum == -2`)**, the 256-symbol alphabet: `nalphabets` and
   `charsize` become 0x100. Touches every place that assumes 26, so it is the
   one scoring item left and the one that is not cheap. `--aamatrix` (a matrix
   read from a file, `n == -1`) rides along with it.

2. **Input modes**, each with its own harness argument vectors: `--add` /
   `--addprofile` (`disttbfast` gets `-K -I <n>`; the aligning work is in
   `addsingle`, still C, but our binaries are on the path); sub-alignments
   (`-H`); `--treein` (`-U`); external anchors (`-l <file>`); `--keeplength`
   (`-Y`); `-q` (sparse pickup); `-r` (one iteration).

3. **Tree output** (`-t`, `-^`): `..._treeout_memsave` is a second copy of the
   UPGMA with Newick emission bolted on; diff it against the one in `tree.rs`
   and add the string bookkeeping. `-y` and `-T` are done.

4. **`-M` (memsave) and the other `alg` values** (`-&`, `-@`): the linear-space
   aligner in `MSalignmm.c` / `partSalignmm.c`. Large, low demand.

5. **10-mer nucleotide tuples** (`--tuplesize 10`): `makepointtable_nuc_dectet`,
   `D10LENFAC*`. Only `--parttree` uses it; do it with 8 if at all.

6. **Optional, separate binaries.** `dvtditr` (`--maxiterate`), `splittbfast`
   (`--parttree`), `addsingle` (`--add`), `makedirectionlist`
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
