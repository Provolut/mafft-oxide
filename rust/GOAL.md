# Goal: feature parity with the C, still byte-identical

Definition of done: every `mafft` invocation that works with the stock
`libexec` produces byte-identical output with `disttbfast-rs` and `tbfast-rs`
swapped in, and nothing is refused but `-C > 0` and the RNA strategies.

The hard constraint does not move: **output must not change.** Anything not
ported is refused loudly. `verify.sh`, `verify_tbfast.sh`, `verify_e2e.sh`,
`cargo test` and `test/script` must stay green, and new argument vectors go
into the harnesses with the feature.

Inventory: issue #2. C quirks to reproduce, not fix: issue #3 and
`README.md`'s "Notes for anyone extending this".

## Work, roughly in order of value

1. **E-INS-i / E-INS-1 (`--genafpair`).** The largest missing strategy. Port
   `genL__align11` from `core/genalign11.c` — the generalised affine-gap local
   aligner — as a third arm of `pairalign.rs` beside `L__align11` and
   `G__align11`. The wrapper selects it with `-N -Z` in `tbfast`'s pairwise
   stage; `bin/tbfast.rs` currently refuses anything but `-L`/`-A`. Self
   contained: one 1-vs-1 DP, and `verify_tbfast.sh` already has the shape.

2. **Scoring configurations.** Cheap, mostly `constants.rs`:
   - other BLOSUM tables (`-b 45`, `-b 80`); `blosum.rs` already parses the
     format, the refusal is conservative rather than structural;
   - `--jtt` / `--tm` (`JTTmtx` in `core/JTT.c`, one flag apart), same shape as
     the DNA branch: a PAM matrix raised to a power, logged, rescaled;
   - `--fmodel` (`-a`): use `calcfreq`/`calcfreq_nuc` over the input instead of
     the table's frequencies;
   - `--text` (`nblosum == -2`): the 256-symbol alphabet, `nalphabets = 0x100`,
     `charsize = 0x100`. Touches every place that assumes 26.

3. **`--legacygapcost` (`-L`).** The pre-7.110 gap scoring, a separate branch
   inside `A__align`.

4. **Input modes**, each with its own harness argument vectors: `--add` /
   `--addprofile` (`disttbfast` gets `-K -I <n>`; the aligning work is in
   `addsingle`, still C, but our binaries are on the path); sub-alignments
   (`-H`); `--treein` (`-U`); external anchors (`-l <file>`); `--keeplength`
   (`-Y`); `-q` (sparse pickup); `-r` (one iteration).

5. **Output-only modes** (`-y` distances, `-t`/`-^` tree, `-T` noalign): the
   runs that skip the alignment. Small, and they unblock `--treeout`.

6. **`-M` (memsave) and the other `alg` values** (`-&`, `-@`): the linear-space
   aligner in `MSalignmm.c` / `partSalignmm.c`. Large, low demand.

7. **10-mer nucleotide tuples** (`--tuplesize 10`): `makepointtable_nuc_dectet`,
   `D10LENFAC*`. Only `--parttree` uses it; do it with 8 if at all.

8. **Optional, separate binaries.** `dvtditr` (`--maxiterate`), `splittbfast`
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
