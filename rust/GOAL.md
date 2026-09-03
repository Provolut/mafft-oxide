# Goal: feature parity with the C, still byte-identical

Done: every `mafft` invocation that works with the stock `libexec` produces
byte-identical output with `disttbfast-rs` and `tbfast-rs` swapped in, and
nothing is refused but `-C > 0` and the RNA strategies.

The constraint does not move: **output must not change.** Anything not ported
is refused loudly. `verify.sh`, `verify_tbfast.sh`, `verify_e2e.sh`,
`cargo test` and `test/script` must stay green, and new argument vectors go
into the harnesses with the feature.

Done: E-INS-i, every built-in BLOSUM table, JTT/TM, `--fmodel`,
`--legacygapcost`, `-y`/`-T`/`-t`/`-^`/`-U`/`-q`, 10-mer nucleotide tuples,
`--add`/`--addprofile`/`--keeplength`/`--mapout`, `--allowshift`, `--aamatrix`, `-r`.
Inventory: issue #2. C quirks to reproduce, not fix: issue #3 and
`README.md`'s "Notes for anyone extending this".

## Work

Each item is sized here from reading the C, not guessed. Size the refusal too,
before porting: `--allowshift` cost nothing -- the warp DP was already in
`warp.rs` and correct in all three aligners, and only a stale up-front
`exit(1)` in `disttbfast.rs` kept it unreachable. And check both copies:
`tbfast`'s `--add` is not `disttbfast`'s (it skips the `samemembern` pass and
divides the added sequences' weights by 100), which cost a day.

1. **Sub-alignments (`-H`).** Needs a *different* guide tree:
   `..._treeout_constrained` keeps each sub-alignment monophyletic and stores
   full member lists, where `tree.rs` stores the memsave single-member form.
   So: a second tree function and a second `topol` representation, plus
   `readsubalignmentstable`.

2. **`--text` / `--textmatrix` (`nblosum == -2`).** The 256-symbol alphabet:
   `nalphabets` and `charsize` become 0x100, so `NALPHABETS` stops being a
   compile-time constant and every bound derived from it follows. Mechanical,
   wide, and every site is a chance to break bit-exactness -- but the
   harnesses would catch it. (`--aamatrix`, `nblosum == -1`, is done: it keeps
   the 26-letter alphabet and only swaps the table.)

3. **External anchors (`-l <file>`)**: `readexternalanchors`, `indexanchors`,
   `Falign_givenanchors`. The C's `case 'l'` falls through into `case 'G'`, so
   `-l` also forces the FFT on.

4. **The generated guide trees `--treein` can ask for** -- `shuffle`,
   `pileup`, `auto`, `compact`, `stepadd`, `youngestlinkage`, `nodepair`, each
   a separate construction. Refused by name today.

5. **`-M` (memsave) and the other `alg` values** (`-&`, `-@`): the linear-space
   aligner in `MSalignmm.c` / `partSalignmm.c`. Largest, lowest demand.

6. **Optional, separate binaries.** `dvtditr` (`--maxiterate`), `splittbfast`,
   `addsingle`, `makedirectionlist`. All work today on our byte-identical
   output; porting them is scope growth, not parity.

## Threads

Not required for parity, and not barred by the output rule: the C's
`disttbfast` and `tbfast` are thread-invariant (39 inputs x 3 repeats, `-C 4`
vs `-C 0`, byte-identical). The `mafft --thread` divergence comes from
`dvtditr`, which we do not port. Parallelising the all-pairs pairwise stage is
therefore legitimate if the harnesses stay green under repeats. Accepting
`-C > 0` without implementing it is not.

## Explicitly out of scope

- **`--qinsi` / `--xinsi`**, which shell out to McCaskill/CONTRAfold.

## Method

Transliterate; do not improve. When output diverges, instrument a scratch copy
of `core/` to dump intermediates as f64 bit patterns, dump the same from Rust
behind an env var, and bisect. Every mismatch so far was found that way, never
by reading.

Performance is not a parity goal (~1.2-1.4x the C today). Do not trade
evaluation order for speed.
