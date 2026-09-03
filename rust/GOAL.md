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
`--add`/`--addprofile`/`--keeplength`/`--mapout`, `--allowshift`, `--aamatrix`, `-r`, `--merge` (`-H`), `--text`/`--textmatrix`, `--treein pileup`/`shuffle`, external anchors (`-l`).
Inventory: issue #2. C quirks to reproduce, not fix: issue #3 and
`README.md`'s "Notes for anyone extending this".

## Work

Size each item from the C, and measure how it is *reached* -- an option nobody
types can still be the default path, as item 3 turned out to be. Size the
refusal too: `--allowshift` cost nothing, only a stale `exit(1)`. And check
both copies: `tbfast`'s `--add` is not `disttbfast`'s.

1. ~~**External anchors (`-l`)**~~ done: `anchors.rs` plus
   `Falign_givenanchors`. 1376 comparisons byte-identical; `verify.sh` block,
   files from `gen_anchors.py`. Two C quirks reproduced and one input refused,
   all in issue #4.

2. ~~**The `compacttree*` builders**~~ done: `compacttree.rs` has
   `preparepartmtx`, the two nearest-neighbour distance passes (k-mer and
   MSA), `compacttreegivendist` with `reformattree`, and
   `compacttree_memsaveselectable` for `youngestlinkage`. All four keywords --
   `compact`, `very compact`, `stepadd`, `youngestlinkage` -- byte-identical
   over 3096 comparisons including `infile.tree`.

3. **`-M` (memsave) and the other `alg` values** (`-&`, `-@`): the linear-space
   aligner in `MSalignmm.c` / `partSalignmm.c`, ~5600 live lines. An earlier
   draft called this "lowest demand". Measuring says the opposite: the C
   switches to it *by itself* whenever a row passes 30000
   columns (`disttbfast.c:2836`, `tbfast.c:1096`), `alg` is a global and never
   switches back, and no `--memsave` is needed. It is not an equivalent
   optimisation either -- forcing it with `-M -B` changes the alignment on 17
   of the 86 corpus inputs. The port now refuses at that threshold instead of
   silently continuing on `A__align`; `-N` (`--nomemsave`) suppresses the
   switch in the C and so also the refusal. Until this is ported, long input
   is a refusal, so this is the item that decides the goal.

4. **Optional, separate binaries.** `dvtditr` (`--maxiterate`), `splittbfast`,
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
