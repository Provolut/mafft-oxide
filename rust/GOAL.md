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

Size each item from the C, and measure how it is *reached* -- an option
nobody types can still be the default path, as item 3 turned out to be. Strip
the dead `#if`s before sizing, and check which switches are really 0.

1. ~~**External anchors (`-l`)**~~ done: `anchors.rs`,
   `Falign_givenanchors`. 1376 comparisons; two C quirks reproduced and one
   input refused, all in issue #4.

2. ~~**The `compacttree*` builders**~~ done: `compacttree.rs`. All four
   keywords byte-identical over 3096 comparisons, `infile.tree` included.

3. **The `alg` values.** `-M` done: `msalign.rs`, `Falign_udpari_long` as a
   mode of `falign.rs`, and the C's 30000-column switch in `treebase`. 688
   comparisons and `mafft --memsave` end to end; peak RSS on four
   31000-column sequences 6.4GB -> 55.9MB (C: 32.7MB). `-&` done: 344
   comparisons. Two upstream defects found here, in issue #4.
   Left: **`-@` (`--sp`)**, `D__align`: the generalized-affine DP, carrying a
   per-column distribution of gap *lengths* rather than two running scores.
   Half done: `dalign.rs` has the `Gaplen` layer (`gaplencount`,
   `fillgaplen`, the `compactx` helpers, all three live `calcpfac`), bit-exact
   over 203 fuzzed group pairs -- `verify_gaplen.sh` links a harness against a
   built `core/` and compares f64 bit patterns.
   Left: the DP (~640 live lines) plus `Dalignmm.c`'s `match_calc` and
   `Atracking`. `gen_dharness.c` calls `D__align` directly for a fast oracle.
   Three notes: columns never alias (`(Gaplen *)1` is a barrier stopping the
   partial free and duplicate loops); `copygaplencompactx` writes into an
   existing destination without checking its length, so each call site must
   size it with `duplicategaplencompactx` first; and `WMCHECK` is **1**, so
   the tail of `D__align` is live -- it recomputes the score naively and, on
   a mismatch over 0.01, dumps both groups to *stdout* and `exit(1)`s.

4. **Optional, separate binaries.** `dvtditr` (`--maxiterate`),
   `splittbfast`, `addsingle`, `makedirectionlist`. All work today on our
   byte-identical output; porting them is scope growth, not parity.

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
