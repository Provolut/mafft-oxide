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
`--add`/`--addprofile`/`--keeplength`/`--mapout`, `--allowshift`, `--aamatrix`, `-r`, `--merge` (`-H`), `--text`/`--textmatrix`, `--treein pileup`/`shuffle`.
Inventory: issue #2. C quirks to reproduce, not fix: issue #3 and
`README.md`'s "Notes for anyone extending this".

## Work

Each item is sized here from reading the C, not guessed. Size the refusal too,
before porting: `--allowshift` cost nothing -- the warp DP was already in
`warp.rs` and correct in all three aligners, and only a stale up-front
`exit(1)` in `disttbfast.rs` kept it unreachable. And check both copies:
`tbfast`'s `--add` is not `disttbfast`'s (it skips the `samemembern` pass and
divides the added sequences' weights by 100), which cost a day.

1. **External anchors (`-l <file>`)**: `readexternalanchors`, `indexanchors`,
   `Falign_givenanchors` (~670 lines). Blocked, not just big: `anchcomp` and
   `anchscorecomp` break ties on `(int)((void *)p - (void *)q)`, the addresses
   of the two elements *inside the array qsort is permuting*. Matching that
   means reimplementing glibc's qsort, not just the comparator. Decide whether
   to accept a documented divergence before starting.

2. **The `compacttree*` builders**, reached by `--treein very compact`,
   `stepadd` and `youngestlinkage`, and by `compact` only above its RAM
   budget -- below it the C falls back to the ordinary tree, which is ported
   and covered. `compacttree_memsaveselectable` is 637 lines and brings its
   own distance pipeline (`partmtx`, `mindist`, `selfscore`, `skiptable`).
   (`pileup`, `shuffle`, `nodepair` done; `auto`/`test` done but refusing the
   inputs where the C's `varpairscore` reads past its `pickseq` row.)

3. **`-M` (memsave) and the other `alg` values** (`-&`, `-@`): the linear-space
   aligner in `MSalignmm.c` / `partSalignmm.c`. Largest, lowest demand.

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
