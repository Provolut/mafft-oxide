# mafft-oxide conda channel

A conda channel, not a source tree. Served from the `gh-pages` branch of
[Provolut/mafft-oxide](https://github.com/Provolut/mafft-oxide); the source
lives on `main`.

```
pixi workspace channel add https://provolut.github.io/mafft-oxide/
pixi add mafft-oxide
```

This provides `mafft-oxide`, not `mafft`, so a conda-forge `mafft` in the same
environment is untouched and the two can be compared against each other.

Only reproducible builds are served here: every artifact is built by
`recipe/` on `main` from a fixed archive pinned by sha256.

The channel is append-only. A `pixi.lock` records a URL *and* a sha256, so
removing or rewriting a build 404s every locked install; new builds are added
alongside the old ones instead.

Current build: `0.3.0-h1250391_3`. Build 2 is still served for anything locked
against it, and produces the same alignments to the byte -- it is the speed
that differs. On a 790-sequence `--auto` run against the stock C: build 2 at
1.35x of the C's time, build 3 at 0.87x, or 0.64x with `--thread 8`. The two
were cut from different commits, which `recipe/meta.yaml` records.
