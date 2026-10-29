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
`recipe/` on `main` from a tagged release tarball pinned by sha256.
