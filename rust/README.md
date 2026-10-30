# mafft-oxide: MAFFT's alignment algorithms, in Rust

[![CI](https://github.com/Provolut/mafft-oxide/actions/workflows/ci.yml/badge.svg)](https://github.com/Provolut/mafft-oxide/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://github.com/Provolut/mafft-oxide/blob/main/rust/LICENSE)

A transliteration of [MAFFT](https://mafft.cbrc.jp/alignment/software/) 7.526's
multiple-sequence-alignment algorithms from C into Rust, written to produce
**bit-identical output**. Where the C and a cleaner formulation disagree, the C
wins.

An independent port of MAFFT by Kazutaka Katoh, not affiliated with or endorsed
by the MAFFT authors.

## What you get

Nine binaries, each a drop-in replacement for its `libexec/mafft/` original —
same argument vector, same input, same output files:

| | |
|---|---|
| `disttbfast-rs` | FFT-NS-1/2, NW-NS-1/2, `--add`, `--merge`, `--treein` (plain `mafft`, `fftns`, `nwns`) |
| `tbfast-rs` | L-INS-1/i, G-INS-1/i, E-INS-1/i (`linsi`, `ginsi`, `einsi`) |
| `dvtditr-rs` | iterative refinement, the `-i` variants |
| `dndpre-rs` | the `hat2` distance matrix for the second progressive cycle |
| `countlen-rs`, `version-rs` | what the driver reads to pick a strategy |
| `replaceu-rs`, `restoreu-rs` | alphabet folding and restoration |
| `f2cl-rs` | the output writer: FASTA, clustal, phylip |

Those are **all nine binaries `mafft --auto` can execute**, so an `--auto` run
executes no C process at all. Protein and nucleotide input throughout.

No dependencies — `Cargo.lock` has one entry, this crate — and
`#![deny(unsafe_code)]`, lifted in exactly one module for an AVX-512 kernel
that dispatches at runtime.

## Use

```
cargo build --release
./target/release/disttbfast-rs <same args as disttbfast> < in.fa > out.fa
```

Or as a conda package, which installs the nine under their C names together
with the upstream driver, as `mafft-oxide`. It ships no `bin/mafft`, so a
`mafft` in the same environment stays stock and the two can be compared:

```
pixi add mafft-oxide
```

## Envelope

`mafft --auto` is covered and byte-identical up to 20,000 sequences. Outside
that the port **refuses rather than approximating**, so an unsupported run
fails instead of quietly returning a different alignment:

- `--allowshift` with iterations — `dvtditr -s`, the specificity scoring path
- `--dash` — `f2cl -d`
- `--parttree`, and inputs past 20,000 sequences — both need `splittbfast`
- `-@` (`--sp`) and the RNA structure strategies

`EXTENDING.md` has the full list, with the measured reason for each.

## Speed

It depends on the shape of the work, in both directions:

| workload | port against the C |
|---|---|
| L-INS-i refinement (ESCOTT's `--auto`, 12 ProteinGym payloads) | **2.32x faster** than `gcc -O3`; 1.38x faster than `-march=native` |
| FFT-NS-2 build path, 790 sequences | **1.14x faster**; **1.59x** with `--thread 8` |

The refinement win is two AVX-512 kernels in `L__align11` plus removed work in
`dvtditr`'s hot loops. The build-path win is in the O(n²) distance matrix that
`--auto`'s second cycle recomputes over the aligned rows: the port reads its gap
skip tables only at a column gapped in both rows, which is the only place their
jump can be non-zero, instead of on every column as the C does.

The build path was 1.37x *slower* until it was measured properly rather than
profiled by function. Everything outside that one matrix is at parity to the
millisecond, so three of the four sites an earlier profile had implicated were
in the half that already matched.

`--thread n` is honoured for the three all-pairs matrices and nothing else.
Each of their cells is a pure function of its own `(i, j)`, so the alignment is
the same for any thread count -- `verify/verify_threads.sh` compares 2,752 runs
across eight `-C` specs against the C's sequential output. The progressive step
stays single-threaded whatever `--threadtb` says, and deliberately:
`core/disttbfast.c:1929` writes the global `alg` from inside a worker with no
lock, and upstream's own comment on the next line ("another thread may have set
it to M") is the reason -- that global picks between two different recurrences,
so above 30,000 columns the C's threaded progressive step has no fixed answer
to be identical to. `EXTENDING.md` has the detail.

The C's kernel is compiled in, so a build carrying it crashes without AVX-512
and a portable build never gets it — which is why the conda-forge 7.526 is the
slowest build in that first row. `bench/OPTIMISATION-LOG.md` carries every
measurement, including the hypotheses that turned out wrong.

## Verification

Bit-exactness is measured, not argued.

| | |
|---|---|
| differential against a compiled `core/` | **12,018 + 1,690 + 2,752** comparisons, all byte-identical |
| unit and property tests, no C needed | **245**, 0 ignored |
| region coverage | **91.37%** |
| upstream's own `test/script` | 8/9 — the ninth is `--allowshift`, refused rather than wrong |

The harnesses live in `verify/`; each script's header says what it stages and
how to run it. `verify_golden.sh` replays the two big ones with no C toolchain
present.

## Extending it

`EXTENDING.md` first. Fifteen things about the C are load-bearing and easy to
"clean up" by accident — scratch buffers that are never zeroed, deliberate
off-by-one reads of the gap tables, a comparison against the wrong variable
that decides which anchors survive. Most were found by a corpus case rather
than by reading the code, and tidying any of them changes alignments.

Two rules above the rest: **never reassociate f64 arithmetic**, and
**transliterate, do not improve**.

## License and citing

BSD 3-Clause. `LICENSE` holds two notices and **both apply**: Provolut GmbH's
for the port, and MAFFT's (Copyright (c) 2009 Kazutaka Katoh) for everything
the port is transliterated from. Keep both with any copy, source or binary; for
binaries BSD requires them in the accompanying documentation. Under clause 3 of
MAFFT's license, Katoh's name may not be used to endorse or promote this port
or anything built on it.

Every alignment method here is MAFFT's. The port takes its programs from the
corresponding `core/*.c` files, the scoring data as MAFFT ships it (the BLOSUM
matrices, the JTT and transmembrane counts, Miyata's volume and polarity
scales, the nucleotide alphabet and the default parameters), and the FFT,
including the sine-table recurrence `core/fft.c` credits to Haruhiko Okumura.
MAFFT's RNA extensions are not part of the port and none of their code is in
it.

**Citing.** A mafft-oxide alignment is a MAFFT alignment, byte for byte, so
cite MAFFT — Katoh & Standley 2013, *Molecular Biology and Evolution*
30:772-780, plus the paper for the strategy used (Katoh et al. 2002, *NAR*
30:3059-3066 for FFT-NS-\*; Katoh et al. 2005, *NAR* 33:511-518 for the
`-INS-i` family). Name mafft-oxide as well only if the implementation matters
to your methods, such as when reporting runtimes.

## More

- `EXTENDING.md` — the refusals, the C's load-bearing quirks, and the rules
- `bench/OPTIMISATION-LOG.md` — every performance measurement
