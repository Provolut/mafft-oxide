# MAFFT 7.526, Provolut fork

[![CI](https://github.com/Provolut/mafft-oxide/actions/workflows/ci.yml/badge.svg)](https://github.com/Provolut/mafft-oxide/actions/workflows/ci.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](LICENSE)

This repository is **MAFFT 7.526** by Kazutaka Katoh
(https://mafft.cbrc.jp/alignment/software/), imported unmodified in commit
`9660455`, plus work by Provolut GmbH to make it faster without changing a
single byte of its output:

- **faster C**: performance changes to `core/`, measured and verified
  byte-identical against the unmodified build; see `provolut/README.md`.
- **mafft-oxide**: a Rust port of the aligner binaries, verified
  byte-identical against the C; see `rust/README.md`.

This is an independent fork. It is not affiliated with or endorsed by the
MAFFT authors, and all credit for the alignment method is theirs. If you use
it in published work, cite MAFFT (see "CITING" below), not this fork.

The rest of this file is upstream's README, with notes where the fork differs.

# MAFFT version 7.526
Multiple sequence alignment program
<br>
https://mafft.cbrc.jp/alignment/software/

## COMPILE
     % cd core
     % make clean
     % make
     % cd ..

**This fork does not ship the compiled programs** -- `make` above is not
optional, it is how `binaries/` gets populated. They were committed once, and
the cost was daily: 61 executables and 15 MB, every one showing as modified
after any build so `git status` was noise, and nothing checked that the
committed copies had been built from the committed source. They had not.

If you have the `./extensions` directory, which is for RNA alignments,

     % cd extensions
     % make clean
     % make
     % cd ..


## INSTALL (select a or b below)
###  a. Install to /usr/local/ using root account
     # cd core
     # make install
     # cd ..

If you have the `./extensions` directory,

     # cd extensions 
     # make install
     # cd ..

By this procedure (a), programs are installed into `/usr/local/bin/`. Some binaries, which are not directly used by a user, are installed into `/usr/local/libexec/mafft/`.

If the MAFFT_BINARIES environment variable is set to `/somewhare/else/`, the binaries in this directory are used, instead of those in `/usr/local/libexec/mafft/`.

### b. Install to non-default location (root account is not necessary)
     % cd core/
          Edit the first line of Makefile 
          From:
          PREFIX = /usr/local
          To:
          PREFIX = /home/your_home/somewhere

          Edit the third line of Makefile 
          From:
          BINDIR = $(PREFIX)/bin
          To:
          BINDIR = /home/your_home/bin 
                   (or elsewhere in your command-search path)
     % make clean
     % make
     % make install

If you have the `./extensions` directory,

     % cd ../extensions/
          Edit the first line of Makefile 
          From:
          PREFIX = /usr/local
          To:
          PREFIX = /home/your_home/somewhere
     % make clean
     % make
     % make install

The `MAFFT_BINARIES` environment variable *must not be* set.

If the `MAFFT_BINARIES` environment variable is set to `/somewhare/else/`, it overrides the setting of `PREFIX` (`/home/your_home/somewhere/` in the above example) in Makefile.

## CHECK
     % cd test
     % rehash                                                   # if necessary
     % mafft sample > test.fftns2                               # FFT-NS-2
     % mafft --maxiterate 100  sample > test.fftnsi             # FFT-NS-i
     % mafft --globalpair sample > test.gins1                   # G-INS-1 
     % mafft --globalpair --maxiterate 100  sample > test.ginsi # G-INS-i 
     % mafft --localpair sample > test.lins1                    # L-INS-1 
     % mafft --localpair --maxiterate 100  sample > test.linsi  # L-INS-i 
     % diff test.fftns2 sample.fftns2
     % diff test.fftnsi sample.fftnsi
     % diff test.gins1 sample.gins1
     % diff test.ginsi sample.ginsi
     % diff test.lins1 sample.lins1

If you have the `./extensions` directory,

     % mafft-qinsi samplerna > test.qinsi                       # Q-INS-i
     % mafft-xinsi samplerna > test.xinsi                       # X-INS-i
     % diff test.qinsi samplerna.qinsi
     % diff test.xinsi samplerna.xinsi

If you use the multithread version, the results of iterative refinement methods (`*-*-i`) are not always identical.  So try this test in the single-thread mode (`--thread 0`).


## INPUT FORMAT
Fasta format.

The type of input sequences (nucleotide or amino acid) is automatically recognized based on the frequency of A, T, G, C, U and N.


##  USAGE
     % /usr/local/bin/mafft input > output

See also https://mafft.cbrc.jp/alignment/software/


## UNINSTALL
     # rm -r /usr/local/libexec/mafft
     # rm /usr/local/bin/mafft
     # rm /usr/local/bin/fftns
     # rm /usr/local/bin/fftnsi
     # rm /usr/local/bin/nwns
     # rm /usr/local/bin/nwnsi
     # rm /usr/local/bin/linsi
     # rm /usr/local/bin/ginsi
     # rm /usr/local/bin/mafft-*
     # rm /usr/local/share/man/man1/mafft*


## CITING
The MAFFT homepage lists these as the papers describing it; for version 7:

- Katoh, Standley 2013. MAFFT multiple sequence alignment software version 7:
  improvements in performance and usability. *Molecular Biology and Evolution*
  30:772-780.

and, for the strategies used:

- Katoh, Misawa, Kuma, Miyata 2002. MAFFT: a novel method for rapid multiple
  sequence alignment based on fast Fourier transform. *Nucleic Acids Res.*
  30:3059-3066. (FFT-NS-1, FFT-NS-2, FFT-NS-i)
- Katoh, Kuma, Toh, Miyata 2005. MAFFT version 5: improvement in accuracy of
  multiple sequence alignment. *Nucleic Acids Res.* 33:511-518. (G-INS-i,
  L-INS-i, E-INS-i)

The homepage lists further papers for the other strategies.


## LICENSE
`./LICENSE` is the BSD 3-Clause license, with two copyright holders:

- **Kazutaka Katoh**, for MAFFT itself: everything imported in `9660455`.
  That includes the wrapper scripts in `scripts/` and `core/` (`mafft*`,
  `*.rb`, `*_premafft.pl`), which are built from upstream's `.tmpl` files.
- **Provolut GmbH**, for this fork's own work: `rust/` (mafft-oxide),
  `provolut/`, `.githooks/`, `scripts/check-encoding.sh`, and the changes to
  `core/` made after `9660455`. A changed upstream file is under both.

`./LICENSE` uses the standard BSD 3-Clause wording so that tools, GitHub
included, recognise it. MAFFT's own notice, whose wording differs slightly but
whose three conditions and disclaimer are the same, is kept verbatim in
`./licenses/MAFFT` (upstream's `license` file), and it is what applies to
Katoh's code.

**Exception: `./extensions` is not BSD.** The RNA extensions (Vienna RNA
package, MXSCARNA, ProbConsRNA) each have their own terms, in
`./licenses/MAFFT-extensions` (upstream's `license.extensions` file). In
particular, the Vienna RNA package may not be redistributed for a fee, and
asks commercial users to contact its authors. Neither the Rust port nor the
fork's C changes use any of it.

`rust/LICENSE` is the Rust port's own copy of both notices, Provolut's and
MAFFT's verbatim, so the crate carries them when it is distributed apart from
this repository.
