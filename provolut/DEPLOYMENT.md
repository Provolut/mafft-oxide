# Getting a faster MAFFT to where it matters

The speedups in this repo are worth nothing until the binary the pipeline
actually executes is this one. Today it is not: `mafft` comes from conda-forge,
resolved through `$PATH` in the `pml` environment, and
`pml.bioseq.gemme.pipeline._run_aligner` shells out to whatever that finds.

That is the whole deployment question. It is a real decision with tradeoffs,
not a packaging detail, so it is written down rather than assumed.

## What is on the table

Measured on a 40-subset corpus drawn from a real 4,927-sequence alignment:

| build | corpus | vs production |
|---|---|---|
| conda-forge mafft 7.526 (production today) | ~223s (extrapolated) | 1.00x |
| this source, plain `gcc -O2`, no patches | 152.5s | ~1.5x |
| + the AVX-512 kernel and the accumulator | 103.6s | ~2.2x |

Two separate things are bundled there, and they can ship separately:

- **Roughly a third of the win needs no patch at all.** Simply compiling v7.526
  with a current GCC beats the conda-forge binary. That is a packaging change
  with no code risk.
- **The rest is the patches in this repo**, which carry the risk of any
  hand-vectorised numerical kernel, mitigated by the verification in the README.

## The options, and what each costs

**1. Do nothing.** Defensible. The trace cache already takes a repeat round's
scoring from ~47s to ~1.2s, so this only helps the first score of a parent
nobody has scored before. If cold-path latency is not a customer-visible
problem, none of this is urgent.

**2. Vendor a built binary.** Fastest to ship, worst to live with: a binary blob
in the tree, built for one architecture, and `-march=native` is not portable
across the machines that might run it. Rejected unless someone wants it argued
again.

**3. Build from source during environment setup.** `scripts/setup-base-env.sh` grows a step that
fetches v7.526, applies this repo's patches and installs into the conda prefix
ahead of the conda-forge binary. Honest about what it is: the project now
maintains a patched third-party aligner, which means tracking upstream releases
and re-verifying each one. Needs a compiler on every machine that sets up the
environment.

**4. Build it once, distribute as a conda package.** The clean version of (3).
A `provolut-mafft` package built in CI, pinned in `environment.yml`, so
workstations and the compute host get an identical binary and nobody needs a
toolchain. The recipe itself is small -- roughly forty lines pointing at a
pinned tag of this repo, a build script that runs `make` in `core/`, and the
toolchain as a build dependency. Three constraints make it more than that, and
they are the reason this is not simply the obvious answer. See "What packaging
actually costs" below.

**5. Upstream the patches.** The right long-term home. They are bit-exact
against upstream's own test suite, so they are contributable, and if accepted
everyone's MAFFT gets faster and the fork disappears. Slowest path, and not
mutually exclusive with any of the above.

## Why a hand-installed binary does not survive

Worth stating before the options are weighed, because it rules out the
tempting shortcut. `environment.yml` declares `- mafft` as a conda dependency,
and `scripts/setup-conda-env.sh` handles an existing environment two ways:

- **Default**: `conda env update -n pml -f environment.yml`. The declared
  dependency is reinstalled, overwriting a replacement dropped into
  `$CONDA_PREFIX/bin`.
- **`RECREATE_ENV=true`**: `conda remove -n pml --all -y`. The prefix is
  deleted outright.

So copying a faster binary over the conda one is undone by `scripts/setup-base-env.sh`, by
`scripts/setup-conda-env.sh`, or by a colleague refreshing their environment -- and it is
undone *silently*. Results stay correct; scoring just quietly takes twice as
long again, with no error and no log line. Installing outside the prefix and
ahead of it on `PATH` survives, but only on the one machine that was set up
that way, which is a benchmarking arrangement rather than a deployment.

## What packaging actually costs

**1. Strict channel priority would ignore the package.** `scripts/setup-conda-env.sh` sets
`CONDA_CHANNEL_PRIORITY=strict` and `environment.yml` lists conda-forge first.
Under strict priority conda takes a package from the highest-priority channel
that has it and does not look lower -- so a package *named* `mafft` in a private
channel below conda-forge is never installed, with no error. Raising the private
channel above conda-forge is the wrong fix: this file's own header warns that a
mixed solve splits the native OpenMP stack and leaves compiled extensions
unloadable (CHANGELOG v1.0.2). The fix is to name it `provolut-mafft`, have it
provide `bin/mafft`, drop `- mafft` from `environment.yml` and depend on the new
package instead. No collision, no priority games.

**2. `-march=native` cannot ship, and most of the speedup is behind it.** The
AVX-512 kernel is guarded by `#if defined(__AVX512F__) && defined(__AVX512VL__)`,
resolved at *compile* time. A portable package built for a generic baseline
compiles the scalar fallback and keeps only the ~1.5x that comes from building
with a current compiler. The choice:

| build target | runs on | keeps the kernel |
|---|---|---|
| baseline (x86-64-v2/v3) | anything | no -- most of the win is gone |
| require AVX-512 (x86-64-v4) | recent hardware only | yes |
| runtime dispatch | anything | yes |

Runtime dispatch -- the kernel compiled twice into separate translation units,
selected by CPUID at startup -- is what numerical libraries do and is the only
option that is both portable and fast. It is also real additional work on top
of the packaging, and it is the item most likely to double the estimate.

**3. macOS.** `scripts/setup-conda-env.sh` applies `environment.yml` on non-Linux hosts too
(it only strips the Linux-only `libgomp` pin), and a conda env file cannot carry
platform selectors -- this repo already documents that limitation. A Linux-only
package therefore fails the macOS solve outright. Either build for osx as well,
or find a way to make the dependency conditional that `environment.yml` cannot
currently express.

Ongoing, beyond the build: somewhere to host the channel (anaconda.org private,
or any indexed directory served over HTTP), CI to build and publish on tag, and
a named owner to re-verify against each upstream MAFFT release -- the patches
are against v7.526 specifically.

## What has to be decided, by a person

- **Is cold-path scoring latency actually a problem?** If the answer is no,
  option 1 is correct and the rest of this document is premature.
- **If it is: (4) or (5)?** (3) is the tempting middle and the worst of both --
  every machine builds its own copy, and drift between them is invisible.
- **Is a portable package worth runtime dispatch?** Without it, packaging
  keeps the compiler win and loses the kernel -- which may be the right
  trade, but it should be chosen rather than discovered.
- **Who re-verifies on the next MAFFT release?** A patched aligner that nobody
  re-checks is a correctness risk pretending to be a speedup. The verification
  is scripted (`provolut/run_upstream_tests.sh` plus the corpus diff), so this
  is a question of ownership, not effort.

## One thing that is not optional

Whatever route is taken, the deployed binary must be verified on the machine
that will run it. The kernel is guarded by
`#if defined(__AVX512F__) && defined(__AVX512VL__)` so the source compiles
anywhere, but that guard resolves at compile time, not run time: a binary built
with `-march=native` on a machine with AVX-512 will not run on one without it.


# Measured, 2026-10-28: the components are already enough

A question this document left open -- whether the Rust port can serve a real
`mafft` invocation -- turns out to be answerable by measurement rather than by
inventory, and the answer is yes for the only invocation this pipeline makes.

## What the pipeline calls

`mafft --auto --preservecase <input.fasta>`, from
`src/pml/bioseq/gemme/pipeline.py:366` and `src/pml/bioseq/msa.py:687`.

## What `--auto` selects, and what it then runs

`core/mafft:1290` picks the strategy from sequence count and length. Each
component was wrapped in a logging shim and the real driver run, so this is
what fired rather than what the script could reach:

```
branch                                    components invoked
nlen<3000 & nseq<100   (local, it=1000)   countlen replaceu tbfast dvtditr restoreu f2cl version
nlen<1000 & nseq<200   (local, it=2)      same seven
nlen<10000 & nseq<500  (ktuples, cyc=2)   the above, plus dndpre and disttbfast
```

**All of those are in `mafft-oxide` 0.3.0**, and have been since 0.2.0. `pairlocalalign` is never
invoked: on the local branch the pairwise stage runs inside `tbfast -L`, which
is why the port's `sample.linsi` profiling case is a `tbfast` invocation and
not a `pairlocalalign` one. A component inventory suggests 26 are missing; for
this invocation the number that matters is zero.

What does still fall through to C: the `parttree` branches, which need
`splittbfast` (unported). Those are reached at roughly 200k sequences or very
long input -- nothing this pipeline feeds.

## Byte-identical through the real driver

The C driver with all seven (or nine) components replaced by their Rust
counterparts, against pure C, on `mafft --auto --preservecase`:

```
25 inputs, 25 byte-identical, 0 differing
  test/sample (36 protein seqs), aa.fa, dup.fa, tandem.fa
  21 generated sub-alignments: 2-34 sequences x 25%/50%/100% truncation
```

This is the composition the deployment options below contemplate, and it works
today with no porting at all: install the C `mafft` as now, point
`MAFFT_BINARIES` at a directory of Rust components. The driver's own contract
(`$MAFFT_BINARIES`, falling back to its `libexec`) is what makes that possible,
and nothing in the port had to change to support it.

## What it buys

Best of three, same driver, same inputs:

```
input              C(ms)  Rust(ms)  ratio
sub_34_100.fa       1452       798   0.55x
sub_21_100.fa        640       395   0.62x
sub_13_100.fa        334       234   0.70x
sub_08_100.fa        281       228   0.81x
test/sample         1522       813   0.53x
TOTAL               4229      2468   0.58x
```

1.7x on the shape the pipeline actually feeds, rising with sequence count.
That is wall clock through the real driver, so it includes the process
start-ups and the shell, not just the aligner.

## The driver does not block packaging

A second trace, run independently against a conda-forge 7.526 libexec with
synthetic protein input at three shapes, reached the same nine binaries -- so
two methods agree that `--auto --preservecase` needs nothing outside the crate.

That trace concluded the missing `mafft` driver is therefore "the whole gap",
blocking a conda package because any package would still need conda's `mafft`
script. It does not, and the reason is in this repository's own Makefile:

```
mafft: mafft.tmpl mltaln.h
	sed "s:_LIBDIR:$(LIBDIR):" mafft.tmpl  > mafft
```

The driver is `mafft.tmpl` plus one substitution, it is BSD-3-Clause, and it is
already here. A `provolut-mafft` package runs that same `sed` and ships the
result next to the Rust components -- which is what conda-forge's own recipe
does. No dependency on another package's script, and nothing ported.

It needs no interpreter either: on the whole file the only external commands
are `cat awk grep tr head sed wc mktemp basename uname tail sort dirname`, and
both `perl` references (`mafft:2366`, `:2370`, the `seekquencer` branch) are
commented out.

What the absent driver does block is a **standalone** Rust MAFFT -- `cargo
install mafft-oxide` giving a working `mafft` with no shell script. That is a
real gap and a reasonable product goal. It is not an adoption blocker, and it
is worth weighing against the cost: `--auto`'s strategy selection is part of
the algorithm (the branch table above shows it switching aligners as input
grows), so a driver that chooses differently changes output even with
bit-exact components. Composing with the upstream script gets that selection
right by construction; reimplementing it is a risk taken for packaging
convenience.

## What this does not establish

- One machine, one invocation shape. The consumer's own parity test (135 real
  sub-alignments from ProteinGym) is the gate that matters, and it should be
  run against a composed directory before anything is adopted.
- Only `--auto --preservecase` was traced. Another flag can reach a component
  that is not ported; the branch table above is the map, not a guarantee.
- No `mafft` driver is ported, and on this evidence none needs to be. Porting
  the driver would mean reproducing its strategy selection exactly -- a
  different choice produces different output even with bit-exact components --
  for no benefit over composing with the C script.
