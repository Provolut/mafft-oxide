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
