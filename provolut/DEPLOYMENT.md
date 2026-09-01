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
toolchain. Costs a build pipeline and somewhere to host the artifact.

**5. Upstream the patches.** The right long-term home. They are bit-exact
against upstream's own test suite, so they are contributable, and if accepted
everyone's MAFFT gets faster and the fork disappears. Slowest path, and not
mutually exclusive with any of the above.

## What has to be decided, by a person

- **Is cold-path scoring latency actually a problem?** If the answer is no,
  option 1 is correct and the rest of this document is premature.
- **If it is: (4) or (5)?** (3) is the tempting middle and the worst of both --
  every machine builds its own copy, and drift between them is invisible.
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
