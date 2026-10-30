# How far the port is from needing no C at all

Written 2026-10-28, in answer to a direct request: make the Rust at feature
parity with the C so no C code is needed, and have a driver as fast as the C
driver.

The second half is already true and is measured below. The first half is not,
it is not a day's work, and this file is the map rather than the claim. Nothing
here is argued from reading the code alone; every row was checked against
`core/` on this machine.

## The driver: done, and it is the C's driver

`bin/mafft-oxide`, as the conda package installs it, is `core/mafft.tmpl` with
one line substituted -- the same single substitution `core/Makefile:168-169`
makes. 3069 lines in, 3069 lines out, and a diff against the C's own generated
`mafft` is empty once that line is normalised.

So "as fast as the C driver" holds by construction, and measuring it only
confirms there is no overhead hiding in the substitution. Same components, best
of 7:

| input | C-generated driver | packaged `mafft-oxide` |
|---|---|---|
| 2 sequences x 4 residues | 29 ms | 30 ms |
| `test/sample`, 36 x 509 | 78 ms | 78 ms |

Identical within noise, which is the expected result for the same script.

What that does *not* mean is that the driver is fast. It is a shell script, and
its fixed cost is the ~29 ms visible in the first row, where the alignment
itself is negligible. An independent measurement on a consumer's production
shapes put it at 39-46 ms, flat with input size, which at ~100 calls per
protein is several seconds of shell per protein. That cost belongs to upstream's
script and is inherited, not introduced.

Removing it means a native driver, and that is the one item here with a
documented argument against it rather than just a cost.
`provolut/DEPLOYMENT.md` makes it: `--auto`'s strategy selection is part of the
algorithm, not a wrapper detail. The table at `core/mafft.tmpl:1294-1325`
switches aligners as the input grows, so a driver that chooses differently
produces a different alignment from bit-exact components. Composing with the
upstream script gets that selection right by construction; reimplementing it is
a correctness risk taken for a packaging convenience. It should be done
deliberately, with a differential over the whole branch table, or not at all.

## The components: 9 of 32, and 13 of the gap is reachable

`core/` builds 32 binaries (plus two Perl scripts). Nine are ported. Of the 23
that are not, 13 are reachable from the driver and 10 are standalone tools it
never invokes -- so "no C needed" for `mafft` itself is 13 binaries, not 23.
Counts are `grep -c "prefix/<name>" core/mafft.tmpl`, which is what the driver
can actually execute.

| unported | driver refs | reached by | what it is |
|---|---|---|---|
| `pairlocalalign` | 14 | `--lara`, `--slara`, `--mxscarna`, `--foldalign`, RNA pairwise | the standalone pairwise aligner. **Not reached by `--auto`:** on the local branch the pairwise stage runs inside `tbfast -L`, which is ported. |
| `addsingle` | 8 | `--add` with some option combinations | incremental addition. The `--add` cases in `verify_pkg.sh` pass with nine binaries, so the combinations it covers do not reach this. |
| `nodepair` | 5 | `--merge`, subalignment paths | pairwise within a guide-tree node |
| `splittbfast` | 2 | `--parttree`, and `--auto` past 20000 sequences | the divide-and-conquer progressive aligner. **The practical edge.** |
| `filter` | 2 | `--seed`, `--mafft-homologs` | seed/homolog filtering |
| `setdirection` | 2 | `--adjustdirection` | reverse-complement adjustment |
| `makedirectionlist` | 2 | `--adjustdirection` | the direction table it consumes |
| `mccaskillwrap` | 2 | `--mccaskill` | RNA secondary structure |
| `setcore` | 1 | `--core` | core-region selection |
| `multi2hat3s` | 1 | `--merge` | constraint table conversion |
| `dndfast7` | 1 | legacy `--fastaguide` paths | the old guide-tree builder |
| `dndblast` | 1 | legacy BLAST-distance path | ditto, from BLAST output |
| `contrafoldwrap` | 1 | `--contrafold` | RNA structure, CONTRAfold |

Never invoked through `$prefix` by the driver, so irrelevant to `mafft`
parity: `sextet5`, `mafft-distance`, `mafft-profile`, `seq2regtable`,
`regtable2seq`, `getlag`, `maffttext2hex`, `hex2maffttext`, `pairash`, `score`.
They are separate entry points or used by the `.rb`/`.pl` wrappers.

### Two refusals inside binaries that are ported

These are not missing binaries; they are paths inside shipped ones that decline
rather than approximate.

- **`dvtditr-rs` refuses `-s` (the specificity scoring path).** Reached by
  `--allowshift` with iterations. Needs `scoringmatrices`, `classifypairs` and
  `partA__align_variousdist`. This is the substantial one -- a scoring family,
  not a formatter.
- **`f2cl-rs` refuses `-d` (`--dash` record exclusion).** It looks like the
  cheapest item on this page and it is not. `commongappick` is already in the
  lib (`src/addfunctions.rs`), so the filtering itself is a few dozen lines.
  The problem is the order in `core/f2cl.c:311-362`: `fillspace` pads every
  record to `nlenmax`, then the DASH filter builds a reordered subset, then
  `commongappick` *shortens* those records in place -- and `setmark_clustal`
  and `clustalout_pointer` are then both called with the original `nlenmax`.
  So the marking loop reads `seq[j][i]` past each shortened record's NUL, into
  whatever the buffer held before `commongappick` moved the terminator. The
  conservation row therefore depends on stale heap contents, and reproducing
  it means keeping that tail rather than truncating. That is the same class as
  items 1 and 2 of `rust/EXTENDING.md`'s notes, and it needs a differential
  built for it rather than an afternoon.

The clustal and phylip writers were the third refusal until today's commit and
are now compared byte-for-byte against the C (860 cases). That is the shape the
remaining two should follow: port, then a differential over the corpus, then
the refusal goes away. It is also the evidence for how long each takes --
clustal and phylip were ~165 lines of C and still turned up one load-bearing
quirk (`nweaker = 2`, README item 14) that only a differential found.

## What is actually true today, stated precisely

- `mafft --auto` is covered for inputs up to 20000 sequences. Verified
  byte-identical through the real driver with all nine components swapped,
  against both this repo's `core/` build and a hash-verified conda-forge
  7.526 libexec.
- Outside that envelope the package exits non-zero naming the missing
  component. Verified, not assumed: a nine-only libexec under `--parttree`
  dies with `splittbfast: No such file or directory`, exit 1, no output. It
  does not produce a different alignment.
- For the one consumer this was packaged for, a grep of its `src/`, `test/`
  and `scripts/` finds zero uses of `--clustalout`, `--phylipout`,
  `--allowshift` or `--parttree`; the only MAFFT flags it passes anywhere are
  `--auto` and `--preservecase`. (Measured by that consumer's own session.)
- **The flags are provably safe for it; the inputs are safe only by current
  data.** One of its entry points carries `max_sequences=None`, and its cached
  alignments run 3831-9655 sequences. Nothing caps the uncapped path, so a
  protein with more than 20000 homologues would select the parttree branch and
  need `splittbfast`. That failure is loud, but it is data-dependent rather
  than configuration-dependent, which is why it is written here.

## Evidence from outside this repository

Recorded here rather than in `rust/README.md`'s harness table, because that
table is for runs reproducible from this tree and these payloads are not in it.
Measured by the consumer's own session, on its machine, with `verify_pkg.sh`
from this commit and a conda-forge 7.526 libexec whose nine components were
checked against the upstream package manifest. Not run here, so attributed
rather than claimed.

    100 real production payloads, 51 x 286, captured from the consumer's
    pipeline by wrapping its aligner variable during one real run:

        pkg: pass=3100 fail=0 refused=100

31 passes and one refusal per payload across 32 cases, and the single refusal
on every one of them is `--allowshift --globalpair --maxiterate 100`. Zero
`--clustalout` refusals, which is the check that the writers added in this
commit hold on real input and not only on the corpus.

The production invocation at true size -- `--auto --preservecase` on a
3831 x 334 input, where `disttbfast` is 98% of the call -- is byte-identical to
the C at 4,542,015 bytes.

Re-measured there on 2026-10-29 after the skip-table fix, same payload, best of
2, fresh process, output byte-identical at those same 4,542,015 bytes:

    C                  10844 ms   1.000
    port, 0.3.0 _2     14682 ms   1.354   the packaged artifact
    port, cde06cc       9054 ms   0.835

So the port leads the C by 1.20x at production size, against 1.14x at 790 here
-- the fix targets the O(n^2) term, so its share grows with n. Attributed, not
claimed: measured on the consumer's machine and payload. The binaries were
staged from this tree's `rust/target/release` under their C names and run
through the shipped driver with `MAFFT_BINARIES`, with an empty-variable canary
failing first so the variable is demonstrably honoured, and `sha256sum -c`
afterwards showing 9 of 9 unchanged during the run. Those binaries are the ones
this repository's own differentials ran against: same sha256, and `cargo build
--release` reports nothing to do at `cde06cc` with a clean tree, so they
correspond to that commit.

**Build 3 carries the fix; build 2 does not.** `0.3.0-h1250391_2` predates it,
so a consumer on that build gets the slow path. `0.3.0-h1250391_3` is built
from `debddb7`, the commit this work merged as, and was installed into a fresh
pixi environment and timed there -- the installed artifact, not this tree's
`target/release` -- on the 790 x ~384 payload, best of 5 interleaved, each arm
byte-identical to the C at 621,040 bytes:

    C, stock                 871 ms   1.000
    packaged _3              756 ms   0.868
    packaged _3 --thread 8   559 ms   0.642
    packaged _2             1175 ms   1.349

The `_2` row is the reason the table exists rather than an aside. Output never
changed, so byte-identity passes for both builds and cannot tell them apart; a
rebuild that had quietly picked up the old source would read 1.349 here and
pass every correctness check in this repository. The clock is the only
discriminator, and it was run as one.

Acceptance, against the installed libexec rather than against this tree:
`verify_pkg.sh --shapes` reports `pass=37 fail=0 refused=1`, the refusal being
the declared `--allowshift --globalpair --maxiterate 100`, with 9 of 9
components differing between the two arms and the empty-`MAFFT_BINARIES`
canary failing first so the variable is demonstrably honoured. Run twice: once
with the stock driver over both libexecs, once with the packaged driver.

Build 3 went out as a build number because the consumer's manifest pinned
`mafft-oxide = "0.3.0.*"`, which accepts one with no `pyproject.toml` edit and
no PR, where a version move would have needed both. The cost of that was a
version string covering two trees: crates.io's 0.3.0 and builds 0 to 2 are the
`v0.3.0` tag, build 3 was cut from `debddb7`.

**0.4.0 closes that, and is what to install.** Same source as build 3, tagged
`v0.4.0`, so the tag, the crate and the package describe one tree again, and
`--thread` is named in the version rather than hidden in a build number --
which it is, since the flag went from exiting 1 to being honoured. Built from
the tag, installed from a fresh pixi environment and measured there, same
payload and method, best of 5 interleaved:

    C, stock                    873 ms   1.000
    packaged 0.4.0              763 ms   0.874
    packaged 0.4.0 --thread 8   562 ms   0.644
    packaged 0.3.0 _3           761 ms   0.872

0.4.0 and 0.3.0 build 3 agreeing to 2 ms is the expected result and worth
recording as one: they are the same source, so a difference would have meant
the tag and the commit disagree. `verify_pkg.sh --shapes` reports
`pass=37 fail=0 refused=1` against the installed 0.4.0 libexec, again run with
both drivers; and `disttbfast` accepts `-C 8-8` returning its own `-C 0-0`
bytes, which is the structural half of the check that needs no clock.

**What the consumer took, and the verb that moves a build number.** They stayed
on `0.3.0.*` and took build 3: no manifest edit and no PR, where `0.4.0.*`
needs both -- and in that repository a dependency change may not ride a
feature branch, so it would have waited for a `chore(deps)` PR. 0.4.0 is still
the better thing to install once such a PR is going out anyway. Same source,
same speed, so the choice was paperwork and not performance.

The command is `pixi update mafft-oxide`, not `pixi lock`. Measured on a
replica rather than reasoned about: with the pin still satisfied by `_2`,
`pixi lock` prints "Lock-file was already up-to-date" and changes nothing,
because a lock that satisfies the manifest is already correct. Only an upgrade
verb moves within a satisfied constraint. A *changed* constraint is different
-- with the pin edited to `0.4.0.*`, plain `pixi lock` moved `_2` to 0.4.0.
Both sessions advising this had the verb wrong before one of them tested it.

**At production size, measured by pml0 on the installed package:** 3831 x 334,
1,202,735 bytes of real pipeline payload, `--auto --preservecase`, best of 2
interleaved per arm, through the environment's own `mafft` and `mafft-oxide`
on PATH, wall clock under `LC_ALL=C`. No cache on either side: the aligner is
a binary invocation and the pipeline's realignment and trace caches sit above
it, untouched by a direct run.

    C          10767 ms   1.000
    oxide _3    9067 ms   0.842

Byte-identical at 4,542,015 bytes. Attributed to pml0, as the 9054 was.

9067 against that 9054, measured earlier from a staged `target/release` and
0.14% apart on the same payload, closes independently the question the `_2`
row answers at 790: the package carries the source it claims to. A rebuild
that had quietly picked up the old tree would have read about 14,700 here, as
`_2` did.

**A lock is not an installation.** After the update the lock named `_3` in all
nine of its references, while in pml0's checkout only `full-cuda` held `_3` on
disk -- `ci-cuda` and `ci-cuda-adept` still held `_2`'s binaries, and the
structural check says so without a clock: `full-cuda`'s `disttbfast` accepts
`-C 8-8`, the other two exit 1 with the old `threads change the result`
refusal.

Which environments exist is a property of a checkout rather than of the
project, which is worth stating because both this file and pml2-b4 have given
a single number for it and no single number is right. Measured across the five
clones:

    checkout  materialised environments                                  mafft-oxide
    pml0      ci-adept ci-cuda ci-cuda-adept default full full-cuda      _3 in full-cuda, _2 in ci-cuda and ci-cuda-adept
    pml1      ci ci-cuda-adept default full full-cuda                    _2 in full-cuda
    pml2      ci ci-cuda-adept default full full-cuda full-cuda-probe    _2 in full-cuda
    pml3      full full-cuda                                             _2 in full-cuda
    pml4      full-cuda                                                  _2 in full-cuda

So `ci-cuda` is real in pml0, where the observation was made, and absent in
pml2, where pml2-b4 looked for it -- neither reading was wrong about its own
tree. Seven materialised environments hold the package in all, six of them at
`_2` and one at `_3`. Every one of those builds produces byte-identical
alignments, so a stale environment is slow and never wrong, and no computed
result anywhere in that table needs revisiting. `pixi run` re-syncs an environment before
running, so this closes itself on next use; what it does not survive is a
binary invoked by absolute path, or a shell activated before the update. Check
what executed rather than what the lock records -- the same lesson as the
empty-`MAFFT_BINARIES` canary.

0.3.0 builds 2 and 3 and 0.4.0 all stay on the channel, which is append-only
once anything has locked against it.

The strategy matrix was *not* run at 3831 sequences, deliberately and after
discussion. `--globalpair --maxiterate 100` at that size is on the order of
seven million pairwise alignments, and it would characterise a strategy nobody
runs at that size rather than the package. A 600 x 297 input takes the
identical `--auto` branch -- `core/mafft.tmpl` sends any `nseq` between 500 and
20000 to the `iterate=0 cycle=2` arm, so 600 and 3831 reach the same components
by the same route -- and the full matrix finishes there. It was run, same
harness, same two arms:

    one 600 x 297 input of the shape the consumer's build step produces:

        pkg: pass=31 fail=0 refused=1

Thirty-one passes and one refusal across the same 32 cases, and the refusal is
again `--allowshift --globalpair --maxiterate 100`. That strategy is now the
only one refused on two unrelated input shapes, which makes the refusal a
property of the port -- `dvtditr -s` is unported, see the table above -- rather
than of a payload. One input is evidence about one shape and one strategy set,
nothing wider.

Two bounded claims rather than one unbounded run: the production path verified
at true size, the strategy matrix verified at a size where it can complete.

## The driver reads the operator's locale, and the answer changes

Not the port, and not the package: upstream's shell driver. Recorded here
because it is the one thing measured in this tree that changes an alignment
without any code changing, and because it invalidates comparisons made across
two machines without anyone doing anything wrong.

Same driver, same input, same binaries, nothing varied but `LC_NUMERIC`, on a
123-sequence corpus through `--auto --preservecase`:

    LC_NUMERIC=C            101549 bytes   742 columns
    LC_NUMERIC=de_DE.UTF-8   95768 bytes   696 columns

Stock C from `core/` gives the same two numbers. Within each locale the C and
the port stay byte-identical to each other, so the port is faithful in both and
the one rule is not in question. It is the driver that changes its mind.

The mechanism is one argument. `core/mafft.tmpl:1507` computes the similarity
offset with `awk`, and mawk's string-to-number conversion honours `LC_NUMERIC`
(gawk does not, outside POSIX mode -- so this depends on which `awk` the host
has). Captured with logging shims, exactly one argument differs across the
whole run:

    C locale : tbfast ... -Q 100.0 -h 0.1 -L ...
    de_DE    : tbfast ... -Q 100.0 -h 0,1 -L ...

The component parses that with C `strtod` -- there is no `setlocale` anywhere in
`core/*.c` -- so `0,1` reads as `0` and the similarity offset is silently
dropped. Measured rather than argued: injecting `-h 0,1` under `LC_NUMERIC=C`
and changing nothing else reproduces the comma-locale output to the byte. The
C-locale answer is the intended one, so a comma-decimal shell is the aberrant
arm -- and it is the arm every driver-level run on this machine has taken.

What it does and does not touch:

- **The golden fixtures and `tests/matrix.rs` are unaffected.** They drive
  components with explicit argv, never through the driver, and 134 tracked
  cases already carry `-h 0.1`. The component-side parse is covered at the
  correct value.
- **Every differential harness result stands.** Both arms of every comparison
  ran in one locale, so what they measured is what they claim to measure.
  `verify_pkg.sh --shapes` was re-run under `LC_NUMERIC=C`, which is the
  parameterisation the comma locale had been hiding from every acceptance run
  made here: `pass=37 fail=0 refused=1` against the installed 0.4.0 package,
  and the same `pass=37 fail=0 refused=1` against installed `0.3.0 _2` -- the
  build provolut-ml's lock actually holds, run because acceptance for a
  version the consumer does not install is evidence about the wrong artifact
  (pml2-b4's flag). Same refusal both times,
  `--allowshift --globalpair --maxiterate 100`. The package is clean at both
  parameterisations and on both builds.
- **Absolute timings are locale-bound.** A figure measured in a comma-decimal
  shell was measured at `-h 0`, on a different alignment from the one a
  C-locale host computes. Ratios between two arms in the same shell are
  unaffected; a number carried from one machine to another is not.

There are in fact **two** locale faults in the driver, not one, and they fail
at different points. Found by pml0-51 on the second and separated by
measurement:

    mawk, LC_NUMERIC=de_DE.UTF-8
      print 0.0 + "0.1"      -> 0        coercion: the value never exists
      print 0.1              -> 0,1      formatting: the value survives, the text does not
    gawk 5.4.1, either locale -> 0.1 and 0.1

**Coercion** is the worse one and the narrower. A quoted string coerced to a
number yields `0`, so the value is gone inside `awk`, before anything reaches
`strtod`; no component-side change could recover it. It is what
`mafft.tmpl:904` and `:914` do, and those lines are gated on `--exp` and
`--ep`, so only a run passing those flags reaches it. The `expr` guard
underneath them does not catch it, because `0` is a number.

**Formatting** is what the `--auto` path hits. A numeric literal in the awk
source keeps its value -- mawk's lexer is locale-independent -- and `print`
then emits it with the locale's decimal separator. That is the `-h 0,1` above,
and the injection test shows it accounts for the whole difference on that path.

Which of the two fires on the default path follows from the two argument slots
in `mafft.tmpl:1507-1509`, and pml2-b4 worked it out:

    aof=`awk "BEGIN{print 0.0 + \"$similarityoffset\" + $aof}"`

`$similarityoffset` lands in the quoted **string** slot, so it goes through the
coercion -- but its default is `"0.0"` (:293), which is 0 either way, so the
fault is invisible there. `$aof` is substituted bare, as an awk numeric
**literal**: mawk's lexer reads it locale-independently, so the value survives
and it is `print` that then renders it with a comma. The nonzero defaults are
all in that literal slot -- `defaultlaof="0.100"` (:91), `pgaof="0.10"` (:206);
`defaultaof="0.000"` (:90) reads as `0` in either locale. So comma-on-output is
the form the default path takes, which is why the captured argv shows `-h 0,1`
and not `-h 0`.

**The driver's own `LANG=C` does not protect it.** `mafft.tmpl:5` is
`LANG=C; export LANG`, which looks like exactly this defence and is not one:
`LANG` is only the fallback an unset category falls back to, and an explicit
`LC_NUMERIC` overrides it. Measured on mawk:

    LANG=C, LC_NUMERIC unset        -> 0.1
    LANG=C, LC_NUMERIC=de_DE.UTF-8  -> 0,1
    LANG=de_DE.UTF-8, unset         -> 0,1

This also bounds the `env -i` test above. That test set no locale at all, so
the driver's `LANG=C` won by default and the run took the correct arm -- which
makes it sound evidence for the claim it was making, that the launcher
self-locates without an activation environment, and no evidence at all about
the locale. The realistic unit is worse than it, because **systemd propagates
`/etc/default/locale` into services** and on this machine that file sets
`LC_NUMERIC=de_DE.UTF-8`. Measured with exactly that environment and systemd's
default unit PATH: 95768 bytes, the broken arm. Adding one line,
`Environment=LC_NUMERIC=C`, restores 101549. So the fix is a unit line and it
is verified to work, rather than inferred to.

**Which `awk` runs decides whether either fires, and PATH decides that.**
gawk is locale-independent here, mawk is not. pml0-51's production figure --
0.842x, byte-identical at 4,542,015 bytes -- was measured in a shell whose
PATH put a pixi environment ahead of `/usr/bin`, and that environment ships
gawk 5.4.1, so the driver got gawk and the intended `-0.1`. The figure
describes the correct parameterisation after all; the inference that it was
measured at `-h 0` was wrong, and the thing that decided it was PATH order
rather than the locale.

That is also what joins this to the service-host gap, and it is the reason to
read the two as one finding rather than two. A host that execs the launcher by
absolute path from a systemd unit loads no activation environment, so it gets
`/usr/bin/awk`. On this machine that is mawk. Off the activation environment
you therefore lose two things at once: the oxide binaries, and the similarity
offset -- and only the first is visible. The second changes the alignment and
says nothing.

**This package is exposed to it.** `recipe/meta.yaml` declares no `run:`
section, on the stated grounds that the driver's external commands are
"coreutils/busybox-level". That is true as a claim about *presence* and false
as a claim about *the answer*: the installed environment ships no `awk` at
all, so the driver resolves whichever one the host has, and which one that is
changes the output.

The first reason given here for leaving it unpinned was wrong, and pml2-b4
inverted it from the packaging metadata. It read: a `gawk` run-dependency
would make this package answer differently from a stock `mafft` on a mawk
host. But the comparison that matters is against the *package* this one
substitutes for, and conda-forge's `mafft` pins `gawk`:

    mafft-7.526-h4bc722e_0        depends: __glibc, gawk, libgcc-ng
    mafft-oxide-0.4.0-h1250391_0  depends: __glibc, libgcc

Read from installed metadata, not from a recipe. conda-forge's `gawk 5.4.1`
ships `bin/awk` as well as `bin/gawk`, and that gawk is locale-independent for
both faults -- verified, `0.1` under `de_DE.UTF-8` for coercion and
formatting alike. So installing stock `mafft` from conda-forge *brings* a
locale-proof `awk` into the environment, and the mafft feedstock evidently
knows why. Faithfulness to the thing being replaced therefore argues **for**
the run-dependency; omitting it is what makes this package behave differently
from stock. The exposed case is the one this work is heading for: a lean
service host that installs the port and no stock `mafft`.

One limit on what the run-dependency would buy, because it is not a substitute
for the unit line. The driver calls bare `awk`, resolved through `PATH`. A
`gawk` inside the environment is only found when `<env>/bin` is on `PATH` --
measured: with a gawk-providing bin directory ahead of the system one, the
launcher returns the correct arm under `de_DE.UTF-8`. A unit that execs
`<env>/bin/mafft-oxide` by absolute path and leaves `PATH` alone never searches
that directory and still gets `/usr/bin/awk`. So the run-dependency closes the
activated-environment case and `Environment=LC_NUMERIC=C` closes the unit case,
and neither covers the other.

Declared as of 0.4.0 build 1, after the reasoning above was corrected. No
source change came with it: the same tag, the same nine binaries, and all nine
`.text` sections byte-identical to build 0's -- the only raw difference is four
bytes of padded build-prefix string per binary, which conda rewrites at install
time. What moved is which `awk` the driver finds.

How much that was riding on luck, measured in a consumer environment rather
than argued from the feedstock. In `pml0`'s `full-cuda` env the awk on `PATH`
prints `0.1`, so the runs there were on the correct arm -- but:

    bin/awk -> gawk                       GNU Awk 5.4.1   prints 0.1
    owned by                              gawk-5.4.1-h0a3468a_0
    mafft-7.526-h4bc722e_0       depends: __glibc, gawk, libgcc-ng
    mafft-oxide-0.3.0-h1250391_3 depends: __glibc, libgcc

The locale-proof `awk` keeping those alignments correct is a transitive
dependency of **the stock C `mafft`**. The installed port declares nothing that
holds it there. So the exposed case named above is not hypothetical and not
distant: it is one manifest edit away, and the edit is the obvious one. Dropping
stock `mafft` on the reasoning that the port is at parity takes `gawk` with it,
`/usr/bin/awk` wins the search, and every alignment on the host moves -- 101549
bytes to 95768 on the input measured here, with no error and no warning. That
ordering constraint is the practical content of this section: the build-1 pin
has to land before the stock package is removed, not after.

Nothing in this tree pins the locale -- not a harness, not the recipe, not a
hook -- and nothing here changes that. Pinning `LC_NUMERIC=C` inside the
packaged launcher would make `mafft-oxide` answer differently from stock
`mafft` under a comma locale, which is exactly the divergence this repository
exists to prevent, however much better the pinned answer is. That trade is a
human's to make. What a caller can do without touching either program is pin
`LC_NUMERIC=C` itself, which is the recommendation.

## The parity evidence at 0.4.0, and the pin that did not hold

The differential suite had never run against the threading code. `pre-push`
says so in as many words -- `no recorded store; the full differential did not
run` -- and exits 0 anyway, which is worth stating plainly: a check that
reports it could not run is indistinguishable, at the call site, from one that
passed. It printed that line on every push since the threading work landed.

Run at `v0.4.0`, all four corpus differentials, `LC_NUMERIC=C` pinned and
absolute paths throughout:

    verify.sh         pass=12018  fail=0      88 corpus files
    verify_tbfast.sh  pass=1690   fail=0      40 corpus files
    verify_threads.sh pass=2752   fail=0  skipped=0
    verify_dvtditr.sh pass=32     fail=0      of 32 staged runs
                      -----------------
                      16492 comparisons, no failures

Still pending, and named rather than implied: `verify_e2e.sh` and the four
f64-bit-pattern harnesses (`verify_gaplen.sh`, `verify_flh.sh`, `verify_td.sh`,
`verify_tdi.sh`), which compare at function level and need a `core/` built with
the file-local statics exposed. So this is four of nine harnesses, not "parity
at 0.4.0".

What makes those counts attestable rather than merely plausible, since a clean
run proves nothing about bytes it did not actually execute:

- **The source is the tag.** `git diff v0.4.0 HEAD -- rust/` is empty; the only
  commits since are CI, docs and this recipe.
- **The binaries were built from it, not left over.** All nine post-date the
  tag commit by 7-11 seconds and post-date every source file. A cargo build
  that finds everything fresh does not relink, so a stale artifact would carry
  an mtime *earlier* than its sources; the inequality runs the other way.
- **And they are reproducible.** `cargo clean` (10550 files, 2.1 GiB) then
  `cargo build --locked --release` on the same toolchain and the same absolute
  target directory reproduced all nine sha256 digests exactly. Pinning those
  variables is what makes the comparison mean something: a Rust rebuild is not
  reproducible across toolchain versions, `RUSTFLAGS` or target paths, so
  differing digests would have licensed only "these were not the same build",
  not "the tested binary was stale".
- **The reference arm too.** `core/` is from an older build, which is fine for
  an unmodified upstream tree but only once checked:
  `find core -maxdepth 1 -name '*.c' -newer core/disttbfast` is empty.
- **No case read another's output.** All four harnesses use `mktemp -d` with an
  EXIT trap and `rm -rf` per case, so nothing writes into a populated
  directory.

### The pin that did not hold

The consumer manifest pins `mafft-oxide = "0.3.0.*"`, and its comment explains
that this admits a rebuild of the same version but not a version move, "because
its byte-for-byte agreement with the stock build is established per version".
That reasoning is false, and this recipe is why.

`_2` was built from tag `v0.3.0`. `_3` was built from a bare commit that
already carried the threading work: 566 insertions across 8 files, `score.rs`
rewritten by 231 lines, `parallel.rs` new. Both carry the version string
`0.3.0`. So the 566 insertions were *published* as a rebuild, under a pin whose
whole purpose was to refuse code movement.

How far they actually travelled is narrower than that, and the correction is
worth keeping because the first version of this section overstated it. The
consumer's **committed lock still names `_2`**, so every machine that
materialises `--locked` gets the tag build. Exactly one checkout runs `_3`,
because someone ran `pixi update mafft-oxide` standing in it, and that lock
change is still unstaged. Measured three ways in that checkout: `git show
HEAD:pixi.lock` gives nine `_2` references, the working tree gives nine `_3`,
and the installed `conda-meta` says `_3`. So the exposure was one machine, not
a fleet, and the lock is what held the line -- not the pin, which had already
been crossed at publish time.

That relocates where the remaining risk lives, and it is the part the recipe
rule above cannot reach. The rule stops the *package* being built that way. It
does nothing about a lock that already names such a build, and committing that
lock would have moved the 566 insertions onto every machine in the fleet under
an unchanged `0.3.0.*` manifest line, with no differential ever having run
against them. The shape that would have carried it is a `chore(deps)` PR: nine
build-string lines, no manifest change, nothing on its face for a reviewer to
stop at. A version pin can be defeated through the lock without the manifest
ever being touched.

One traceability gap follows. Nothing in a produced number says which build
made it: the consumer stamps a `code_version` from its own git tag onto run
records, and that string is identical whether the box was on `_2` or `_3`. So
"which numeric core produced this result" cannot be recovered from a run
record after the fact.

An earlier revision of this paragraph said the package should state its own
build identity to close "our half" of that. That was wrong twice over, and the
check that settles it is one file read.

**The identity is already on disk, per environment, and it is better than a
build string.** `<prefix>/conda-meta/mafft-oxide-*.json` carries the build
string, the build number, the exact download URL and the artifact's `sha256` --
so the installed bytes are identifiable by hash, not merely by label, without
running anything and without the package printing anything. Nothing needs to
change here for a caller to record which aligner it used; what is missing is a
stamp in the run record, which is the consumer's surface and not this one.

**And the change it proposed would not have been free.** The driver reads the
`version` component's stdout and compares it against its own baked string --
`versionbin=\`"$prefix/version"\`` then
`expr "$version" : v"$versionbin"` at `core/mafft.tmpl:49-59`, warning
`v$versionbin != $version` on a mismatch. So making `version-rs` print a build
identity on stdout would trip upstream's own consistency check on every run.
The surface that looks inert is load-bearing. If a packaged build identity were
ever wanted anyway, it belongs on a flag stock does not have, or a stream
nothing parses -- not in what `version` prints.

Two things follow. A `<version>.*` pin is only as strong as the publisher's
discipline about what a build number means -- the rule is in `recipe/meta.yaml`
now, and it is that a build number may only repackage the same code. And the
move from `_2` to `0.4.0_1` is a genuine code move rather than a relabel, which
is why the parity run above exists and why the component-level counts are
necessary evidence without being the whole gate: the consumer's own gate is a
pipeline-level comparison of its output through the port against stock `mafft`,
and that one belongs to the consumer, not here.

## If parity is the goal, the order that makes sense

Ranked by what it unblocks per unit of risk, not by size.

1. **`splittbfast`.** The only unported binary `--auto` can reach on its own,
   and the whole of the >20000-sequence envelope. Biggest correctness surface
   of the three big ones.
2. **`dvtditr -s`.** Closes the last refusal that a plausible flag
   (`--allowshift`) reaches, and it is scoring code, so it is where a quiet
   numerical difference is most likely. Deserves the `verify_dvtditr.sh`
   treatment before anything else.
3. **`addsingle`** and **`nodepair`**, together: `--add` and `--merge` are the
   two features most likely to be asked for next, and both already have
   partial coverage through `tbfast`.
4. **`pairlocalalign`.** 14 driver references makes it look urgent and it is
   not: `--auto` never reaches it, and its RNA options need external programs
   (`lara`, `mxscarna`, `foldalign`, CONTRAfold) that are not in this tree at
   all. Only worth doing if someone asks for those strategies.
5. **`setdirection` + `makedirectionlist`.** `--adjustdirection` is a pair of
   small binaries and a self-contained feature.
6. The RNA structure wrappers and the legacy guide-tree builders last. Each
   depends on a third-party program or on a path upstream has superseded.

Every one of these should land the way the clustal writer did: port it, run a
differential over the generated corpus against `core/`, and only then let the
refusal go. A half-ported binary that answers is worse than one that refuses,
and this repository exists because of that distinction.
