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

**The packaged artifact still carries the 1.354x.** `0.3.0-h1250391_2` predates
this fix, so a consumer installing the conda package gets the slow build path
until a new one is published. The consumer's manifest pins
`mafft-oxide = "0.3.0.*"`, which accepts a new *build number* with no manifest
edit -- a `_3` needs only a re-lock -- but refuses a version move, so 0.3.1 or
0.4.0 would need a change to their `pyproject.toml` and therefore a PR.
Publishing this as a 0.3.0 rebuild is the cheaper path for them. Both the
rebuild and the re-lock are the user's to run.

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
