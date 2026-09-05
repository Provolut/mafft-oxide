# What still needs the C, and what porting it would take

`mafft --auto` runs with no C process. Six configurations outside that path do
not, and each one fails loudly, naming the setting — the port never runs them
and quietly returns something else.

This is the successor to `GOAL-C-FREE.md` §5's refusal list. It is ordered by
what a user is most likely to hit, and every size is measured from the source
rather than guessed. The rules in `README.md` — transliterate, never
reassociate f64, verify with a differential before believing anything — apply
to all of it.

---

## 1. `--thread n` where `n > 0` — **not a porting task**

The refusal is `tbfast-rs`'s and `disttbfast-rs`'s (`-C > 0`) and, for
`dvtditr`, the absent `athread`.

**Do not start by porting `athread`.** Measure first, because the C's own
answer moves:

```
$ for i in 1 2 3; do mafft --auto --thread 4 test/sample > out.$i; done
run1 vs run2: DIFFER
run1 vs run3: DIFFER
```

Three runs of one input, three different alignments. The same test with
`--maxiterate 0` — which never reaches `dvtditr` — is identical across runs, so
the nondeterminism is `dvtditr`'s alone, and `--thread 4` also differs from
`--thread 0`.

That makes bit-exactness against the C **impossible in principle** here, not
merely unimplemented. The repo's one rule is output verified byte-identical,
never argued equivalent; there is nothing to verify against when the reference
has no single answer. `BAATARI2` is why: worker threads race to claim branches,
and which worker gets which branch at which generation decides the accept
order.

So the options are, in order of honesty:

1. **Leave it refused** (today's behaviour). `-C 0` is what `scripts/mafft`
   passes by default and what every verified number in this repo was produced
   at.
2. **Implement a threaded refinement of your own** and drop the byte-identity
   claim *for that mode only*, with the ladder redefined for it — e.g. the same
   SP score within a tolerance, not the same bytes. That is a product decision,
   not a transliteration, and it needs saying out loud in `README.md`.
3. **Accept `-C n` and run sequentially.** Wrong for `dvtditr` — the C at `-C n`
   genuinely produces a different alignment, so matching it means matching a
   moving target — but *defensible for `tbfast` / `disttbfast`*, whose threaded
   output is byte-identical to their sequential output (measured: 39 corpus
   inputs x 3 repeats, `pre`, `order`, `hat2` and `hat3` all identical). A user
   passing `--thread 8 --maxiterate 0` would then get the right answer at
   one-core speed instead of an error.

If you want the wall-clock, option 3 for the progressive stage plus a plain
refusal in `dvtditr` gets most of it, and keeps every existing guarantee.

## 2. `--allowshift` — the specificity scoring path

The strategy passes `dvtditr -s 0.8`. `specificityconsideration != 0` switches
`tditeration` onto a second scoring path throughout, so this is the largest of
the remaining items and the only one that touches the DP.

| piece | C source | size |
|---|---|---|
| `partA__align_variousdist` | `partSalignmm.c` | 746 lines |
| `makescoringmatrices`, `calcmaxdistclass` | `tditeration.c`, `mltaln9.c` | ~90 |
| `classifypairs` | `tditeration.c` | 65 |
| `intergroup_score_multimtx` | `mltaln9.c` | ~60 |
| `smalldistmtx` half of `OneClusterAndTheOther_fast` (the `USEDISTONTREE` branch) | `tddis.c` | ~15 |

Where to be careful:

- `partA__align_variousdist` is **not** `partA__align` with a matrix argument.
  Diff the two before assuming otherwise; `partalign.rs`'s existing note about
  `Salignmm` vs `partSalignmm` tracebacks is the warning shape to expect.
- `OneClusterAndTheOther_fast` already computes both member lists; only the
  `smalldistmtx` fill is missing, and it is `distontree[im] + distontree[jm]`
  under `#define USEDISTONTREE 1`, *not* the `distmtx` lookup the `#else`
  branch shows.
- Verification is cheap because the harnesses exist: extend
  `gen_tdiharness.c` / `gen_tdicases.py` with a non-zero `-s` and the same
  126-case corpus runs against it.

## 3. `--add` (`-K n`)

`dvtditr-rs` asserts `nadd == 0`. What is missing is the `skipthisbranch`
computation in `dvtditr.c`'s `main` (~40 lines): build `reftable` from the
first `njob - nadd` sequences, then mark every branch whose member list
overlaps the reference set without matching it, using `overlapmember` and
`samemember` from `mltaln9.c` (~20 lines together).

Small and self-contained. `verify_dvtditr.sh` needs one added input — an
alignment plus a few extra sequences, run through `mafft --add` — and the
staged replay then covers it like everything else.

## 4. Inputs of 200000 sequences or more — `splittbfast`

`--auto` picks the parttree bracket there and runs `splittbfast`, a binary this
port has never covered and which `GOAL-C-FREE.md` scoped out from the start. It
is a separate program with its own tree machinery, not a mode of the ported
ones.

Nothing about the rest of the port blocks it; it is simply a second project of
comparable size. Note the practical ceiling before starting: the corpora here
top out in the thousands, so the differential would need inputs nobody in this
repo currently generates.

## 5. Unconstrained input over 50000 residues — the memsave switch

`dvtditr.c` silently sets `alg = 'M'` above 50000 residues when there is no
`hat3`, which selects the memsave aligner — a *different alignment*, not a
cheaper route to the same one (`README.md`, "Notes for anyone extending this",
item 12). `dvtditr-rs` refuses rather than keep running `'A'`.

**This is probably the cheapest of the five.** `falign` already implements the
memsave path (`alg == b'M'`: `NKOUHO_LONG`, no `break` on an empty lag, every
block through `MSalignmm`), and `msalign.rs` is ported and verified through
`tbfast`. The constrained path never reaches it — the C's switch is guarded by
`!constraint` — so `partA__align` needs nothing.

Likely the whole change: relax `tditeration`'s `alg == b'A'` assertion to
accept `b'M'`, thread it to the existing `falign` call, and mirror `dvtditr`'s
`nlenmax > 50000 && !constraint` switch in `dvtditr-rs`. Then verify: the
existing corpora are far too short, so `gen_tdicases.py` needs a case above the
threshold, and `verify_dvtditr.sh` an input to match. Do not ship it on the
strength of it compiling — that is exactly how `partA__align` sat unverified
for a release.

## 6. `--qinsi` / `--xinsi`

RNA secondary-structure strategies that shell out to external folding programs
(`mccaskillwrap`, `contrafoldwrap`). They fail identically with the stock C
build on a machine without those programs, so there is nothing to match and
nothing to port until the external tools are in the picture.

---

## How to check this list is still true

It was derived by tracing, not by reading, and the same commands re-derive it:

```
# which binaries a strategy actually executes
cp -r binaries logbin && for f in logbin/*; do ... exec wrapper ... done
EXEC_LOG=/tmp/exec.log MAFFT_BINARIES=$PWD/logbin scripts/mafft --auto test/sample

# what dvtditr was really called with
wrap binaries/dvtditr in a script that records "$@" and stdin, then execs the real one

# which C functions a run enters
valgrind --tool=callgrind --callgrind-out-file=cg.out binaries/dvtditr <argv> < stdin
grep -aE '^(fn|cfn)=' cg.out | sed -E 's/^c?fn=//; s/^\([0-9]+\)[[:space:]]*//' | sort -u
```

A configuration that stops appearing in those outputs has left the reachable
surface; one that appears and is not in this file is a gap in the file.
