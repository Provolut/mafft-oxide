#!/bin/bash
# Packaging check: run the real `mafft` wrapper twice -- once against a libexec
# directory of C binaries, once against a packaged libexec in which *all nine*
# ported components carry their C names -- and require byte-identical output.
#
#   verify_pkg.sh <c-libexec-dir> <packaged-libexec-dir> [--shapes] [inputs...]
#
# This is the all-nine sibling of `verify_e2e.sh`, which swaps one component
# (two or three with TBFAST_RS / DVTDITR_RS set) into a copy of the C libexec.
# A conda package swaps every one of them at once and ships its own driver, so
# the thing that needs checking is a whole installed directory rather than a
# single binary. The strategy list, the `--add` cases and the corpus are
# verify_e2e.sh's, deliberately: the question here is whether the *packaging*
# changed an answer, so the inputs must be the ones already known to the
# single-swap harness.
#
# With `--shapes` in place of the input list -- after the two directories, not
# before them -- the harness builds its own corpus
# and adds a second pass that covers the driver's `--auto` strategy table. That
# table is the reason an all-nine check needs more than one input: `--auto`
# chooses from sequence count and length (`core/mafft.tmpl:1294-1325`) and each
# branch runs a different set of components, so a single small input exercises
# one branch and reads as complete. See `gen_shapes` below.
#
# Set MAFFT=<abs path to the driver> to choose which wrapper runs. Both arms run
# the SAME driver -- the comparison is of components, not of scripts. To check a
# packaged driver against the stock one, run this twice with different MAFFT.

set -u
CLIB="${1:?usage: verify_pkg.sh <c-libexec-dir> <packaged-libexec-dir> [inputs...]}"
PLIB="${2:?usage: verify_pkg.sh <c-libexec-dir> <packaged-libexec-dir> [inputs...]}"
shift 2

here=$(cd "$(dirname "$0")" && pwd)
rust=$(cd "$here/.." && pwd)
repo=$(cd "$rust/.." && pwd)

# The wrapper runs each alignment from a scratch directory of its own, so a
# relative `MAFFT_BINARIES` resolves to somewhere that does not exist and every
# reference run fails -- which this harness would report as failures rather than
# as its own misuse. verify_e2e.sh learned that the hard way. Absolute from here
# on.
case "$CLIB" in /*) ;; *) CLIB="$(pwd)/$CLIB" ;; esac
case "$PLIB" in /*) ;; *) PLIB="$(pwd)/$PLIB" ;; esac
[ -d "$CLIB" ] || { echo "verify_pkg.sh: no such libexec dir: $CLIB" 1>&2; exit 2; }
[ -d "$PLIB" ] || { echo "verify_pkg.sh: no such libexec dir: $PLIB" 1>&2; exit 2; }

# The nine components the port provides, under the names the driver calls.
NINE=(disttbfast tbfast dvtditr dndpre countlen f2cl replaceu restoreu version)

# Three guards, and all three exist because a comparison harness that cannot
# notice it is comparing something to itself reports success in exactly the case
# it was written to catch.
#
#  1. Every component must be present and executable in both arms. A missing one
#     makes the driver die with "No such file or directory" (exit 1) on whatever
#     strategy reaches it, which would be counted as a failure rather than named.
#  2. Every component must be a real ELF executable. A shell wrapper that
#     `exec`s the other arm's binary passes a byte comparison trivially and
#     makes both arms the same implementation. That has happened on this
#     machine: nine C binaries in a pixi libexec were replaced in place by
#     two-line `exec` shims, and every measurement taken through them compared
#     mafft-oxide against mafft-oxide.
#  3. The two arms must not be byte-identical. If they are, one directory was
#     passed twice, or the package installed the C binaries, and every case
#     below would pass while proving nothing.
same_count=0
for b in "${NINE[@]}"; do
  for d in "$CLIB" "$PLIB"; do
    [ -x "$d/$b" ] || { echo "verify_pkg.sh: $d/$b missing or not executable" 1>&2; exit 2; }
    case "$(head -c 4 "$d/$b" | od -An -tx1 | tr -d ' \n')" in
      7f454c46) ;;
      *) echo "verify_pkg.sh: $d/$b is not an ELF executable -- a wrapper script here" 1>&2
         echo "               makes both arms the same implementation. Refusing to run." 1>&2
         exit 2 ;;
    esac
  done
  cmp -s "$CLIB/$b" "$PLIB/$b" && same_count=$((same_count+1))
done
if [ "$same_count" -eq "${#NINE[@]}" ]; then
  echo "verify_pkg.sh: all ${#NINE[@]} components are byte-identical between" 1>&2
  echo "               $CLIB and $PLIB -- this would compare an implementation" 1>&2
  echo "               against itself. Refusing to run." 1>&2
  exit 2
fi
echo "verify_pkg.sh: ${#NINE[@]} components present in both arms, $((${#NINE[@]}-same_count)) of them differing"

# The `--auto` shape corpus. Built from this repository's own generator and
# `test/sample` rather than from anything invented here: `gencorpus.py` is what
# verify.sh and verify_tbfast.sh use, and the two large inputs are plain
# concatenations of its `subset_*` files in sorted order, so the whole set is
# reproducible from the repository.
#
# The point is which branch of the driver's table each one lands in:
#
#   input                nseq x nlen   branch taken by --auto
#   test/sample            36 x 509    nlen<3000  & nseq<100   local,   it=1000
#   shapes_123.fa         123 x 509    nlen<1000  & nseq<200   local,   it=2
#   shapes_790.fa         790 x 509    nseq<20000             ktuples, it=0 cyc=2
#
# The first two reach tbfast and dvtditr; the third reaches disttbfast, which is
# the component a production-shaped input spends almost all of its time in. One
# input of the first shape alone would leave the third branch unvisited while
# the tally looked complete.
gen_shapes() {
  local d="$1"
  # gencorpus.py writes straight into the directory it is given and does not
  # create it, exactly as verify.sh's gen_corpus assumes.
  mkdir -p "$d/c"
  python3 "$rust/gen/gencorpus.py" "$repo/test/sample" "$d/c" full >/dev/null
  cat "$d/c"/subset_0[0-5].fa > "$d/shapes_123.fa"
  cat "$d/c"/subset_*.fa      > "$d/shapes_790.fa"
  # An empty shape file turns every case below into "reference run failed",
  # which reads as six divergences rather than as a corpus that was not built.
  for f in "$d/shapes_123.fa" "$d/shapes_790.fa"; do
    n=$(grep -c '^>' "$f" 2>/dev/null || echo 0)
    [ "$n" -ge 100 ] || { echo "verify_pkg.sh: shape corpus $f has $n sequences; gencorpus.py failed" 1>&2; exit 2; }
  done
}

SHAPES=0
if [ "${1:-}" = "--shapes" ]; then SHAPES=1; shift; fi

INPUTS=("$@")
if [ ${#INPUTS[@]} -eq 0 ]; then
  INPUTS=("$repo/test/sample")
fi

# The wrapper itself: $MAFFT if set, else the `mafft` next to the C libexec dir,
# else whatever is on PATH. Refuse to run without one -- a missing wrapper makes
# both runs fail identically, which would otherwise be counted as a pass.
MAFFT="${MAFFT:-}"
if [ -z "$MAFFT" ]; then
  cand="$(cd "$CLIB/../.." 2>/dev/null && pwd)/bin/mafft"
  if [ -x "$cand" ]; then MAFFT="$cand"; else MAFFT="$(command -v mafft || true)"; fi
fi
if [ ! -x "$MAFFT" ]; then
  echo "verify_pkg.sh: cannot find the mafft wrapper; set MAFFT=<path>" 1>&2
  exit 2
fi
case "$MAFFT" in /*) ;; *) MAFFT="$(pwd)/$MAFFT" ;; esac
echo "verify_pkg.sh: driver $MAFFT"

# Guard four, and the only one that is about the *run* rather than the two
# directories. The three above prove the arms hold different binaries; none of
# them can see a driver that never consults MAFFT_BINARIES at all. A driver
# with its component directory wired in some other way -- a packaged launcher
# whose generated prefix line overwrote the `if [ "$MAFFT_BINARIES" ]` branch,
# say -- would run its own components for both arms, and then every case below
# compares an implementation against itself and passes. That is the same false
# pass the ELF guard exists for, reached by a route the ELF guard cannot cover,
# and output alone cannot distinguish it from a correct port: byte-identity is
# the expected result either way.
#
# Pointing the variable at an empty directory settles it in one run. A driver
# that honours it finds no components and fails; a driver that ignores it
# succeeds. No reliance on the two arms differing in any observable way, so
# this keeps working as components stop refusing things.
canary=$(mktemp -d)
# Always test/sample, never "$@": by this point the positional arguments may
# start with `--shapes`, and a flag where a filename belongs makes the driver
# fail for the wrong reason, which would pass this guard vacuously.
MAFFT_BINARIES="$canary" "$MAFFT" --retree 2 --maxiterate 0 "$repo/test/sample" > "$canary/out" 2>/dev/null
cst=$?
csz=$(wc -c < "$canary/out")
rm -rf "$canary"
if [ $cst -eq 0 ] && [ "$csz" -gt 0 ]; then
  echo "verify_pkg.sh: the driver aligned successfully with MAFFT_BINARIES pointing at an" 1>&2
  echo "               empty directory, so it is not resolving components through that" 1>&2
  echo "               variable. Both arms below would run the same components and every" 1>&2
  echo "               case would pass while comparing an implementation against itself." 1>&2
  echo "               Refusing to run." 1>&2
  exit 2
fi
echo "verify_pkg.sh: driver honours MAFFT_BINARIES (empty-dir canary failed as it must)"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# verify_e2e.sh's strategy list, unchanged. Each entry reaches a different part
# of the port, and between them they reach all nine: the plain and `--retree`
# runs go through disttbfast, `--maxiterate` adds dndpre and dvtditr,
# `--localpair` / `--globalpair` / `--genafpair` run tbfast, `--clustalout`
# runs f2cl, and countlen, replaceu, restoreu and version run on every
# invocation the driver makes.
STRATS=(
  "--nofft --retree 2 --maxiterate 0"
  "--nofft --retree 1 --maxiterate 0"
  "--nofft --retree 2 --maxiterate 2 --thread 0"
  "--nofft --retree 2 --maxiterate 0 --reorder"
  "--nofft --retree 2 --maxiterate 0 --clustalout"
  "--nofft --retree 2 --maxiterate 0 --inputorder"
  ""                                    # plain mafft: FFT-NS-2
  "--retree 1 --maxiterate 0"           # FFT-NS-1
  "--maxiterate 2 --thread 0"           # FFT-NS-i
  "--retree 2 --maxiterate 0 --reorder"
  "--retree 2 --maxiterate 0 --clustalout"
  "--retree 2 --maxiterate 0 --clustalout --reorder"
  "--retree 2 --maxiterate 0 --phylipout"
  "--retree 2 --maxiterate 0 --phylipout --reorder"
  "--localpair --maxiterate 0"
  "--localpair --maxiterate 2 --thread 0"
  "--localpair --maxiterate 100 --thread 0"
  "--globalpair --maxiterate 0"
  "--globalpair --maxiterate 2 --thread 0"
  "--globalpair --maxiterate 100 --thread 0"
  "--allowshift --globalpair --maxiterate 100 --thread 0"
  "--genafpair --maxiterate 0"
  "--genafpair --maxiterate 2 --thread 0"
  "--oldgenafpair --maxiterate 0"
  "--auto"
)

# One strategy is refused by the port rather than approximated, and a packaged
# libexec always contains the binary that refuses it -- so unlike verify_e2e.sh,
# where the refusing component may or may not be swapped in, this is
# unconditional. It is counted separately: not as a pass, which would hide it,
# and not as a failure, which would be untrue.
#
#  * `--allowshift` passes `dvtditr -s 0.8`, which selects the specificity
#    scoring path (`scoringmatrices`, `classifypairs`,
#    `partA__align_variousdist`) that `dvtditr-rs` refuses.
# `--clustalout` was the second, until `f2cl-rs` grew the clustal and phylip
# writers. It is now a compared case rather than a refused one, and
# `--clustalout --reorder` / `--phylipout` / `--phylipout --reorder` are here
# with it -- `--reorder` makes the driver pass `-r order`, which is a different
# path through both writers. It is worth recording why that pair was ever
# refused, because it is the reason this harness exists separately:
# verify_e2e.sh keeps the C `f2cl`, so its `--clustalout` cases passed while a
# package shipping `f2cl-rs` under that name could not serve them at all. An
# all-nine harness saw it on its first run; a one-swap harness structurally
# cannot.
#
# The refusal is only accepted as one on the evidence that it *is* a refusal:
# non-zero exit AND an empty stdout. A port that wrote a plausible-looking but
# different alignment and then exited non-zero would otherwise be tallied here
# instead of failing, which is precisely the outcome this repository exists to
# prevent. verify_e2e.sh checks the exit status alone; this is deliberately
# stricter, and the extra check costs nothing.
REFUSED_BY_PKG=("--allowshift")

is_refusable() {
  local s="$1" pat
  for pat in "${REFUSED_BY_PKG[@]}"; do
    [ "${s#*$pat}" != "$s" ] && return 0
  done
  return 1
}

pass=0; fail=0; refused=0; failed=""; refusals=""
for f in "${INPUTS[@]}"; do
  for i in "${!STRATS[@]}"; do
    s="${STRATS[$i]}"
    MAFFT_BINARIES="$CLIB" "$MAFFT" $s "$f" > "$work/o.c" 2>/dev/null
    cst=$?
    MAFFT_BINARIES="$PLIB" "$MAFFT" $s "$f" > "$work/o.p" 2>/dev/null
    pst=$?
    label="$(basename "$f") [$s]"
    if [ $cst -ne 0 ] || [ ! -s "$work/o.c" ]; then
      # The reference run itself failed, so the comparison would be vacuous.
      fail=$((fail+1)); failed="$failed\n  $label (reference run failed: exit $cst)"
    elif [ $pst -ne 0 ] && is_refusable "$s" && [ ! -s "$work/o.p" ]; then
      # Expected refusal, and only for the strategies named above, and only
      # when nothing was written: any other non-zero exit still fails.
      refused=$((refused+1)); refusals="$refusals\n  $label"
    elif [ $cst -ne $pst ]; then
      fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $pst)"
    elif ! cmp -s "$work/o.c" "$work/o.p"; then
      fail=$((fail+1)); failed="$failed\n  $label (output differs)"
    else
      pass=$((pass+1))
    fi
  done
done

# --add takes two files, so it needs its own loop: align with the driver first,
# then split the result and add the tail back. `--keeplength` also writes a
# `.map` file next to the input, which is part of the output and is compared
# too. `--add <file>` takes its filename immediately, so every other option has
# to come first. Each entry is "<mode>|<extra options>".
# --addprofile is deliberately absent, for the reason verify_e2e.sh gives: it
# requires the added block to form a monophyletic cluster in the guide tree,
# which an arbitrary split of an alignment usually does not.
ADDSTRATS=(
  "add|"
  "add|--keeplength"
  "add|--keeplength --mapout"
  "add|--localpair --maxiterate 0"
  "add|--localpair --keeplength --mapout --maxiterate 0"
  "add|--globalpair --maxiterate 0"
  "add|--allowshift --globalpair --maxiterate 0"
)
for f in "${INPUTS[@]}"; do
  n=$(grep -c '^>' "$f")
  [ "$n" -ge 4 ] || continue
  MAFFT_BINARIES="$CLIB" "$MAFFT" --retree 2 --maxiterate 0 "$f" > "$work/aligned.fa" 2>/dev/null || continue
  [ -s "$work/aligned.fa" ] || continue
  nadd=$(( n / 2 ))
  for mode in add; do
    python3 "$rust/gen/gen_add.py" "$work/aligned.fa" "$nadd" "$mode" > "$work/both.fa" || continue
    # Split it back into "the existing alignment" and "the sequences to add".
    python3 - "$work/both.fa" "$nadd" "$work/exist.fa" "$work/new.fa" <<'PYEOF'
import sys
recs, name, seq = [], None, []
for line in open(sys.argv[1]):
    if line.startswith('>'):
        if name is not None:
            recs.append((name, ''.join(seq)))
        name, seq = line.rstrip('\n'), []
    else:
        seq.append(line.strip())
recs.append((name, ''.join(seq)))
nadd = int(sys.argv[2])
for path, part in ((sys.argv[3], recs[:len(recs)-nadd]), (sys.argv[4], recs[len(recs)-nadd:])):
    with open(path, 'w') as fh:
        for nm, sq in part:
            fh.write(nm + '\n' + sq + '\n')
PYEOF
    for i in "${!ADDSTRATS[@]}"; do
      entry="${ADDSTRATS[$i]}"
      [ "${entry%%|*}" = "$mode" ] || continue
      extra="${entry#*|}"
      addopt="--add"
      rm -f "$work/new.fa.map"
      MAFFT_BINARIES="$CLIB" "$MAFFT" $extra $addopt "$work/new.fa" "$work/exist.fa" > "$work/o.c" 2>/dev/null
      cst=$?
      [ -f "$work/new.fa.map" ] && mv "$work/new.fa.map" "$work/map.c"
      rm -f "$work/new.fa.map"
      MAFFT_BINARIES="$PLIB" "$MAFFT" $extra $addopt "$work/new.fa" "$work/exist.fa" > "$work/o.p" 2>/dev/null
      pst=$?
      [ -f "$work/new.fa.map" ] && mv "$work/new.fa.map" "$work/map.p"
      label="$(basename "$f") [$addopt $extra, n=$nadd]"
      if [ $cst -ne 0 ] || [ ! -s "$work/o.c" ]; then
        fail=$((fail+1)); failed="$failed\n  $label (reference run failed: exit $cst)"
      elif [ $cst -ne $pst ]; then
        fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $pst)"
      elif ! cmp -s "$work/o.c" "$work/o.p"; then
        fail=$((fail+1)); failed="$failed\n  $label (output differs)"
      elif { [ -f "$work/map.c" ] || [ -f "$work/map.p" ]; } && ! cmp -s "$work/map.c" "$work/map.p"; then
        fail=$((fail+1)); failed="$failed\n  $label (.map differs)"
      else
        pass=$((pass+1))
      fi
      rm -f "$work/map.c" "$work/map.p"
    done
  done
done

# The `--auto` pass. Only the strategies `--auto` itself can select run here:
# the full STRATS list on a 790-sequence input would mean `--localpair
# --maxiterate 100`, which is L-INS-i on 790 sequences and does not finish in
# any time useful to a gate. Restricting the list is not weakening the check --
# these are the invocations the driver actually makes on these shapes, and the
# small-input pass above has already run the full list.
if [ "$SHAPES" = 1 ]; then
  shapedir="$work/shapes"
  gen_shapes "$shapedir"
  AUTOSTRATS=("--auto" "--retree 2 --maxiterate 0" "")
  for f in "$shapedir/shapes_123.fa" "$shapedir/shapes_790.fa"; do
    n=$(grep -c '^>' "$f")
    for i in "${!AUTOSTRATS[@]}"; do
      s="${AUTOSTRATS[$i]}"
      MAFFT_BINARIES="$CLIB" "$MAFFT" $s "$f" > "$work/o.c" 2>/dev/null
      cst=$?
      MAFFT_BINARIES="$PLIB" "$MAFFT" $s "$f" > "$work/o.p" 2>/dev/null
      pst=$?
      label="$(basename "$f") (n=$n) [$s]"
      if [ $cst -ne 0 ] || [ ! -s "$work/o.c" ]; then
        fail=$((fail+1)); failed="$failed\n  $label (reference run failed: exit $cst)"
      elif [ $pst -ne 0 ] && is_refusable "$s" && [ ! -s "$work/o.p" ]; then
        refused=$((refused+1)); refusals="$refusals\n  $label"
      elif [ $cst -ne $pst ]; then
        fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $pst)"
      elif ! cmp -s "$work/o.c" "$work/o.p"; then
        fail=$((fail+1)); failed="$failed\n  $label (output differs)"
      else
        pass=$((pass+1))
      fi
    done
  done
fi

echo "pkg: pass=$pass fail=$fail refused=$refused"
# Name them. A tally that says "refused=3" without saying which three lets a
# newly-refused strategy hide inside a number that has not changed.
if [ $refused -ne 0 ]; then
  printf 'REFUSED BY THE PORT (C succeeds, the port declines rather than approximating):%b\n' "$refusals"
fi

# A harness reporting pass=0 fail=0 is vacuous, not green. This repository has
# been fooled by that twice, which is why verify.sh counts its corpus and
# verify_e2e.sh refuses to run without a driver. Here the same trap is reached
# by passing an input the C cannot align, or a directory the driver never looks
# at: every loop runs zero times and the tally is silent about it.
if [ $pass -eq 0 ]; then
  echo "verify_pkg.sh: pass=0 -- nothing was actually compared. Refusing to report" 1>&2
  echo "               this as a result." 1>&2
  exit 2
fi

if [ $fail -ne 0 ]; then
  printf 'FAILED:%b\n' "$failed"
  exit 1
fi
