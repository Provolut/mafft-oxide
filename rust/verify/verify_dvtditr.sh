#!/bin/bash
# Differential test: run the C `dvtditr` and `dvtditr-rs` on the same staged
# working directory and require the refined alignment to be identical.
#
#   verify_dvtditr.sh [--regenerate] <libexec-dir-with-C-binaries> <path-to-dvtditr-rs> [inputs...]
#
# With --regenerate, each staged case and the C's answers are written to
# tests/golden/dvtditr/ as well, which is what `tests/dvtditr.rs` replays on a
# machine with no C. Regenerate when the C changes or a case is added -- never
# to accommodate a diff in the port, which would be the port changing its
# answer.
#
# `dvtditr` is not a filter: it reads `pre` from stdin and writes it back, and
# reads `hat2` and `hat3` from the working directory the `mafft` script built.
# So each case is staged by running the real wrapper once with `dvtditr`
# replaced by a script that copies the working directory and its stdin aside --
# the §9 "wrap, trace, then port" method -- and both binaries are then replayed
# against copies of that directory with the argument vector the script really
# passed.
#
# Set MAFFT=<path-to-mafft-wrapper> if it is not next to the libexec directory.
set -u
REGEN=
if [ "${1:-}" = "--regenerate" ]; then REGEN=1; shift; fi
LIBEXEC="${1:?usage: verify_dvtditr.sh [--regenerate] <libexec-dir> <dvtditr-rs> [inputs...]}"
RBIN="${2:?usage: verify_dvtditr.sh <libexec-dir> <dvtditr-rs> [inputs...]}"
shift 2

here=$(cd "$(dirname "$0")" && pwd)
rust=$(cd "$here/.." && pwd)
repo=$(cd "$rust/.." && pwd)
abspath() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s/%s\n' "$(pwd)" "$1" ;; esac; }
LIBEXEC=$(abspath "$LIBEXEC")
RBIN=$(abspath "$RBIN")
[ -x "$RBIN" ] || { echo "verify_dvtditr: $RBIN is not executable" 1>&2; exit 2; }
[ -x "$LIBEXEC/dvtditr" ] || { echo "verify_dvtditr: no dvtditr in $LIBEXEC" 1>&2; exit 2; }

GOLDEN="$rust/tests/golden/dvtditr"
if [ -n "$REGEN" ]; then
  rm -rf "$GOLDEN"
  mkdir -p "$GOLDEN"
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# Two inputs the repo's own test data does not cover, generated here so nothing
# has to be committed:
#
#  * a wide alignment with three strongly conserved blocks, which is what makes
#    `searchAnchors` cut it into four segments. Without it the segment loop only
#    ever runs once, and the concatenation across segments is never exercised.
#  * a two-sequence input, which takes dvtditr's early exit -- no branch to
#    refine, so it writes its input straight back.
python3 - "$work/segments.fa" "$work/pair.fa" <<'GENPY'
import random, sys
random.seed(5)
AA = "ACDEFGHIKLMNPQRSTVWY"
blocks = ["".join(random.choice("WFYLIVM") for _ in range(90)) for _ in range(3)]
rows = []
for _ in range(8):
    parts = []
    for b in blocks:
        parts.append("".join(random.choice(AA) for _ in range(70)))
        parts.append(b)
    parts.append("".join(random.choice(AA) for _ in range(70)))
    rows.append("".join(parts))
with open(sys.argv[1], "w") as fh:
    for i, r in enumerate(rows):
        fh.write(">seg%d\n%s\n" % (i, r))
base = "".join(random.choice(AA) for _ in range(120))
other = list(base)
for p in range(0, 120, 7):
    other[p] = random.choice(AA)
with open(sys.argv[2], "w") as fh:
    fh.write(">a\n%s\n>b\n%s\n" % (base, "".join(other)))
GENPY

INPUTS=("$@")
if [ ${#INPUTS[@]} -eq 0 ]; then
  INPUTS=("$repo/test/sample" "$repo/test/samplerna" "$work/segments.fa" "$work/pair.fa")
fi

MAFFT="${MAFFT:-$repo/scripts/mafft}"
[ -x "$MAFFT" ] || { echo "verify_dvtditr: cannot find the mafft wrapper" 1>&2; exit 2; }

# A libexec copy whose dvtditr stages what it was called with, then runs the C.
cp -r "$LIBEXEC" "$work/stager"
mv "$work/stager/dvtditr" "$work/stager/dvtditr.real"
cat > "$work/stager/dvtditr" <<'WRAP'
#!/bin/bash
mkdir -p "$STAGE_DIR"
cp -r "$PWD"/. "$STAGE_DIR"/ 2>/dev/null
printf '%s\n' "$*" > "$STAGE_DIR/.argv"
cat > "$STAGE_DIR/.stdin"
exec "$(dirname "$0")/dvtditr.real" "$@" < "$STAGE_DIR/.stdin"
WRAP
chmod +x "$work/stager/dvtditr"

# The strategies that reach dvtditr at all. `--allowshift` is deliberately not
# here: it passes `-s 0.8`, the specificity scoring path dvtditr-rs refuses.
STRATS=(
  "--auto"
  "--retree 2 --maxiterate 2 --thread 0"
  "--retree 2 --maxiterate 5 --thread 0"
  "--localpair --maxiterate 2 --thread 0"
  "--localpair --maxiterate 16 --thread 0"
  "--globalpair --maxiterate 16 --thread 0"
  "--genafpair --maxiterate 16 --thread 0"
  "--nofft --retree 2 --maxiterate 2 --thread 0"
)

pass=0; fail=0; staged=0; failed=""
for f in "${INPUTS[@]}"; do
  for i in "${!STRATS[@]}"; do
    s="${STRATS[$i]}"
    label="$(basename "$f") [$s]"
    stage="$work/stage.$staged"
    rm -rf "$stage"
    STAGE_DIR="$stage" MAFFT_BINARIES="$work/stager" "$MAFFT" $s "$f" > /dev/null 2>/dev/null
    if [ ! -f "$stage/.argv" ]; then
      # This strategy never reaches dvtditr; nothing to compare.
      continue
    fi
    staged=$((staged+1))
    argv=$(cat "$stage/.argv")

    rm -rf "$work/c" "$work/r"
    cp -r "$stage" "$work/c"; cp -r "$stage" "$work/r"
    # The staged copy carries the `pre` of whatever ran before; delete it so a
    # comparison cannot pass on a file neither binary wrote.
    rm -f "$work/c/pre" "$work/r/pre"

    ( cd "$work/c" && "$LIBEXEC/dvtditr" $argv < .stdin > /dev/null 2>/dev/null )
    cst=$?
    ( cd "$work/r" && "$RBIN" $argv < .stdin > /dev/null 2>/dev/null )
    rst=$?

    if [ $cst -ne $rst ]; then
      fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
    elif [ ! -s "$work/c/pre" ]; then
      fail=$((fail+1)); failed="$failed\n  $label (the C wrote no pre)"
    elif ! cmp -s "$work/c/pre" "$work/r/pre"; then
      fail=$((fail+1)); failed="$failed\n  $label (refined alignment differs)"
    elif ! cmp -s "$work/c/order" "$work/r/order"; then
      fail=$((fail+1)); failed="$failed\n  $label (order differs)"
    elif ! cmp -s "$work/c/infile.tree" "$work/r/infile.tree"; then
      # The guide tree dvtditr writes; `mafft --treeout` hands it to the user.
      fail=$((fail+1)); failed="$failed\n  $label (infile.tree differs)"
    else
      pass=$((pass+1))
      # Only a case the C and the port already agree on is worth recording; a
      # store written from a failing run would pin the disagreement.
      # Which cases are recorded, and why not all of them. The store is
      # tracked and read on every checkout, so it buys breadth per kilobyte
      # rather than breadth for its own sake: every strategy on `pair.fa`,
      # which is the smallest input and covers each argument vector and the
      # two-sequence early exit, plus `--auto` on the others for the shapes
      # `pair.fa` cannot have -- the segment loop, and nucleotide input. The
      # `test/sample` cases are 300 KB each and add no path the others miss.
      # The wide 32-case differential against a live C is this script's
      # ordinary mode and is unaffected by the selection.
      if [ -n "$REGEN" ] && { [ "$(basename "$f")" = pair.fa ] || [ "$i" -eq 0 ]; } && \
         [ "$(basename "$f")" != sample ]; then
        d=$(printf '%s/%02d' "$GOLDEN" "$staged")
        mkdir -p "$d/in" "$d/out"
        # Everything the working directory held before the run. `pre` is
        # excluded on purpose: it is the output, and the harness deletes it
        # above so nothing can pass on a file neither binary wrote.
        for g in "$stage"/*; do
          b=$(basename "$g")
          [ "$b" = pre ] && continue
          [ -f "$g" ] && cp "$g" "$d/in/$b"
        done
        cp "$stage/.argv" "$d/argv"
        cp "$stage/.stdin" "$d/stdin"
        cp "$work/c/pre" "$d/out/pre"
        cp "$work/c/order" "$d/out/order"
        cp "$work/c/infile.tree" "$d/out/infile.tree"
        printf '%s\n' "$cst" > "$d/exit"
        printf '%s\n' "$label" > "$d/label"
      fi
    fi
  done
done

if [ "$staged" -lt 8 ]; then
  echo "verify_dvtditr: only $staged cases reached dvtditr; refusing" 1>&2
  exit 2
fi
echo "dvtditr: pass=$pass fail=$fail (of $staged staged runs)"
if [ $fail -ne 0 ]; then printf 'FAILED:%b\n' "$failed" | head -20; exit 1; fi
