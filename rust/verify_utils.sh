#!/bin/bash
# Differential test for the small utilities: run the C `countlen` / `version`
# and their ports on the same inputs and require byte-identical stdout, stderr
# and exit status.
#
#   verify_utils.sh <path-to-core-dir> [corpus-dir]
#   verify_utils.sh --regenerate <path-to-core-dir>
#
# With no corpus directory, the committed fixtures in tests/golden/utils/in are
# used, plus a set of generated edge cases. Point it at a directory of FASTA
# files to widen it -- that is the whole reason this exists next to
# tests/utils.rs, which pins the same behaviour without needing a C at all.
#
# `--regenerate` rewrites tests/golden/utils/*.expected from the C. Use it when
# the C changes or a case is added, never to accommodate a diff in the port:
# a port whose output moved is the finding, not the fixture.

set -u
here=$(cd "$(dirname "$0")" && pwd)

if [ "${1:-}" = "--regenerate" ]; then
  REGEN=1; shift
else
  REGEN=0
fi

CORE="${1:?usage: verify_utils.sh [--regenerate] <core-dir> [corpus-dir]}"
CORPUS="${2:-}"
RSDIR="$here/target/release"
GOLDEN="$here/tests/golden/utils"

for b in "$CORE/countlen" "$CORE/version" "$RSDIR/countlen-rs" "$RSDIR/version-rs"; do
  [ -x "$b" ] || { echo "verify_utils: $b is not executable" 1>&2; exit 2; }
done

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# Edge cases the committed corpus cannot carry as files: an empty input, a
# record with no trailing newline, and an input past the 100000-site sampling
# cap, which is where the C's N tally starts carrying the last sampled record
# forward (see issue #8 -- the port reproduces it deliberately).
gen_edge_cases() {
  local d="$1"
  printf '' > "$d/edge_empty.fa"
  printf '>a\nMKVLTA' > "$d/edge_nonewline.fa"
  printf '>a\n\n>b\nMKV\n' > "$d/edge_blankrecord.fa"
  printf '>a\n----\n>b\nMKVL\n' > "$d/edge_allgap.fa"
  printf '>a\nACGTNNNNACGT\n>b\nACGTACGTACGT\n' > "$d/edge_dna_n.fa"
  awk 'BEGIN{for(i=0;i<400;i++){printf(">s%d\n",i); for(j=0;j<400;j++) printf("%s","A"); printf("\n")}}' \
    > "$d/edge_over_cap.fa"
}

inputs="$work/inputs"
mkdir -p "$inputs"
cp "$GOLDEN"/in/* "$inputs/" 2>/dev/null
gen_edge_cases "$inputs"
if [ -n "$CORPUS" ]; then
  n=0
  for f in "$CORPUS"/*; do
    [ -f "$f" ] || continue
    cp "$f" "$inputs/corpus_$(basename "$f")"
    n=$((n+1))
  done
  echo "verify_utils: added $n file(s) from $CORPUS"
fi

if [ "$REGEN" = "1" ]; then
  : > "$GOLDEN/countlen.expected"
  for f in $(ls "$GOLDEN/in" | sort); do
    printf '%s\t%s\n' "$f" "$("$CORE/countlen" -i "$GOLDEN/in/$f" </dev/null 2>/dev/null)" \
      >> "$GOLDEN/countlen.expected"
  done
  "$CORE/version" > "$GOLDEN/version.expected"
  echo "verify_utils: regenerated $GOLDEN/*.expected from $CORE"
  exit 0
fi

pass=0; fail=0; indet=0

# Inputs whose *reference* output is undefined, so there is nothing to be
# bit-exact against. `countATGCandN` returns early for an empty sequence body
# and sets `*total` without setting `*countN`, so `getnumlen_nogap_countn`'s
# uninitialised `nN` reaches the `nfreq` division -- valgrind names the stack
# allocation, and the value it finds is not a contract (issue #8). The port
# initialises to zero and so prints `-nan` where this C happens to print `inf`.
#
# Listed here rather than dropped from the corpus: the run still executes them
# and still reports what each side said, it just does not score the comparison.
# Anything not on this list is a real failure.
INDETERMINATE="edge_blankrecord.fa"

is_indeterminate() {
  case " $INDETERMINATE " in *" $1 "*) return 0 ;; *) return 1 ;; esac
}

# Compare one argument vector on both binaries: stdout, stderr and status.
cmp_run() {
  local label="$1" cbin="$2" rbin="$3"; shift 3
  local co ce cs ro re rs
  co=$("$cbin" "$@" </dev/null 2>"$work/ce"); cs=$?
  ro=$("$rbin" "$@" </dev/null 2>"$work/re"); rs=$?
  ce=$(cat "$work/ce"); re=$(cat "$work/re")
  if [ "$co" = "$ro" ] && [ "$ce" = "$re" ] && [ "$cs" = "$rs" ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL $label"
    [ "$co" = "$ro" ] || { echo "  stdout C : $co"; echo "  stdout RS: $ro"; }
    [ "$ce" = "$re" ] || { echo "  stderr C : $ce"; echo "  stderr RS: $re"; }
    [ "$cs" = "$rs" ] || echo "  status C: $cs  RS: $rs"
  fi
}

for f in "$inputs"/*; do
  base=$(basename "$f")
  if is_indeterminate "$base"; then
    indet=$((indet+1))
    echo "INDETERMINATE $base -- C: $("$CORE/countlen" -i "$f" </dev/null 2>/dev/null) | RS: $("$RSDIR/countlen-rs" -i "$f" </dev/null 2>/dev/null)"
    continue
  fi
  cmp_run "countlen -i $base" "$CORE/countlen" "$RSDIR/countlen-rs" -i "$f"
done

# The stdin path is a separate code path in both, and the script uses it twice.
for f in "$inputs"/*; do
  is_indeterminate "$(basename "$f")" && continue
  co=$("$CORE/countlen" < "$f" 2>/dev/null); cs=$?
  ro=$("$RSDIR/countlen-rs" < "$f" 2>/dev/null); rs=$?
  if [ "$co" = "$ro" ] && [ "$cs" = "$rs" ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL countlen < $(basename "$f")"
    echo "  C : $co (status $cs)"
    echo "  RS: $ro (status $rs)"
  fi
done

# Argument handling, where the C's loop is easy to get subtly wrong.
cmp_run "countlen -Z"            "$CORE/countlen" "$RSDIR/countlen-rs" -Z
cmp_run "countlen -Zi file"      "$CORE/countlen" "$RSDIR/countlen-rs" -Zi "$inputs/edge_dna_n.fa"
cmp_run "countlen positional"    "$CORE/countlen" "$RSDIR/countlen-rs" "$inputs/edge_dna_n.fa"
cmp_run "countlen -i (no arg)"   "$CORE/countlen" "$RSDIR/countlen-rs" -i
cmp_run "countlen -i f extra"    "$CORE/countlen" "$RSDIR/countlen-rs" -i "$inputs/edge_dna_n.fa" extra
cmp_run "countlen --"            "$CORE/countlen" "$RSDIR/countlen-rs" --
cmp_run "countlen -i /nope"      "$CORE/countlen" "$RSDIR/countlen-rs" -i /nope
cmp_run "version"                "$CORE/version"  "$RSDIR/version-rs"

# A run that compared nothing is vacuous, not green.
if [ "$pass" -eq 0 ]; then
  echo "verify_utils: pass=0 -- nothing was compared, which is a harness failure" 1>&2
  exit 2
fi

echo "pass=$pass fail=$fail indeterminate=$indet"
[ "$fail" -eq 0 ] || exit 1
