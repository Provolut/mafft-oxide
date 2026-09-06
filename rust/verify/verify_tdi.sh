#!/bin/bash
# Differential test: run `TreeDependentIteration` in the C and
# `tditeration::tree_dependent_iteration` in the port on the same staged
# alignment, tree and constraint table, and require the refined alignments to
# be identical.
#
# This one compares a whole driver rather than a function: hundreds of
# accept/reject decisions, each resting on a score, feed the next branch's
# input. A single differing bit anywhere shows up as a different alignment.
#
#   verify_tdi.sh [--regenerate] <path-to-core-dir>
#
# With --regenerate, the C's answers are written to tests/golden/tdi/ instead
# of compared, so `cargo test` can replay them with no C present. Regenerate
# only when the C changes or a case is added -- never to accommodate a diff in
# the port, which would be the port changing its answer.
#
# The core directory is a built `core/`; its objects are linked into
# `gen_tdiharness.c`. Scores are compared as raw f64 bit patterns, so a
# formatting difference cannot hide a value difference.
set -u
REGEN=0
if [ "${1:-}" = "--regenerate" ]; then REGEN=1; shift; fi
COREDIR="${1:?usage: verify_tdi.sh [--regenerate] <core-dir>}"
here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT

# Everything except the objects that carry a `main`.
objs=$(ls "$COREDIR"/*.o 2>/dev/null | grep -vE '/(disttbfast|tbfast|dvtditr|splittbfast|addsingle|pairlocalalign|pairlocalalignmain|multi2hat3s|nodepair|setcore|f2cl|countlen|version|replaceu|restoreu|dndpre|filter|score|getlag|seq2regtable|regtable2seq|setdirection|makedirectionlist|maffttext2hex|hex2maffttext|mccaskillwrap|contrafoldwrap|dndfast7|dndblast|sextet5|mafft-distance|pairash|mafft-profile)\.o$')
if [ -z "$objs" ]; then echo "verify_tdi: no objects in $COREDIR" 1>&2; exit 2; fi

# `-Denablemultithread` matches how the objects were built: without it the
# `TLS` storage class of the shared globals disagrees and the link fails.
CC=${CC:-/usr/bin/gcc-13}
$CC -std=c99 -Denablemultithread -Denableatomic -O0 -g -I "$COREDIR" \
    -o "$work/tdiharness" "$here/gen_tdiharness.c" $objs -lm -lpthread || exit 2
if [ $REGEN -eq 0 ]; then
  ( cd "$here" && cargo build --release --example tdiref ) || exit 2
fi
RREF="$here/target/release/examples/tdiref"
if [ $REGEN -eq 0 ] && [ ! -x "$RREF" ]; then
  echo "verify_tdi: $RREF is not executable" 1>&2; exit 2
fi
golden="$here/tests/golden/tdi"
if [ $REGEN -eq 1 ]; then mkdir -p "$golden" && rm -f "$golden"/*.args "$golden"/*.expected; fi

# The committed fixtures hold the small trees only; the live differential runs
# the whole corpus. See gen_tdicases.py.
# The committed fixtures hold the small runs only; the live differential runs
# the whole corpus. See gen_tdicases.py.
if [ $REGEN -eq 1 ]; then
  python3 "$here/gen_tdicases.py" --small > "$work/cases.txt"
else
  python3 "$here/gen_tdicases.py" > "$work/cases.txt"
fi
n=$(wc -l < "$work/cases.txt")
# A harness that compares nothing reports pass=0 fail=0 and looks green; this
# floor is why it cannot.
if [ "$n" -lt 50 ]; then echo "verify_tdi: only $n cases; refusing" 1>&2; exit 2; fi

cd "$work" || exit 2
pass=0; fail=0; failed=""
i=0
while IFS= read -r line; do
  i=$((i+1))
  IFS=$'\t' read -r -a args <<< "$line"
  "$work/tdiharness" "${args[@]}" > "$work/c.txt" 2>/dev/null || {
    fail=$((fail+1)); failed="$failed\n  case $i (C failed)"; continue; }
  if [ ! -s "$work/c.txt" ]; then
    fail=$((fail+1)); failed="$failed\n  case $i (C printed nothing)"; continue; fi
  if [ $REGEN -eq 1 ]; then
    key=$(printf '%03d' "$i")
    printf '%s\n' "$line" > "$golden/$key.args"
    cp "$work/c.txt" "$golden/$key.expected"
    pass=$((pass+1))
    continue
  fi
  "$RREF" "${args[@]}" > "$work/r.txt" 2>/dev/null || {
    fail=$((fail+1)); failed="$failed\n  case $i (port failed)"; continue; }
  if cmp -s "$work/c.txt" "$work/r.txt"; then pass=$((pass+1)); else
    fail=$((fail+1)); failed="$failed\n  case $i"
  fi
done < "$work/cases.txt"

if [ $REGEN -eq 1 ]; then
  echo "refinement: wrote $pass fixtures to $golden"
else
  echo "refinement: pass=$pass fail=$fail"
fi
if [ $fail -ne 0 ]; then printf 'FAILED:%b\n' "$failed" | head -20; exit 1; fi
