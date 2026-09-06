#!/bin/bash
# Differential test: run `Falign_localhom` in the C and `falign_localhom` in
# the port on the same staged pair of groups and require bit-identical
# results.
#
#   verify_flh.sh [--regenerate] <path-to-core-dir>
#
# With --regenerate, the C's answers are written to tests/golden/flh/ instead
# of compared, so `cargo test` can replay them with no C present. Regenerate
# only when the C changes or a case is added -- never to accommodate a diff in
# the port, which would be the port changing its answer.
#
# The core directory is a built `core/`; its objects are linked into
# `gen_flhharness.c`. Scores are compared as raw f64 bit patterns, so a
# formatting difference cannot hide a value difference.
set -u
REGEN=0
if [ "${1:-}" = "--regenerate" ]; then REGEN=1; shift; fi
COREDIR="${1:?usage: verify_flh.sh [--regenerate] <core-dir>}"
here=$(cd "$(dirname "$0")" && pwd)
rust=$(cd "$here/.." && pwd)
repo=$(cd "$rust/.." && pwd)
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT

# Everything except the objects that carry a `main`.
objs=$(ls "$COREDIR"/*.o 2>/dev/null | grep -vE '/(disttbfast|tbfast|dvtditr|splittbfast|addsingle|pairlocalalign|pairlocalalignmain|multi2hat3s|nodepair|setcore|f2cl|countlen|version|replaceu|restoreu|dndpre|filter|score|getlag|seq2regtable|regtable2seq|setdirection|makedirectionlist|maffttext2hex|hex2maffttext|mccaskillwrap|contrafoldwrap|dndfast7|dndblast|sextet5|mafft-distance|pairash|mafft-profile)\.o$')
if [ -z "$objs" ]; then echo "verify_flh: no objects in $COREDIR" 1>&2; exit 2; fi

# `-Denablemultithread` matches how the objects were built: without it the
# `TLS` storage class of the shared globals disagrees and the link fails.
CC=${CC:-/usr/bin/gcc-13}
$CC -std=c99 -Denablemultithread -Denableatomic -O0 -g -I "$COREDIR" \
    -o "$work/flhharness" "$rust/gen/gen_flhharness.c" $objs -lm -lpthread || exit 2
if [ $REGEN -eq 0 ]; then
  ( cd "$rust" && cargo build --release --example flhref ) || exit 2
fi
RREF="$rust/target/release/examples/flhref"
if [ $REGEN -eq 0 ] && [ ! -x "$RREF" ]; then
  echo "verify_flh: $RREF is not executable" 1>&2; exit 2
fi
golden="$rust/tests/golden/flh"
if [ $REGEN -eq 1 ]; then mkdir -p "$golden" && rm -f "$golden"/*.args "$golden"/*.expected; fi

python3 "$rust/gen/gen_flhcases.py" > "$work/cases.txt"
n=$(wc -l < "$work/cases.txt")
# A harness that compares nothing reports pass=0 fail=0 and looks green; this
# floor is why it cannot.
if [ "$n" -lt 50 ]; then echo "verify_flh: only $n cases; refusing" 1>&2; exit 2; fi

pass=0; fail=0; failed=""
i=0
while IFS= read -r line; do
  i=$((i+1))
  IFS=$'\t' read -r -a args <<< "$line"
  "$work/flhharness" "${args[@]}" > "$work/c.txt" 2>/dev/null || {
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

# Single calls say nothing about the buffers a call leaves behind, and a real
# run makes hundreds of calls in one process. This stage replays a coherent run
# -- constant njob and width, as tditeration keeps them -- through one process
# on the C side and one state object on the port's.
python3 "$rust/gen/gen_flhcases.py" --sequence > "$work/seq.txt"
seqn=$(wc -l < "$work/seq.txt")
if [ "$seqn" -lt 20 ]; then echo "verify_flh: sequence has only $seqn calls; refusing" 1>&2; exit 2; fi
"$work/flhharness" --sequence < "$work/seq.txt" > "$work/seq_c.txt" 2>/dev/null || {
  echo "verify_flh: the C failed on the sequence" 1>&2; exit 1; }
if [ ! -s "$work/seq_c.txt" ]; then
  echo "verify_flh: the C printed nothing for the sequence" 1>&2; exit 1; fi
if [ $REGEN -eq 1 ]; then
  cp "$work/seq.txt" "$golden/sequence.cases"
  cp "$work/seq_c.txt" "$golden/sequence.expected"
  seqres="wrote $seqn-call sequence"
else
  "$RREF" --sequence < "$work/seq.txt" > "$work/seq_r.txt" 2>/dev/null || {
    echo "verify_flh: the port failed on the sequence" 1>&2; exit 1; }
  if cmp -s "$work/seq_c.txt" "$work/seq_r.txt"; then
    seqres="sequence of $seqn calls identical"
  else
    seqres="SEQUENCE DIFFERS"; fail=$((fail+1))
  fi
fi

if [ $REGEN -eq 1 ]; then
  echo "falign_localhom: wrote $pass fixtures to $golden, $seqres"
else
  echo "falign_localhom: pass=$pass fail=$fail, $seqres"
fi
if [ $fail -ne 0 ]; then printf 'FAILED:%b\n' "$failed" | head -20; exit 1; fi
