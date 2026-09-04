#!/bin/bash
# Per-function instruction counts for one case, C and port side by side.
#
#   ./hotspots.sh [case-tag] [top-n]
#
# `profile.tsv` says *how much* slower the port is on a case; this says *where*.
# Both columns are callgrind `Ir`, so they are directly comparable and neither
# moves with machine load.
#
# The two lists do not line up name for name -- the port has no `A__align`,
# it has `AlignState::a_align`, and Rust's inlining differs from the C's -- so
# read them as two profiles of the same work rather than as a diff. What
# transfers is the shape: which phase dominates, and whether the port spends
# instructions on something the C does not do at all (allocation, bounds
# checks, memcpy of buffers the C keeps `static`).

set -u
here=$(cd "$(dirname "$0")" && pwd)
CORE="$here/../../core"
RSDIR="$here/../target/profiling"
SAMPLE="$here/../../test/sample"
MATRIX="$here/../tests/golden/matrix"

TAG="${1:-sample.fftns2}"
TOPN="${2:-20}"

D="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0"
TAIL="-X 0.1 -x 100"
TB_P="-u 0.0 -l 2.7 -C 0 -b 62 -g -0.100 -f -2.00 -Q 100.0 -h 0.1"
TB_T="-W 0.00001 -V -1.53 -s 0.0 -C 0 -b 62 -f -1.53 -Q 100.0 -h 0 -F -l 2.7 -X 0.1"

case "$TAG" in
  sample.fftns2)  BIN=disttbfast; IN=$SAMPLE;            ARGV="$D -F $TAIL" ;;
  sample.nwns2)   BIN=disttbfast; IN=$SAMPLE;            ARGV="$D $TAIL" ;;
  sample.memsave) BIN=disttbfast; IN=$SAMPLE;            ARGV="$D -F $TAIL -M -B" ;;
  sample.linsi)   BIN=tbfast;     IN=$SAMPLE;            ARGV="_ $TB_P -L _ -+ 2 $TB_T -O" ;;
  aa.fftns2)      BIN=disttbfast; IN=$MATRIX/aa.fa;      ARGV="$D -F $TAIL" ;;
  tandem.fftns2)  BIN=disttbfast; IN=$MATRIX/tandem.fa;  ARGV="$D -F $TAIL" ;;
  *) echo "hotspots: unknown case '$TAG'" 1>&2; exit 2 ;;
esac

command -v callgrind_annotate >/dev/null || {
  echo "hotspots: callgrind_annotate not found" 1>&2; exit 2; }

# The `profiling` profile is release codegen plus line tables. Building it
# here rather than expecting it to exist, because the binary it produces is
# not the one `profile.sh` measures.
( cd "$here/.." && cargo build --profile profiling --quiet ) || exit 2

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

run_one() {
  local label="$1" bin="$2"
  # shellcheck disable=SC2086
  ( cd "$work" && valgrind --tool=callgrind --callgrind-out-file="$work/$label.out" \
      "$bin" $ARGV -i "$IN" >/dev/null 2>/dev/null )
  echo
  echo "=== $label: top $TOPN functions by Ir ==="
  callgrind_annotate --threshold=99 "$work/$label.out" 2>/dev/null \
    | sed -n '/Ir *file:function/,$p' \
    | grep -E '^ *[0-9,]+' \
    | head -"$TOPN" \
    | sed 's/  */ /g'
}

echo "case: $TAG   binary: $BIN"
run_one c "$CORE/$BIN"
run_one rust "$RSDIR/$BIN-rs"
