#!/bin/bash
# End-to-end check: run the real `mafft` wrapper twice -- once against the
# stock libexec directory, once against a copy in which `disttbfast` has been
# replaced by `disttbfast-rs` -- and require byte-identical output.
#
#   verify_e2e.sh <libexec-dir-with-C-binaries> <path-to-disttbfast-rs> [inputs...]
#
# Set TBFAST_RS=<path-to-tbfast-rs> to swap that binary in as well, which
# brings the L-INS-i / G-INS-i strategies into the comparison.

set -u
LIBEXEC="${1:?usage: verify_e2e.sh <libexec-dir> <disttbfast-rs> [inputs...]}"
RBIN="${2:?usage: verify_e2e.sh <libexec-dir> <disttbfast-rs> [inputs...]}"
shift 2

here=$(cd "$(dirname "$0")" && pwd)
INPUTS=("$@")
if [ ${#INPUTS[@]} -eq 0 ]; then
  INPUTS=("$here/../test/sample")
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cp -r "$LIBEXEC" "$work/rs"
cp "$RBIN" "$work/rs/disttbfast"
chmod +x "$work/rs/disttbfast"
if [ -n "${TBFAST_RS:-}" ]; then
  cp "$TBFAST_RS" "$work/rs/tbfast"
  chmod +x "$work/rs/tbfast"
fi

# Strategies whose alignment work goes through the ported disttbfast path.
# NW-NS-i also runs dvtditr, which consumes disttbfast's output, so an
# identical result there exercises the port through the refinement stage too.
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
  # Consistency-based strategies: these run tbfast, whose pairwise stage and
  # constrained progressive stage are the other half of the port.
  "--localpair --maxiterate 0"
  "--localpair --maxiterate 2 --thread 0"
  "--localpair --maxiterate 100 --thread 0"
  "--globalpair --maxiterate 0"
  "--globalpair --maxiterate 2 --thread 0"
  "--globalpair --maxiterate 100 --thread 0"
  "--allowshift --globalpair --maxiterate 100 --thread 0"
  "--auto"
)

pass=0; fail=0; failed=""
for f in "${INPUTS[@]}"; do
  for i in "${!STRATS[@]}"; do
    s="${STRATS[$i]}"
    MAFFT_BINARIES="$LIBEXEC" mafft $s "$f" > "$work/o.c" 2>/dev/null
    cst=$?
    MAFFT_BINARIES="$work/rs" mafft $s "$f" > "$work/o.r" 2>/dev/null
    rst=$?
    label="$(basename "$f") [$s]"
    if [ $cst -ne $rst ]; then
      fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
    elif ! cmp -s "$work/o.c" "$work/o.r"; then
      fail=$((fail+1)); failed="$failed\n  $label (output differs)"
    else
      pass=$((pass+1))
    fi
  done
done

echo "pass=$pass fail=$fail"
if [ $fail -ne 0 ]; then
  printf 'FAILED:%b\n' "$failed"
  exit 1
fi
