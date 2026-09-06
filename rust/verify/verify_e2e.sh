#!/bin/bash
# End-to-end check: run the real `mafft` wrapper twice -- once against the
# stock libexec directory, once against a copy in which `disttbfast` has been
# replaced by `disttbfast-rs` -- and require byte-identical output.
#
#   verify_e2e.sh <libexec-dir-with-C-binaries> <path-to-disttbfast-rs> [inputs...]
#
# Set TBFAST_RS=<path-to-tbfast-rs> to swap that binary in as well, which
# brings the L-INS-i / G-INS-i strategies into the comparison, and
# DVTDITR_RS=<path-to-dvtditr-rs> for the refinement stage. With all three set,
# every `--maxiterate` strategy below runs with no C binary in it at all.

set -u
LIBEXEC="${1:?usage: verify_e2e.sh <libexec-dir> <disttbfast-rs> [inputs...]}"
RBIN="${2:?usage: verify_e2e.sh <libexec-dir> <disttbfast-rs> [inputs...]}"
shift 2

here=$(cd "$(dirname "$0")" && pwd)

# The wrapper runs each alignment from a scratch directory of its own, so a
# relative `MAFFT_BINARIES` resolves to somewhere that does not exist and every
# reference run fails -- which this harness reports as 22 failures rather than
# as its own misuse. Absolute from here on.
case "$LIBEXEC" in /*) ;; *) LIBEXEC="$(pwd)/$LIBEXEC" ;; esac
case "$RBIN" in /*) ;; *) RBIN="$(pwd)/$RBIN" ;; esac
[ -d "$LIBEXEC" ] || { echo "verify_e2e.sh: no such libexec dir: $LIBEXEC" 1>&2; exit 2; }
[ -x "$RBIN" ] || { echo "verify_e2e.sh: not executable: $RBIN" 1>&2; exit 2; }

INPUTS=("$@")
if [ ${#INPUTS[@]} -eq 0 ]; then
  INPUTS=("$here/../test/sample")
fi

# The wrapper itself: $MAFFT if set, else the `mafft` next to the libexec dir,
# else whatever is on PATH. Refuse to run without one -- a missing wrapper makes
# both runs fail identically, which would otherwise be counted as a pass.
MAFFT="${MAFFT:-}"
if [ -z "$MAFFT" ]; then
  cand="$(cd "$LIBEXEC/../.." 2>/dev/null && pwd)/bin/mafft"
  if [ -x "$cand" ]; then MAFFT="$cand"; else MAFFT="$(command -v mafft || true)"; fi
fi
if [ ! -x "$MAFFT" ]; then
  echo "verify_e2e.sh: cannot find the mafft wrapper; set MAFFT=<path>" 1>&2
  exit 2
fi
case "$MAFFT" in /*) ;; *) MAFFT="$(pwd)/$MAFFT" ;; esac

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cp -r "$LIBEXEC" "$work/rs"
cp "$RBIN" "$work/rs/disttbfast"
chmod +x "$work/rs/disttbfast"
if [ -n "${TBFAST_RS:-}" ]; then
  cp "$TBFAST_RS" "$work/rs/tbfast"
  chmod +x "$work/rs/tbfast"
fi
if [ -n "${DVTDITR_RS:-}" ]; then
  cp "$DVTDITR_RS" "$work/rs/dvtditr"
  chmod +x "$work/rs/dvtditr"
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
  # E-INS-i: the generalised affine-gap pairwise aligner
  "--genafpair --maxiterate 0"
  "--genafpair --maxiterate 2 --thread 0"
  "--oldgenafpair --maxiterate 0"
  "--auto"
)

# `--allowshift` passes `dvtditr -s 0.8`, which selects the specificity scoring
# path (`scoringmatrices`, `classifypairs`, `partA__align_variousdist`) that
# `dvtditr-rs` refuses rather than approximates. With DVTDITR_RS set, that
# strategy is expected to be refused, and is counted separately -- not as a
# pass, which would hide it, and not as a failure, which would be untrue.
REFUSED_BY_DVTDITR_RS="--allowshift"

pass=0; fail=0; refused=0; failed=""
for f in "${INPUTS[@]}"; do
  for i in "${!STRATS[@]}"; do
    s="${STRATS[$i]}"
    MAFFT_BINARIES="$LIBEXEC" "$MAFFT" $s "$f" > "$work/o.c" 2>/dev/null
    cst=$?
    MAFFT_BINARIES="$work/rs" "$MAFFT" $s "$f" > "$work/o.r" 2>/dev/null
    rst=$?
    label="$(basename "$f") [$s]"
    if [ $cst -ne 0 ] || [ ! -s "$work/o.c" ]; then
      # The reference run itself failed, so the comparison would be vacuous.
      fail=$((fail+1)); failed="$failed\n  $label (reference run failed: exit $cst)"
    elif [ -n "${DVTDITR_RS:-}" ] && [ $rst -ne 0 ] \
         && [ "${s#*$REFUSED_BY_DVTDITR_RS}" != "$s" ]; then
      # Expected refusal, and only for the strategy named above: any other
      # non-zero exit still fails.
      refused=$((refused+1))
    elif [ $cst -ne $rst ]; then
      fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
    elif ! cmp -s "$work/o.c" "$work/o.r"; then
      fail=$((fail+1)); failed="$failed\n  $label (output differs)"
    else
      pass=$((pass+1))
    fi
  done
done

# --add / --addprofile / --keeplength take two files, so they need their own
# loop: align with the C wrapper first, then split the result and add the tail
# back. `--keeplength` also writes a `.map` file next to the input, which is
# part of the output and is compared too.
# `--add <file>` takes its filename immediately, so every other option has to
# come first. Each entry is "<mode>|<extra options>".
# --addprofile is deliberately absent: it requires the added block to form a
# monophyletic cluster in the guide tree, which an arbitrary split of an
# alignment usually does not, so the C refuses most cases here. verify.sh and
# verify_tbfast.sh cover it at the binary level instead.
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
  MAFFT_BINARIES="$LIBEXEC" "$MAFFT" --retree 2 --maxiterate 0 "$f" > "$work/aligned.fa" 2>/dev/null || continue
  [ -s "$work/aligned.fa" ] || continue
  nadd=$(( n / 2 ))
  for mode in add; do
    python3 "$here/gen_add.py" "$work/aligned.fa" "$nadd" "$mode" > "$work/both.fa" || continue
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
      MAFFT_BINARIES="$LIBEXEC" "$MAFFT" $extra $addopt "$work/new.fa" "$work/exist.fa" > "$work/o.c" 2>/dev/null
      cst=$?
      [ -f "$work/new.fa.map" ] && mv "$work/new.fa.map" "$work/map.c"
      rm -f "$work/new.fa.map"
      MAFFT_BINARIES="$work/rs" "$MAFFT" $extra $addopt "$work/new.fa" "$work/exist.fa" > "$work/o.r" 2>/dev/null
      rst=$?
      [ -f "$work/new.fa.map" ] && mv "$work/new.fa.map" "$work/map.r"
      label="$(basename "$f") [$addopt $extra, n=$nadd]"
      if [ $cst -ne 0 ] || [ ! -s "$work/o.c" ]; then
        fail=$((fail+1)); failed="$failed\n  $label (reference run failed: exit $cst)"
      elif [ $cst -ne $rst ]; then
        fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
      elif ! cmp -s "$work/o.c" "$work/o.r"; then
        fail=$((fail+1)); failed="$failed\n  $label (output differs)"
      elif { [ -f "$work/map.c" ] || [ -f "$work/map.r" ]; } && ! cmp -s "$work/map.c" "$work/map.r"; then
        fail=$((fail+1)); failed="$failed\n  $label (.map differs)"
      else
        pass=$((pass+1))
      fi
      rm -f "$work/map.c" "$work/map.r"
    done
  done
done

echo "e2e: pass=$pass fail=$fail refused=$refused"
if [ $fail -ne 0 ]; then
  printf 'FAILED:%b\n' "$failed"
  exit 1
fi
