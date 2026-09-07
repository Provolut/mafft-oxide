#!/bin/bash
# Differential test: the AVX-512 kernels against the scalar paths they replace.
#
#   verify_lalign_kernel.sh [inputs...]
#
# `l_align11` and `partA__align` each have two implementations now, and the
# vector ones are the reason the port carries any `unsafe` at all. This runs the
# same alignments twice through the same binary -- once with
# `MAFFT_RS_NO_AVX512=1`, which makes `lalign_avx512::available()` answer false
# and so switches off both kernels -- and requires byte-identical output. It is
# the check the other harnesses cannot make: they compare the port against the
# C, so they would pass if *both* of the port's paths were wrong in the same
# way.
#
# Both kernels are covered by the strategy list below. The `--localpair
# --maxiterate` runs and `--auto` reach the iterative refinement stage, which is
# where `dvtditr` and `partalign_avx512` run; the pairwise stage of every
# strategy reaches `l_align11`. The script keeps its name because `README.md`
# and `bench/OPTIMISATION-LOG.md` refer to it by that name.
#
# `pre-push` does NOT run this -- it runs the two golden replays, which compare
# against the C and so cannot tell the port's two paths apart. Run it by hand
# when either kernel changes.
#
# On a CPU without AVX-512 the two runs take the same path and the comparison is
# vacuous, so that case is reported rather than counted as a pass.
set -u
here=$(cd "$(dirname "$0")" && pwd)
rust=$(cd "$here/.." && pwd)
repo=$(cd "$rust/.." && pwd)
MAFFT="${MAFFT:-$repo/scripts/mafft}"
[ -x "$MAFFT" ] || { echo "verify_lalign_kernel: cannot find the mafft wrapper" 1>&2; exit 2; }

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# A libexec of the Rust binaries.
cp -r "$repo/binaries" "$work/rs"
for b in countlen version f2cl tbfast disttbfast dvtditr dndpre replaceu restoreu; do
  [ -f "$rust/target/release/$b-rs" ] && cp "$rust/target/release/$b-rs" "$work/rs/$b"
done

if ! "$rust/target/release/tbfast-rs" --help 2>&1 | head -0; then :; fi
# Does this machine even have the kernel? Ask the binary, not /proc.
if ! grep -qE '^flags.*\bavx512f\b' /proc/cpuinfo || ! grep -qE '^flags.*\bavx512vl\b' /proc/cpuinfo; then
  echo "verify_lalign_kernel: this CPU has no avx512f+avx512vl; both runs would take"
  echo "verify_lalign_kernel: the scalar path, so there is nothing to compare. SKIPPED."
  exit 0
fi

INPUTS=("$@")
if [ ${#INPUTS[@]} -eq 0 ]; then
  INPUTS=("$repo/test/sample" "$repo/test/samplerna")
  # A wide, deeply conserved input, which is where the kernel's row is longest
  # and its tail shortest.
  python3 - "$work/wide.fa" <<'GENPY'
import random, sys
random.seed(17)
AA = "ACDEFGHIKLMNPQRSTVWY"
base = "".join(random.choice(AA) for _ in range(600))
rows = []
for k in range(12):
    r = list(base)
    for p in range(0, 600, 11 + k):
        r[p] = random.choice(AA)
    rows.append("".join(r))
with open(sys.argv[1], "w") as fh:
    for i, r in enumerate(rows):
        fh.write(">w%d\n%s\n" % (i, r))
GENPY
  INPUTS+=("$work/wide.fa")
fi

# Strategies whose pairwise stage runs L__align11, and -- for those that
# iterate -- whose refinement stage runs partA__align.
STRATS=(
  "--auto"
  "--localpair --maxiterate 0"
  "--localpair --maxiterate 2 --thread 0"
  "--localpair --maxiterate 16 --thread 0"
  "--auto --preservecase"
)

pass=0; fail=0; failed=""
for f in "${INPUTS[@]}"; do
  for s in "${STRATS[@]}"; do
    label="$(basename "$f") [$s]"
    MAFFT_BINARIES="$work/rs" "$MAFFT" $s "$f" > "$work/on.fa" 2>/dev/null
    on=$?
    MAFFT_RS_NO_AVX512=1 MAFFT_BINARIES="$work/rs" "$MAFFT" $s "$f" > "$work/off.fa" 2>/dev/null
    off=$?
    if [ $on -ne 0 ] || [ $off -ne 0 ]; then
      fail=$((fail+1)); failed="$failed\n  $label (exit $on vs $off)"
    elif [ ! -s "$work/on.fa" ]; then
      fail=$((fail+1)); failed="$failed\n  $label (no output)"
    elif cmp -s "$work/on.fa" "$work/off.fa"; then
      pass=$((pass+1))
    else
      fail=$((fail+1)); failed="$failed\n  $label (kernel and scalar differ)"
    fi
  done
done

echo "avx512 kernels vs scalar: pass=$pass fail=$fail"
if [ $fail -ne 0 ]; then printf 'FAILED:%b\n' "$failed" | head -20; exit 1; fi
