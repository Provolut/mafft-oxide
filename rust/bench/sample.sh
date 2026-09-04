#!/bin/bash
# Sampling profile of the port, cheap enough to run on every push.
#
#   ./sample.sh [budget-percent]        # default 8
#
# `hotspots.sh` uses callgrind, which counts every instruction exactly and
# costs about 50x. That is the right tool for deciding whether a change helped,
# and the wrong one for running on each push. `gprofng` samples instead:
# measured at +3.8% here, against the 8% budget this script enforces.
#
# The budget is checked, not assumed. Each case is timed with and without the
# profiler, and if the profiler costs more than the budget the profile is
# discarded and the script says so rather than recording a number taken under
# conditions that no longer resemble a normal run. That matters more than it
# sounds: profiling that distorts what it measures has already caught this
# repository out once, when line tables on the release profile read as a 9%
# memory regression.
#
# Exit status is 0 whenever the profile was produced or skipped for a reason
# the caller cannot fix (no gprofng). It is 1 only if the profiler is over
# budget, which is a real finding.

set -u
here=$(cd "$(dirname "$0")" && pwd)
BUDGET="${1:-8}"
RSDIR="$here/../target/release"
SAMPLE="$here/../../test/sample"
OUT="$here/sample-profile.txt"

if ! command -v gprofng >/dev/null; then
  echo "sample: gprofng not found, skipping the sampling profile" 1>&2
  exit 0
fi
[ -x "$RSDIR/disttbfast-rs" ] || { echo "sample: build the release binaries first" 1>&2; exit 0; }

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

D="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0"
TAIL="-X 0.1 -x 100"
TB_P="-u 0.0 -l 2.7 -C 0 -b 62 -g -0.100 -f -2.00 -Q 100.0 -h 0.1"
TB_T="-W 0.00001 -V -1.53 -s 0.0 -C 0 -b 62 -f -1.53 -Q 100.0 -h 0 -F -l 2.7 -X 0.1"

# Only the cases long enough to collect samples from. A 5 ms run yields none,
# and its overhead ratio is process startup rather than the profiler.
CASES=$(cat <<EOF
memsave	disttbfast-rs	$D -F $TAIL -M -B
linsi	tbfast-rs	_ $TB_P -L _ -+ 2 $TB_T -O
EOF
)

best_ms() {
  local best=9999999 s e ms i
  for i in 1 2 3; do
    s=$(date +%s%N); ( cd "$work" && "$@" >/dev/null 2>&1 ); e=$(date +%s%N)
    ms=$(( (e - s) / 1000000 ))
    [ "$ms" -lt "$best" ] && best=$ms
  done
  echo "$best"
}

# gprofng refuses to reuse an experiment directory, so a repeated run would
# fail in milliseconds and best-of-N would report it as the fastest -- which
# reads as the profiler making the program 99% faster. Each iteration gets its
# own directory, and the last one is kept for the report.
best_ms_sampled() {
  local dir="$1"; shift
  local best=9999999 s e ms i
  for i in 1 2 3; do
    rm -rf "$dir"
    s=$(date +%s%N)
    ( cd "$work" && gprofng collect app -o "$dir" "$@" >/dev/null 2>&1 )
    e=$(date +%s%N)
    [ -d "$dir" ] || { echo "gprofng produced no experiment" 1>&2; echo -1; return; }
    ms=$(( (e - s) / 1000000 ))
    [ "$ms" -lt "$best" ] && best=$ms
  done
  echo "$best"
}

status=0
{
  echo "# Sampling profile (gprofng), refreshed by bench/sample.sh"
  echo "# Overhead budget: ${BUDGET}%. Measured per case below."
  echo
} > "$OUT.tmp"

printf '%s\n' "$CASES" | while IFS=$'\t' read -r tag bin argv; do
  [ -n "$tag" ] || continue
  # shellcheck disable=SC2086
  set -- "$RSDIR/$bin" $argv -i "$SAMPLE"
  plain=$(best_ms "$@")
  sampled=$(best_ms_sampled "$work/$tag.er" "$@")
  if [ "$sampled" -lt 0 ]; then
    { echo "== $tag: gprofng failed, no profile"; echo; } >> "$OUT.tmp"
    continue
  fi
  over=$(awk -v a="$plain" -v b="$sampled" 'BEGIN{ printf "%.1f", a? (b/a-1)*100 : 0 }')

  {
    echo "== $tag: ${plain}ms plain, ${sampled}ms sampled, overhead ${over}%"
    if awk -v o="$over" -v b="$BUDGET" 'BEGIN{ exit !(o > b) }'; then
      echo "   OVER BUDGET -- profile discarded; the numbers a profiler this"
      echo "   expensive produces are not about the program any more."
      echo "over-budget" > "$work/failed"
    else
      gprofng display text -functions "$work/$tag.er" 2>/dev/null \
        | sed -n '/Excl\./,$p' | head -14 | sed 's/^/   /'
    fi
    echo
  } >> "$OUT.tmp"
done

[ -f "$work/failed" ] && status=1
mv "$OUT.tmp" "$OUT"
cat "$OUT"
exit $status
