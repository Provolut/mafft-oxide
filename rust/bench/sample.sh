#!/bin/bash
# Sampling profile of the port, cheap enough to run on every push.
#
#   ./sample.sh [budget-percent]        # default 8
#
# `hotspots.sh` uses callgrind, which counts every instruction exactly and
# costs about 50x. That is the right tool for deciding whether a change helped,
# and the wrong one for running on each push. `gprofng` samples instead, and
# the sampling itself measures at under 1% here, against the 8% budget this
# script enforces.
#
# The budget is checked, not assumed. Each case is timed with and without the
# profiler, and if the profiler costs more than the budget the profile is
# discarded and the script says so rather than recording a number taken under
# conditions that no longer resemble a normal run. That matters more than it
# sounds: profiling that distorts what it measures has already caught this
# repository out once, when line tables on the release profile read as a 9%
# memory regression.
#
# What counts as "the profiler costing something" is the subtle part; see the
# note above `calibrate` below. Briefly: `gprofng` also costs a flat ~13ms per
# run to start its collector, that cost is the same on a 130ms case as on a
# 570ms one, and it distorts nothing -- so it is measured and subtracted rather
# than charged to the program.
#
# A short case limits how finely any of this can be judged. Where the fixed
# cost is a noticeable share of the case, the report says so instead of
# implying a precision the measurement does not have.
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

# What the budget is actually about, and why the ratio alone will not do.
#
# `gprofng collect` costs a fixed amount per run -- forking the collector,
# creating the experiment directory, writing its headers, waiting for the flush
# at exit. Measured here that is ~11 ms, and it is ~11 ms whether the program
# then runs for 130 ms or 570 ms. It also does not distort anything: it happens
# around the program rather than inside it, so no sample lands in it and no
# function's share moves because of it.
#
# Charging it to the program is what broke this check. The budget is meant to
# ask "is the profile still about the program?", but a plain sampled/plain ratio
# also answers "how short is the case?" -- and it got the wrong answer here for
# the best possible reason: `tbfast` got 2.6x faster, the linsi case fell to
# ~128 ms, and the same unchanged 11 ms became 9% of it. Every push then warned
# and threw the profile away.
#
# So measure the fixed cost against a run that does nothing, take it off the
# top, and let the budget see only what is left -- the sampling distortion,
# which is the thing worth refusing a profile over.
calibrate() {
  local plain sampled c
  plain=$(best_ms "$RSDIR/version-rs")
  sampled=$(best_ms_sampled "$work/calibrate.er" "$RSDIR/version-rs")
  [ "$sampled" -lt 0 ] && { echo 0; return; }
  c=$(( sampled - plain ))
  # Never let noise turn into a credit against the cases.
  [ "$c" -lt 0 ] && c=0
  echo "$c"
}
FIXED_MS=$(calibrate)

status=0
{
  echo "# Sampling profile (gprofng), refreshed by bench/sample.sh"
  echo "# Overhead budget: ${BUDGET}%, counting only the sampling itself."
  echo "# The collector's fixed cost per run, measured against version-rs just"
  echo "# now, is ${FIXED_MS}ms; it is subtracted before the budget is applied,"
  echo "# because it happens around the program rather than inside it."
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
  over=$(awk -v a="$plain" -v b="$sampled" -v c="$FIXED_MS" \
    'BEGIN{ d = b - c - a; if( d < 0 ) d = 0; printf "%.1f", a? d / a * 100 : 0 }')

  {
    echo "== $tag: ${plain}ms plain, ${sampled}ms sampled (${FIXED_MS}ms of that is" \
         "the collector starting up), sampling overhead ${over}%"
    # Run-to-run the fixed cost moves by a millisecond or two. On a long case
    # that is nothing; on a short one it is most of the answer, and quoting a
    # tenth of a percent from it would be inventing precision.
    if awk -v c="$FIXED_MS" -v a="$plain" 'BEGIN{ exit !(a && c * 100 / a > 5) }'; then
      res=$(awk -v a="$plain" 'BEGIN{ printf "%.0f", 200 / a }')
      echo "   NOTE -- at ${plain}ms this case is short enough that the fixed cost is" \
           "$(awk -v c="$FIXED_MS" -v a="$plain" 'BEGIN{ printf "%.0f", c * 100 / a }')% of it;" \
           "the figure above resolves to about +-${res}%."
    fi
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
