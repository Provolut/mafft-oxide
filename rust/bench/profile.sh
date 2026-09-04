#!/bin/bash
# Build or refresh `profile.tsv`: the C and the port side by side, per
# (input, argument vector), on compute and memory.
#
#   ./profile.sh [core-dir] [rust-target-dir] > profile.tsv
#
# Compute is callgrind's `Ir` -- instructions retired -- rather than wall clock,
# because it is deterministic: the same binary on the same input gives the same
# number on a loaded machine, so a diff in this file is a real change and not
# noise. Wall clock is recorded too, best of five, for the cases where the two
# disagree (they do: instruction count overstates the port's gap because its
# extra work is cache-friendly).
#
# Memory is peak RSS read from /proc/<pid>/VmHWM while the process runs, which
# is what the caller actually has to provision. massif and DHAT are available
# for explaining a row; they are not run here because they are slow and their
# output is not one number.
#
# Runtime is a few minutes: callgrind costs roughly 50x, so the inputs are the
# small ones from the golden matrix plus `test/sample` for a realistic size.

set -u
here=$(cd "$(dirname "$0")" && pwd)
CORE="${1:-$here/../../core}"
RSDIR="${2:-$here/../target/release}"
MATRIX="$here/../tests/golden/matrix"
SAMPLE="$here/../../test/sample"

for b in "$CORE/disttbfast" "$CORE/tbfast" "$RSDIR/disttbfast-rs" "$RSDIR/tbfast-rs"; do
  [ -x "$b" ] || { echo "profile: $b is not executable" 1>&2; exit 2; }
done
command -v valgrind >/dev/null || { echo "profile: valgrind not found" 1>&2; exit 2; }

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

D="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0"
TAIL="-X 0.1 -x 100"
TB_P="-u 0.0 -l 2.7 -C 0 -b 62 -g -0.100 -f -2.00 -Q 100.0 -h 0.1"
TB_T="-W 0.00001 -V -1.53 -s 0.0 -C 0 -b 62 -f -1.53 -Q 100.0 -h 0 -F -l 2.7 -X 0.1"

# tag <TAB> binary <TAB> input <TAB> argv
CASES=$(cat <<EOF
sample.fftns2	disttbfast	$SAMPLE	$D -F $TAIL
sample.nwns2	disttbfast	$SAMPLE	$D $TAIL
sample.memsave	disttbfast	$SAMPLE	$D -F $TAIL -M -B
sample.linsi	tbfast	$SAMPLE	_ $TB_P -L _ -+ 2 $TB_T -O
aa.fftns2	disttbfast	$MATRIX/aa.fa	$D -F $TAIL
nuc.fftns2	disttbfast	$MATRIX/nuc.fa	$D -F $TAIL
tandem.fftns2	disttbfast	$MATRIX/tandem.fa	$D -F $TAIL
dup.nwns2	disttbfast	$MATRIX/dup.fa	$D $TAIL
EOF
)

# Peak RSS in KB, by polling VmHWM: the kernel's own high-water mark, so the
# sampling interval does not matter.
#
# Reported as the median of several runs, not one. A single sample of this once
# read as a 7.6% memory regression from a change that only ever removed
# allocations, and as "no change" from one that removed 7 MB of arrays.
#
# The spread that caused it -- 31% on `sample.linsi`, clustered at 11.0, 13.4
# and 14.4 MB -- was the harness, not the program: every run shared one working
# directory, and these binaries write `hat2`, `hat3`, `order` and `pre` there
# and read them back. Each run now gets its own empty directory (see
# `scratch`), and the same binary then repeats inside 1%. The median is kept
# because it costs nothing and a stale conclusion costs a day.
# Every measured run gets an empty working directory of its own. Sharing one
# made the numbers meaningless: `tbfast` and `disttbfast` write `hat2`, `hat3`,
# `order` and `pre` into the cwd and read them back, so the second run of a
# case was not doing the same work as the first, and later cases inherited
# earlier ones' files. That, not glibc, is where the 31% spread in peak RSS
# came from -- in a fresh directory the same binary on the same input repeats
# inside 1%.
scratch() { mktemp -d "$work/run.XXXXXX"; }

peak_rss_once() {
  local d; d=$(scratch)
  ( cd "$d" && exec "$@" ) >/dev/null 2>&1 &
  local pid=$! peak=0 v
  while kill -0 "$pid" 2>/dev/null; do
    v=$(awk '/VmHWM/{print $2}' "/proc/$pid/status" 2>/dev/null)
    [ -n "${v:-}" ] && [ "$v" -gt "$peak" ] && peak=$v
    sleep 0.02
  done
  wait "$pid" 2>/dev/null
  rm -rf "$d"
  echo "$peak"
}

peak_rss() {
  local vals i
  vals=$(for i in 1 2 3 4 5; do peak_rss_once "$@"; done | sort -n)
  echo "$vals" | awk 'NR==3'
}

ir_count() {
  local out="$work/cg.$$" d; d=$(scratch)
  ( cd "$d" && valgrind --tool=callgrind --callgrind-out-file="$out" "$@" ) \
    >/dev/null 2>"$work/cg.log"
  grep -oP 'refs:\s+\K[0-9,]+' "$work/cg.log" | tr -d ,
  rm -f "$out"; rm -rf "$d"
}

# A run that failed is not a fast run. Without this a crashing port records as
# "0.000x", which is how a profiling database turns into a lie: check both
# binaries exit 0 and produce the same bytes before timing either of them.
# `disttbfast` and `tbfast` write different things. `tbfast` puts the alignment
# in `pre`; `disttbfast` puts it on stdout and separately creates an *empty*
# `pre` plus a `trace` log, from `initFiles()` in io.c -- neither is output, and
# `trace` starts with the pid, so it could not be reproduced even in principle.
# Comparing those two files is what made every disttbfast row read "differs".
outputs_for() {
  case "$1" in
    tbfast) echo "pre order hat2 hat3" ;;
    *)      echo "stdout.fa order hat2 infile.tree" ;;
  esac
}

check_pair() {
  local bin="$1" cbin="$2" rbin="$3"; shift 3
  local cd rd
  cd=$(mktemp -d); rd=$(mktemp -d)
  ( cd "$cd" && "$cbin" "$@" > stdout.fa 2>/dev/null ); local cst=$?
  ( cd "$rd" && "$rbin" "$@" > stdout.fa 2>/dev/null ); local rst=$?
  local verdict=ok
  if [ $cst -ne 0 ] || [ $rst -ne 0 ]; then
    verdict="exit:$cst/$rst"
  else
    for f in $(outputs_for "$bin"); do
      if [ -f "$cd/$f" ] || [ -f "$rd/$f" ]; then
        cmp -s "$cd/$f" "$rd/$f" || { verdict="differs:$f"; break; }
      fi
    done
  fi
  rm -rf "$cd" "$rd"
  echo "$verdict"
}

best_ms() {
  local best=999999 s e ms i d
  for i in 1 2 3 4 5; do
    d=$(scratch)
    s=$(date +%s%N); ( cd "$d" && "$@" ) >/dev/null 2>&1; e=$(date +%s%N)
    rm -rf "$d"
    ms=$(( (e - s) / 1000000 ))
    [ "$ms" -lt "$best" ] && best=$ms
  done
  echo "$best"
}

printf 'tag\tbinary\tc_ir\tr_ir\tir_ratio\tc_rss_kb\tr_rss_kb\trss_ratio\tc_ms\tr_ms\tms_ratio\tnote\n'
printf '%s\n' "$CASES" | while IFS=$'\t' read -r tag bin input argv; do
  [ -n "$tag" ] || continue
  cbin="$CORE/$bin"
  rbin="$RSDIR/$bin-rs"
  cd "$work" || exit 1
  # shellcheck disable=SC2086
  set -- $argv -i "$input"
  verdict=$(check_pair "$bin" "$cbin" "$rbin" "$@")
  if [ "$verdict" != ok ]; then
    printf '%s\t%s\tSKIPPED\t%s\t-\t-\t-\t-\t-\t-\t-\n' "$tag" "$bin" "$verdict" 1>&2
    printf '%s\t%s\t-\t-\t-\t-\t-\t-\t-\t-\t-\t%s\n' "$tag" "$bin" "$verdict"
    continue
  fi
  cir=$(ir_count "$cbin" "$@")
  rir=$(ir_count "$rbin" "$@")
  crss=$(peak_rss "$cbin" "$@")
  rrss=$(peak_rss "$rbin" "$@")
  cms=$(best_ms "$cbin" "$@")
  rms=$(best_ms "$rbin" "$@")
  awk -v t="$tag" -v b="$bin" -v ci="$cir" -v ri="$rir" -v cr="$crss" -v rr="$rrss" \
      -v cm="$cms" -v rm="$rms" 'BEGIN{
        printf "%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%.3f\t%s\t%s\t%.3f\tok\n",
               t, b, ci, ri, (ci?ri/ci:0), cr, rr, (cr?rr/cr:0), cm, rm, (cm?rm/cm:0)
      }'
done
