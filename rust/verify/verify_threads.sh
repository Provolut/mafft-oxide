#!/bin/bash
# Differential test: the port's output must not depend on how many threads ran.
#
#   verify_threads.sh <path-to-C-disttbfast> <path-to-disttbfast-rs> [corpus-dir]
#
# `disttbfast-rs` threads the three all-pairs matrices and nothing else. Each
# cell of those is a pure function of its own (i, j) and of data the loop does
# not write, so the result is thread-count independent by construction. This
# harness is here because "by construction" is an argument, and this repository
# does not ship arguments: it ships byte comparisons.
#
# Two questions, not one, and the second is the one that is easy to forget:
#
#   1. Does the port give the same bytes at `-C 8-8` as at `-C 0-0`?
#   2. Are those bytes still the C's?
#
# A port that was self-consistently wrong at every thread count would pass (1)
# and is exactly what (2) catches, so every thread count is compared against
# the *C's* sequential run rather than against the port's own.
#
# The thread counts are not all powers of two on purpose. Rows are dealt
# round-robin across workers, so a count that does not divide the row count
# leaves one worker with an extra row, and a count larger than the row count
# leaves some workers with none; 3, 5 and 64 reach those. `-C 8-0` and `-C 0-8`
# check that the two halves of `-C a-b` are read the right way round -- a
# transposed parse would still pass every symmetric spec.
set -u
CBIN="${1:?usage: verify_threads.sh <c-disttbfast> <rust-disttbfast-rs> [corpus-dir]}"
RBIN="${2:?usage: verify_threads.sh <c-disttbfast> <rust-disttbfast-rs> [corpus-dir]}"
CORPUS="${3:-}"

abspath() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s/%s\n' "$(pwd)" "$1" ;; esac; }
CBIN=$(abspath "$CBIN")
RBIN=$(abspath "$RBIN")
for b in "$CBIN" "$RBIN"; do
  [ -x "$b" ] || { echo "verify_threads: $b is not executable" 1>&2; exit 2; }
done

here=$(cd "$(dirname "$0")" && pwd)
rust=$(cd "$here/.." && pwd)
repo=$(cd "$rust/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

if [ -z "$CORPUS" ]; then
  CORPUS="$work/corpus"
  mkdir -p "$CORPUS"
  python3 "$rust/gen/gencorpus.py" "$repo/test/sample" "$CORPUS" full >/dev/null
fi
n=$(ls -1 "$CORPUS"/*.fa 2>/dev/null | wc -l)
[ "$n" -gt 0 ] || { echo "verify_threads: no inputs in $CORPUS" 1>&2; exit 2; }

# The strategies `--auto` can reach with threads. `-E 2 -F` is the production
# shape (ktuples, iterate=0, cycle=2, which is what any input from 500 to 20000
# sequences selects); the others vary the pieces that decide which matrices run
# at all, since `-E 1` skips the recompute over the aligned rows entirely.
BASE="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -b 62 -g 0 -f -1.53 -Q 100.0 -h 0"
TAIL="-X 0.1 -x 100"
SHAPES="fftns2:$BASE -F|nwns2:$BASE|memsave:$BASE -F -M -B|cycle1:-q 0 -E 1 -V -1.53 -s 0.0 -W 6 -O -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F"
SPECS="1-1 2-2 3-3 5-5 8-8 64-64 8-0 0-8"

pass=0
fail=0
skip=0
ran_threaded=0

run() { # run <dir> <bin> <argv...>; prints exit status
  local d=$1; shift
  local b=$1; shift
  mkdir -p "$d"
  ( cd "$d" && "$b" "$@" < "$INPUT" > stdout.fa 2> stderr.txt )
  echo $?
}

for f in "$CORPUS"/*.fa; do
  INPUT=$f
  base=$(basename "$f" .fa)
  echo "$SHAPES" | tr '|' '\n' | while IFS= read -r shape; do
    label=${shape%%:*}
    argv=${shape#*:}
    d="$work/$base.$label"
    rm -rf "$d"

    # The C's sequential run is the reference for every spec below.
    # shellcheck disable=SC2086
    cst=$(run "$d/ref" "$CBIN" $argv $TAIL -C 0-0)
    if [ "$cst" -ne 0 ]; then
      echo "SKIP $base [$label]: the C exits $cst"
      continue
    fi
    # A reference that produced nothing would make every comparison below
    # vacuous, which is how this repository's harnesses have been fooled twice.
    if [ ! -s "$d/ref/stdout.fa" ]; then
      echo "VACUOUS $base [$label]: the C wrote no alignment"
      continue
    fi

    for spec in $SPECS; do
      # shellcheck disable=SC2086
      rst=$(run "$d/t$spec" "$RBIN" $argv $TAIL -C "$spec")
      if [ "$rst" -ne 0 ]; then
        echo "FAIL $base [$label] -C $spec: the port exits $rst"
        sed 's/^/    /' "$d/t$spec/stderr.txt" | head -5
        echo "__FAIL__"
        continue
      fi
      if ! cmp -s "$d/t$spec/stdout.fa" "$d/ref/stdout.fa"; then
        echo "FAIL $base [$label] -C $spec: differs from the C's -C 0-0"
        echo "__FAIL__"
        continue
      fi
      if [ -f "$d/ref/order" ] && ! cmp -s "$d/t$spec/order" "$d/ref/order"; then
        echo "FAIL $base [$label] -C $spec: \`order\` differs from the C's"
        echo "__FAIL__"
        continue
      fi
      echo "__PASS__"
    done
    rm -rf "$d"
  done
done > "$work/log"

pass=$(grep -c '^__PASS__$' "$work/log" || true)
fail=$(grep -c '^__FAIL__$' "$work/log" || true)
skip=$(grep -c -E '^(SKIP|VACUOUS) ' "$work/log" || true)
grep -vE '^__(PASS|FAIL)__$' "$work/log" | grep -E '^(FAIL|SKIP|VACUOUS)' | head -40

echo "threads: pass=$pass fail=$fail skipped=$skip"

# Guard against the two ways this could report success having tested nothing:
# an empty corpus (checked above) and a loop that never reached a threaded run.
if [ "$((pass + fail))" -eq 0 ]; then
  echo "verify_threads: compared nothing -- vacuous, not green" 1>&2
  exit 2
fi
nspecs=$(printf '%s\n' $SPECS | wc -l)
if [ "$pass" -lt "$nspecs" ]; then
  echo "verify_threads: only $pass passes for $nspecs thread specs; at least one" \
       "spec never ran to completion anywhere" 1>&2
  exit 2
fi
[ "$fail" -eq 0 ]
