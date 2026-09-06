#!/bin/bash
# Run one of the corpus differentials across every core, without changing what
# it checks.
#
#   verify_parallel.sh [-j N] <harness.sh> <c-binary> <rust-binary>
#
#   verify_parallel.sh verify.sh        ../core/disttbfast target/release/disttbfast-rs
#   verify_parallel.sh verify_tbfast.sh ../core/tbfast     target/release/tbfast-rs
#
# `verify.sh` takes 8m36s and `verify_tbfast.sh` 3m41s, and the time is real
# work: a case costs about 43ms, of which 42ms is the two binaries actually
# aligning something and well under a millisecond is the shell. So there is
# nothing to shave -- the only thing wrong with these harnesses is that they use
# one of twelve cores.
#
# Every case in both harnesses sits inside a `for f in "$CORPUS"/*` loop, and no
# case reads more than one corpus file. So the corpus can be dealt out among
# several copies of the harness and the results added up: same inputs, same
# argument vectors, same comparisons, same count. Nothing is sampled, skipped or
# weakened -- this is the identical run, dealt out.
#
# The corpus comes from the harness itself, via VERIFY_EMIT_CORPUS, rather than
# being generated a second time here. A driver that grew its own idea of what
# the corpus is would drift from the harness the moment either changed, and the
# drift would look like a passing run.
#
# Slices are at least 10 files because that is the floor each harness refuses to
# run below -- a guard against a vacuous run, which is exactly what slicing
# could otherwise manufacture. That caps the workers at corpus/10 (8 on the 88
# file corpus), which is the honest limit rather than a tuning choice.

set -u
here=$(cd "$(dirname "$0")" && pwd)

JOBS=$(nproc 2>/dev/null || echo 4)
if [ "${1:-}" = "-j" ]; then
  JOBS="${2:?usage: verify_parallel.sh [-j N] <harness.sh> <c-binary> <rust-binary>}"
  shift 2
fi
HARNESS="${1:?usage: verify_parallel.sh [-j N] <harness.sh> <c-binary> <rust-binary>}"
CBIN="${2:?usage: verify_parallel.sh [-j N] <harness.sh> <c-binary> <rust-binary>}"
RBIN="${3:?usage: verify_parallel.sh [-j N] <harness.sh> <c-binary> <rust-binary>}"

case "$HARNESS" in /*) ;; *) HARNESS="$here/$HARNESS" ;; esac
[ -x "$HARNESS" ] || { echo "verify_parallel: no such harness: $HARNESS" 1>&2; exit 2; }
# The harness resolves these itself, but it does so relative to its own cwd, and
# the workers run from here. Resolve once, up front.
case "$CBIN" in /*) ;; *) CBIN="$(pwd)/$CBIN" ;; esac
case "$RBIN" in /*) ;; *) RBIN="$(pwd)/$RBIN" ;; esac
for b in "$CBIN" "$RBIN"; do
  [ -x "$b" ] || { echo "verify_parallel: $b is not executable" 1>&2; exit 2; }
done

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

VERIFY_EMIT_CORPUS="$work/corpus" "$HARNESS" "$CBIN" "$RBIN" > "$work/emit.log" 2>&1 || {
  echo "verify_parallel: the harness would not produce its corpus:" 1>&2
  cat "$work/emit.log" 1>&2
  exit 2
}
n=$(find "$work/corpus" -maxdepth 1 -type f | wc -l)
[ "$n" -ge 10 ] || { echo "verify_parallel: corpus has $n files; refusing" 1>&2; exit 2; }

# At least 10 files per slice, so no worker trips its own vacuity guard.
max_jobs=$(( n / 10 ))
[ "$max_jobs" -lt 1 ] && max_jobs=1
[ "$JOBS" -gt "$max_jobs" ] && JOBS=$max_jobs
[ "$JOBS" -lt 1 ] && JOBS=1

# Dealt round-robin rather than in blocks: `gencorpus.py` emits related inputs
# next to each other, and a block split would hand one worker all the big ones.
i=0
for f in "$work/corpus"/*; do
  [ -f "$f" ] || continue
  d="$work/slice.$(( i % JOBS ))"
  mkdir -p "$d"
  cp "$f" "$d/"
  i=$(( i + 1 ))
done

# Each harness has one block of cases that build their own input instead of
# reading the corpus (`gen_text.py`'s). Those would otherwise run once per
# worker, multiplying the pass count by the worker count while checking nothing
# new -- so worker 0 runs them and the rest skip them. The check that this is
# right is arithmetic: the parallel tally must equal the serial one exactly, and
# it does (12018 and 1690).
echo "verify_parallel: $(basename "$HARNESS") over $n corpus files, $JOBS workers"
for j in $(seq 0 $(( JOBS - 1 ))); do
  skip=1
  [ "$j" = 0 ] && skip=""
  ( VERIFY_SKIP_STANDALONE="$skip" "$HARNESS" "$CBIN" "$RBIN" "$work/slice.$j" \
      > "$work/out.$j" 2>&1; echo $? > "$work/rc.$j" ) &
done
wait

pass=0
fail=0
failed=""
broke=0
for j in $(seq 0 $(( JOBS - 1 ))); do
  rc=$(cat "$work/rc.$j" 2>/dev/null || echo 99)
  line=$(grep -E '^pass=[0-9]+ fail=[0-9]+$' "$work/out.$j" | tail -1)
  if [ -z "$line" ]; then
    # No tally at all means the worker died before finishing -- which must never
    # read as "nothing failed".
    broke=1
    failed="$failed\n  worker $j produced no result (exit $rc):\n$(sed 's/^/    /' "$work/out.$j" | tail -5)"
    continue
  fi
  p=${line#pass=}; p=${p%% *}
  f=${line##*fail=}
  pass=$(( pass + p ))
  fail=$(( fail + f ))
  if [ "$f" -ne 0 ]; then
    failed="$failed\n$(sed -n '/^FAILED:/,$p' "$work/out.$j" | tail -n +2)"
  fi
done

echo "pass=$pass fail=$fail"
if [ "$fail" -ne 0 ] || [ "$broke" -ne 0 ]; then
  printf 'FAILED:%b\n' "$failed"
  exit 1
fi
