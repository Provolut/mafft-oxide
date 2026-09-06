#!/bin/bash
# Run a differential harness without the C, against a recording of what the C
# did.
#
#   verify_golden.sh record <harness.sh> <c-binary> <rust-binary> <store.tar.gz>
#   verify_golden.sh replay <harness.sh> <store.tar.gz> <rust-binary>
#
#   ./verify_golden.sh record verify.sh ../core/disttbfast \
#       target/release/disttbfast-rs tests/golden/disttbfast.tar.gz
#   ./verify_golden.sh replay verify.sh tests/golden/disttbfast.tar.gz \
#       target/release/disttbfast-rs
#
# The harnesses take the C binary as an argument and otherwise care nothing
# about it, so `goldenc` -- which reproduces the recorded bytes -- takes its
# place. Neither harness is modified, neither knows, and both make exactly the
# comparisons they always made and report the same tallies. What goes away is
# the half of the run that recomputes an answer nobody has changed.
#
# Recording happens during a real differential run, C and port side by side, so
# a store can only be written by a run that passed. There is no mode that
# records without checking.
#
# WHAT THIS IS AND IS NOT
#
# Replaying is a comparison against a recording, and a recording is only as good
# as its provenance. Two things follow, and both are enforced rather than
# documented:
#
#  * The store records which C produced it -- the binary's size and content
#    hash. Replaying with a `core/` that no longer matches is refused, because
#    the thing that would otherwise happen is the expensive one: a green run
#    that checked the port against a previous upstream's answers.
#  * A case with no recording is a hard error, never a skip and never a fall
#    back to running something. A store that has lost cases must stop the run.
#
# So this is the fast gate, and `verify.sh` against a real `core/` remains the
# authority. Re-record after any change to the C, to `gencorpus.py`, or to a
# harness's argument vectors -- all three change what the cases are.

set -u
here=$(cd "$(dirname "$0")" && pwd)
MODE="${1:?usage: verify_golden.sh record|replay ...}"
HARNESS="${2:?usage: verify_golden.sh record|replay <harness.sh> ...}"
case "$HARNESS" in /*) ;; *) HARNESS="$here/$HARNESS" ;; esac
[ -x "$HARNESS" ] || { echo "verify_golden: no such harness: $HARNESS" 1>&2; exit 2; }

GOLDENC="$here/target/release/examples/goldenc"
if [ ! -x "$GOLDENC" ]; then
  echo "verify_golden: building goldenc" 1>&2
  ( cd "$here" && cargo build --release --example goldenc ) || exit 2
fi

abspath() { case "$1" in /*) printf '%s\n' "$1" ;; *) printf '%s/%s\n' "$(pwd)" "$1" ;; esac; }

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# One tally line from a full parallel run of the harness against $1 as the C.
tally() {
  local o
  o=$("$here/verify_parallel.sh" "$(basename "$HARNESS")" "$1" "$RBIN" 2>&1) || {
    printf '%s\n' "$o" 1>&2
    return 1
  }
  printf '%s\n' "$o" | grep -E '^pass=[0-9]+ fail=[0-9]+$' | tail -1
}

# `verify_tbfast.sh` runs `$(dirname "$CBIN")/disttbfast` to stage an alignment,
# so the stand-in needs its neighbours as well as itself -- see the invariant in
# the record path.
make_wrappers() {
  local mode="$1" realdir=""
  [ -n "${2:-}" ] && realdir=$(dirname "$2")
  local b
  for b in cbin disttbfast tbfast dvtditr; do
    local target=""
    if [ "$mode" = record ]; then
      case "$b" in
        cbin) target="$CBIN" ;;
        *) target="$realdir/$b"; [ -x "$target" ] || continue ;;
      esac
      cat > "$work/bin/$b" <<EOF
#!/bin/sh
exec "$GOLDENC" record "$store" "$target" "\$@"
EOF
    else
      cat > "$work/bin/$b" <<EOF
#!/bin/sh
exec "$GOLDENC" replay "$store" "\$@"
EOF
    fi
    chmod +x "$work/bin/$b"
  done
}

# The identity of the C a store was made from. Cheap, and enough to catch the
# case that matters: a rebuilt or reimported `core/`.
c_identity() {
  printf '%s %s\n' "$(wc -c < "$1" | tr -d ' ')" "$(cksum < "$1" | awk '{print $1}')"
}

case "$MODE" in

record)
  CBIN=$(abspath "${3:?record needs the C binary}")
  RBIN=$(abspath "${4:?record needs the Rust binary}")
  OUT=$(abspath "${5:?record needs an output archive}")
  [ -x "$CBIN" ] || { echo "verify_golden: $CBIN is not executable" 1>&2; exit 2; }

  store="$work/store"
  mkdir -p "$store" "$work/bin"
  make_wrappers record "$CBIN"

  echo "verify_golden: reference run against the real C"
  ref=$(tally "$CBIN") || { echo "verify_golden: the differential fails against the real C; nothing to record" 1>&2; exit 1; }

  echo "verify_golden: recording $(basename "$HARNESS") from $CBIN"
  # Through the parallel driver: recording is a full differential run, and there
  # is no reason to wait nine minutes for it. `goldenc` writes each case under a
  # content-derived name and renames it into place, so the workers cannot
  # collide destructively.
  got=$(tally "$work/bin/cbin") || { echo "verify_golden: the differential failed while recording" 1>&2; exit 1; }

  # The invariant that matters, and it is not theoretical: `verify_tbfast.sh`
  # looks for `disttbfast` beside the C binary, so the first version of this
  # script -- whose wrapper lived alone in a temp directory -- silently skipped
  # a 592-case section and recorded a run that was 1098 cases instead of 1690.
  # It reported a clean pass while checking two thirds of what it claimed. So
  # the recording run is required to have the same tally as a run against the
  # real C, and anything else stops here.
  if [ "$ref" != "$got" ]; then
    echo "verify_golden: the recording run did not check the same cases" 1>&2
    echo "verify_golden:   against the real C: $ref" 1>&2
    echo "verify_golden:   through goldenc:    $got" 1>&2
    echo "verify_golden: something the harness reaches for is missing beside the" 1>&2
    echo "verify_golden: stand-in. Refusing to record a partial store." 1>&2
    exit 1
  fi

  n=$(find "$store" -maxdepth 1 -mindepth 1 -type d | wc -l)
  [ "$n" -gt 0 ] || { echo "verify_golden: nothing was recorded" 1>&2; exit 2; }
  c_identity "$CBIN" > "$store/C-IDENTITY"
  printf '%s\n' "$ref" > "$store/TALLY"
  mkdir -p "$(dirname "$OUT")"
  ( cd "$store" && tar -cf - . ) | gzip -9 > "$OUT"
  echo "verify_golden: $n cases recorded, tally $ref -> $OUT ($(du -h "$OUT" | cut -f1))"
  ;;

replay)
  ARCHIVE=$(abspath "${3:?replay needs the store archive}")
  RBIN=$(abspath "${4:?replay needs the Rust binary}")
  [ -f "$ARCHIVE" ] || { echo "verify_golden: no store at $ARCHIVE" 1>&2; exit 2; }

  store="$work/store"
  mkdir -p "$store"
  gzip -dc "$ARCHIVE" | ( cd "$store" && tar -xf - ) || {
    echo "verify_golden: cannot read $ARCHIVE" 1>&2; exit 2; }

  # If a `core/` is present, it must be the one this store was recorded from.
  # Absent, replaying is exactly the point and nothing is checked.
  guess="$here/../core/$(basename "$HARNESS" | sed 's/^verify_//; s/^verify$/disttbfast/; s/\.sh$//')"
  if [ -f "$store/C-IDENTITY" ] && [ -x "$guess" ]; then
    have=$(c_identity "$guess")
    want=$(cat "$store/C-IDENTITY")
    if [ "$have" != "$want" ]; then
      echo "verify_golden: $guess does not match the C this store was recorded from" 1>&2
      echo "verify_golden:   store: $want" 1>&2
      echo "verify_golden:   local: $have" 1>&2
      echo "verify_golden: re-record, or move that binary aside to replay anyway." 1>&2
      exit 2
    fi
  fi

  mkdir -p "$work/bin"
  make_wrappers replay ""

  out=$("$here/verify_parallel.sh" "$(basename "$HARNESS")" "$work/bin/cbin" "$RBIN" 2>&1)
  rc=$?
  printf '%s\n' "$out"
  got=$(printf '%s\n' "$out" | grep -E '^pass=[0-9]+ fail=[0-9]+$' | tail -1)

  # A store that has lost cases, or a harness whose case list has moved on,
  # would otherwise show up as a smaller green run -- the failure mode this
  # whole approach is most exposed to.
  if [ -f "$store/TALLY" ]; then
    want=$(cat "$store/TALLY")
    wantn=${want#pass=}; wantn=${wantn%% *}
    gotn=${got#pass=}; gotn=${gotn%% *}
    gotf=${got##*fail=}
    if [ "$rc" -eq 0 ] && [ "$(( gotn + gotf ))" != "$(( wantn ))" ]; then
      echo "verify_golden: this run checked $(( gotn + gotf )) cases; the store was recorded from $wantn" 1>&2
      echo "verify_golden: the store and the harness no longer agree on what the cases are." 1>&2
      exit 1
    fi
  fi
  exit $rc
  ;;

*)
  echo "verify_golden: unknown mode '$MODE'" 1>&2
  exit 2
  ;;
esac
