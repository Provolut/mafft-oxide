#!/bin/bash
# Region and line coverage of `cargo test`, as `GOAL.md` defines it.
#
#   ./coverage.sh [--show FILE]
#
# Two things make this easy to get wrong, and both under-report badly:
#
#  * The golden tests drive `disttbfast-rs` and `tbfast-rs` as subprocesses, so
#    the binaries must be passed to `llvm-cov` as `--object` alongside the test
#    harnesses. Leaving them out reported 56% where the answer was 74%.
#  * Every `cargo` call in the sequence must carry the same `RUSTFLAGS`.
#    One without it silently rebuilds the binaries uninstrumented, and then
#    `llvm-cov` fails with "no coverage data found" -- or worse, succeeds
#    against a stale profile.
#
# So the flags are set once here and every cargo invocation inherits them.

set -eu
here=$(cd "$(dirname "$0")" && pwd)
cd "$here/.."

prof=$(mktemp -d)
trap 'rm -rf "$prof"' EXIT

export RUSTFLAGS="-C instrument-coverage"
export LLVM_PROFILE_FILE="$prof/%p-%m.profraw"

bindir=$(rustc --print sysroot)/lib/rustlib/$(rustc -vV | sed -n 's/host: //p')/bin
PROFDATA="$bindir/llvm-profdata"
LLVMCOV="$bindir/llvm-cov"
for t in "$PROFDATA" "$LLVMCOV"; do
  [ -x "$t" ] || { echo "coverage: $t missing (rustup component add llvm-tools-preview)" 1>&2; exit 2; }
done

cargo build --tests --bins --quiet
# Collect the object list before running, so no cargo call happens afterwards
# with different flags.
objects=$(cargo test --tests --no-run --message-format=json 2>/dev/null | python3 -c '
import json, sys
for line in sys.stdin:
    try:
        m = json.loads(line)
    except ValueError:
        continue
    if m.get("profile", {}).get("test") and m.get("executable"):
        print(m["executable"])
')
cargo test --tests --quiet > /dev/null

# Every binary the tests drive as a subprocess has to be listed. Adding one to
# the workspace without adding it here does not lower the total -- it reports
# the new binary at 0.00% while its tests are in fact exercising it, which is
# the same under-report the header describes, pointed the other way.
args=()
for o in $objects \
         target/debug/disttbfast-rs target/debug/tbfast-rs \
         target/debug/countlen-rs target/debug/version-rs; do
  args+=(--object "$o")
done

# A directory of 500+ profraw files overflows the argument list; -f takes a
# file of paths instead.
ls "$prof"/*.profraw > "$prof/list"
"$PROFDATA" merge -sparse -f "$prof/list" -o "$prof/all.profdata"

if [ "${1:-}" = "--show" ]; then
  # Uncovered lines of one file, for finding what to test next.
  "$LLVMCOV" show "${args[@]}" --instr-profile="$prof/all.profdata" "${2:?--show needs a file}" \
    | awk -F'|' '$2 ~ /^ *0 *$/ {print $1 $3}'
  exit 0
fi

"$LLVMCOV" report "${args[@]}" --instr-profile="$prof/all.profdata" \
  --ignore-filename-regex='(/rustc/|cargo/registry)' \
  | awk '{printf "%-26s regions %8s  lines %8s\n", $1, $4, $10}' \
  | grep -E '^src|^tests|^TOTAL'
