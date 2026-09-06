#!/bin/bash
# Differential test for the `Gaplen` half of `-@` (D__align): run
# `gaplencount` and the three live `calcpfac` functions in both the C and the
# port on the same groups and require bit-identical results.
#
#   verify_gaplen.sh <path-to-core-objdir>
#
# The object directory is a built `core/` whose Dalignmm.c has had `static`
# removed from those four functions -- they are file-local in the shipped
# source. Everything is compared as raw f64 bit patterns, so a formatting
# difference cannot hide a value difference.
set -u
OBJDIR="${1:?usage: verify_gaplen.sh <core-objdir>}"
here=$(cd "$(dirname "$0")" && pwd)
rust=$(cd "$here/.." && pwd)
repo=$(cd "$rust/.." && pwd)
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT

objs=$(ls "$OBJDIR"/*.o | grep -v 'disttbfast\.o$')
if [ -z "$objs" ]; then echo "verify_gaplen: no objects in $OBJDIR" 1>&2; exit 2; fi
gcc -O0 -g -o "$work/glharness" "$rust/gen/gen_glharness.c" $objs -lm -lpthread || exit 2
( cd "$rust" && cargo build --release --example glref ) || exit 2
RREF="$rust/target/release/examples/glref"

python3 "$rust/gen/gen_glcases.py" > "$work/cases.txt"
n=$(wc -l < "$work/cases.txt")
if [ "$n" -lt 50 ]; then echo "verify_gaplen: only $n cases; refusing" 1>&2; exit 2; fi

pass=0; fail=0; failed=""
while IFS= read -r line; do
  IFS=$'\t' read -r -a args <<< "$line"
  "$work/glharness" "${args[@]}" > "$work/c.txt" 2>/dev/null
  "$RREF" "${args[@]}" > "$work/r.txt" 2>/dev/null || { fail=$((fail+1)); failed="$failed\n  ${args[*]:0:6} (port failed)"; continue; }
  if cmp -s "$work/c.txt" "$work/r.txt"; then pass=$((pass+1)); else
    fail=$((fail+1)); failed="$failed\n  ${args[*]:0:6}"
  fi
done < "$work/cases.txt"

echo "gaplen: pass=$pass fail=$fail"
if [ $fail -ne 0 ]; then printf 'FAILED:%b\n' "$failed" | head -20; exit 1; fi
