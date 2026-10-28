#!/bin/bash
# Smoke test run by conda-build inside the installed environment.
#
# This is deliberately *not* the acceptance test. Byte-identity against the C
# is what makes this package correct, and that needs a compiled `core/` to
# compare against, which a package test does not have. The gate for that is
# `rust/verify/verify_pkg.sh`. What this checks is that the package is wired up
# at all: the nine are present as real executables under their C names, the
# driver resolves them with no environment variable set, and it produces an
# alignment.
set -euo pipefail

LIBEXEC="$PREFIX/libexec/mafft-oxide"
NINE="disttbfast tbfast dvtditr dndpre countlen f2cl replaceu restoreu version"

for c in $NINE; do
  test -x "$LIBEXEC/$c" || { echo "missing component: $LIBEXEC/$c" >&2; exit 1; }
  # Must be a real ELF/Mach-O executable. A wrapper script here would mean the
  # package is dispatching somewhere else entirely.
  head -c 4 "$LIBEXEC/$c" | od -An -tx1 | tr -d ' \n' | grep -qE '^(7f454c46|cffaedfe|cafebabe)' \
    || { echo "$LIBEXEC/$c is not a native executable" >&2; exit 1; }
done
echo "ok: 9 components present and native"

# The driver must find them with nothing exported. Unset the variable rather
# than trusting it to be absent, so the test cannot pass by inheriting a
# pointer to some other directory.
unset MAFFT_BINARIES

printf '>a\nMKVLATSW\n>b\nMKILATSW\n>c\nMKVLGTSW\n' > in.fa

mafft-oxide --retree 2 --maxiterate 0 in.fa > out.fa 2> err.txt
test -s out.fa || { echo "driver produced no output" >&2; cat err.txt >&2; exit 1; }
test "$(grep -c '^>' out.fa)" = 3 || { echo "expected 3 records" >&2; cat out.fa >&2; exit 1; }
echo "ok: driver aligns with no MAFFT_BINARIES set"

# The version warning the driver emits on a mismatch goes to stderr, but a
# mismatch means the libexec and the driver disagree about what they are.
grep -q '!=' err.txt && { echo "driver/component version mismatch" >&2; cat err.txt >&2; exit 1; }
echo "ok: driver and components agree on version"
