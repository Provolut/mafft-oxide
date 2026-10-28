#!/bin/bash
# Build mafft-oxide and install it as a self-contained MAFFT.
#
# Three things are installed, and the shape of them is the whole point:
#
#  1. The nine ported components into a private libexec, **under the names the
#     C driver calls** -- `tbfast`, not `tbfast-rs`. `core/mafft.tmpl` contains
#     zero occurrences of "-rs", so renaming on install is what lets the
#     upstream driver run the port with no patch to the driver at all.
#  2. The driver itself, generated the way `core/Makefile:168-169` generates it:
#     one `sed` replacing `_LIBDIR` with the directory from (1). That is the
#     entire difference between `mafft.tmpl` and a working `mafft`.
#  3. Nothing else. In particular no `etc/conda/activate.d` script and no
#     `MAFFT_BINARIES` export: the driver's `_LIBDIR` already points at (1),
#     so the launcher is self-contained. `MAFFT_BINARIES` is how MAFFT's driver
#     locates its components, so exporting it environment-wide would silently
#     redirect *every* `mafft` invocation in the environment, including the
#     "stock" arm of anyone's differential test -- which then compares this
#     port against itself and passes unconditionally. That failure has already
#     happened on a developer machine here. The launcher is installed under its
#     own name instead, so a stock `mafft` in the same environment stays stock.
set -euxo pipefail

LIBEXEC="$PREFIX/libexec/mafft-oxide"
mkdir -p "$LIBEXEC" "$PREFIX/bin"

# The crate has no dependencies, so `--locked --offline` is honest here rather
# than aspirational: there is nothing to resolve and nothing to fetch.
pushd "$SRC_DIR/rust"
cargo build --release --locked --offline
popd

# conda's rust activation exports CARGO_BUILD_TARGET, which puts the output
# under `target/<triple>/release/` rather than `target/release/`. Handle both,
# so the recipe also builds outside conda-build (where it is unset).
RELDIR="$SRC_DIR/rust/target/${CARGO_BUILD_TARGET:+$CARGO_BUILD_TARGET/}release"
test -d "$RELDIR" || { echo "no release dir at $RELDIR" >&2; exit 1; }

# Rust binary -> the name the driver expects.
for c in disttbfast tbfast dvtditr dndpre countlen f2cl replaceu restoreu version; do
  install -m 755 "$RELDIR/${c}-rs" "$LIBEXEC/$c"
done

# Nine, and only nine. A tenth name here would be a component the driver calls
# that is not in fact ported.
test "$(ls -1 "$LIBEXEC" | wc -l)" = 9

# The driver. `conda-build` notices the build prefix in a text file and records
# it in `has_prefix`, so the path is rewritten to the real prefix at install
# time; the launcher therefore works from whatever directory the environment
# ends up in, with no environment variable.
sed "s:_LIBDIR:$LIBEXEC:" "$SRC_DIR/core/mafft.tmpl" > "$PREFIX/bin/mafft-oxide"
chmod 755 "$PREFIX/bin/mafft-oxide"

# The driver compares its own version string against what `$prefix/version`
# prints and warns on stdout if they disagree -- which would land in the
# alignment. Fail the build instead of shipping that.
want=$(sed -n 's/^version="v\([0-9.]*\).*/\1/p' "$SRC_DIR/core/mafft.tmpl")
got=$("$LIBEXEC/version")
test -n "$want"
test "$want" = "$got" || { echo "driver wants v$want, version-rs prints $got" >&2; exit 1; }
