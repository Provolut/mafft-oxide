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

# The driver, generated the way the upstream Makefile generates it -- one
# substitution over `_LIBDIR` -- but substituting a path the script works out
# from its own location rather than a baked absolute one.
#
# An absolute `_LIBDIR` is what upstream and conda-build expect to work:
# conda-build finds the build prefix in a text file, records it in `has_prefix`,
# and rewrites it at install time, which was verified to happen. The reason not
# to rely on it is what happens when a prefix moves *after* install. A copied or
# relocated environment keeps the old absolute path, which does not fail -- it
# silently resolves to the components of the *original* installation, so the
# launcher in environment B runs environment A's binaries. Measured: a `cp -a`
# of the installed env produced a working launcher still pointing at the source
# env's libexec. Silently running a different installation is a worse failure
# than not starting, and this port exists not to produce quietly-different
# answers.
#
# Self-locating costs one `cd`. `$0` is the path the shell resolved, including
# via PATH, and this line runs before the driver changes directory, so a
# relative invocation still resolves. If the directory is absent the
# substitution yields an empty prefix and the driver fails at its first
# component rather than guessing.
#
# awk on an exact line match rather than sed: the replacement text is full of
# characters sed would reinterpret (`&`, `$`, quotes), and matching the whole
# line guarantees nothing else in a 3069-line script is touched.
REPL='"$(cd "$(dirname "$0")/../libexec/mafft-oxide" && pwd)"'
awk -v repl="$REPL" '
  $0 == "\tprefix=_LIBDIR" { print "\tprefix=" repl; n++; next }
  { print }
  END { if (n != 1) { print "expected exactly one prefix=_LIBDIR line, saw " n > "/dev/stderr"; exit 1 } }
' "$SRC_DIR/core/mafft.tmpl" > "$PREFIX/bin/mafft-oxide"
chmod 755 "$PREFIX/bin/mafft-oxide"

# That substitution is the one piece of this recipe that edits executable logic
# rather than copying it, so check it landed and that nothing else changed.
grep -qF 'prefix="$(cd "$(dirname "$0")/../libexec/mafft-oxide" && pwd)"' \
  "$PREFIX/bin/mafft-oxide" \
  || { echo "the _LIBDIR substitution did not produce the expected line" >&2; exit 1; }
d=$(diff <(sed 's:^\tprefix=.*:PREFIXLINE:' "$SRC_DIR/core/mafft.tmpl") \
         <(sed 's:^\tprefix=.*:PREFIXLINE:' "$PREFIX/bin/mafft-oxide") | wc -l)
test "$d" = 0 \
  || { echo "generated driver differs from mafft.tmpl beyond the prefix line" >&2; exit 1; }

# It must be a valid shell script. A substitution that produced something
# syntactically broken would otherwise only show up at run time.
bash -n "$PREFIX/bin/mafft-oxide"

# The driver compares its own version string against what `$prefix/version`
# prints and warns on stdout if they disagree -- which would land in the
# alignment. Fail the build instead of shipping that.
want=$(sed -n 's/^version="v\([0-9.]*\).*/\1/p' "$SRC_DIR/core/mafft.tmpl")
got=$("$LIBEXEC/version")
test -n "$want"
test "$want" = "$got" || { echo "driver wants v$want, version-rs prints $got" >&2; exit 1; }
