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

# conda-forge asks every Rust recipe to vendor its dependencies' licence text.
# This crate has none -- `rust/Cargo.toml` declares no `[dependencies]` at all
# -- so what this writes describes only the crate itself. Generated anyway:
# the requirement is that the artifact be self-describing, and an empty answer
# that was produced by the tool beats an absent one backed by an argument.
cargo-bundle-licenses --format yaml --output "$SRC_DIR/THIRDPARTY.yml"
popd

test -s "$SRC_DIR/THIRDPARTY.yml" || { echo "cargo-bundle-licenses wrote nothing" >&2; exit 1; }

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
# substitution over `_LIBDIR` -- plus one added fallback line.
#
# `_LIBDIR` becomes the absolute libexec path, which is what upstream and
# conda-build expect: conda-build finds the build prefix in a text file,
# records it in `has_prefix`, and rewrites it at install time. Verified to
# happen. Two reasons that literal matters rather than being merely
# traditional. It is how the conda-forge C `mafft` package resolves its own
# components, so this package resolves by the same mechanism as the package it
# substitutes for. And it is *statically readable*: a tool that wants to know
# which binaries a driver will run can read the line. A self-locating
# `$(cd ...)` is opaque to that -- anything containing `$` is a shell expansion
# an inspector should refuse to evaluate rather than guess at -- so a driver
# with no literal anywhere makes such a check fail closed. The consumer's
# build-parity test does exactly that check and did exactly that.
#
# The fallback covers what the literal cannot. A copied or relocated
# environment keeps the old absolute path, which does not fail -- it silently
# resolves to the components of the *original* installation, so the launcher in
# environment B runs environment A's binaries. Measured: a `cp -a` of the
# installed env produced a working launcher still pointing at the source env's
# libexec. So `[ -d "$prefix" ]` guards the literal and self-location takes
# over only when that directory is gone, which is precisely the moved-or-copied
# case -- and the case where there is nothing static to read anyway. Normal
# installs never reach it, so what an inspector reads is what runs.
#
# Order matters and the obvious order is wrong. Literal-as-fallback, with
# self-location first, would leave an inspector reading a literal the driver
# does not use -- reporting a component directory that is not the one in use.
# Confidently wrong beats refusing only in the sense that it is quieter.
#
# `$0` is the path the shell resolved, including via PATH, and this runs before
# the driver changes directory, so a relative invocation still resolves. Note
# the whole block sits in the `else` of `if [ "$MAFFT_BINARIES" ]`, so an
# explicit `MAFFT_BINARIES` still wins over both -- which is what lets this
# package be compared against another build.
#
# awk on an exact line match rather than sed: the replacement text is full of
# characters sed would reinterpret (`&`, `$`, quotes), and matching the whole
# line guarantees nothing else in a 3069-line script is touched.
FALLBACK='\t[ -d "$prefix" ] || prefix="$(cd "$(dirname "$0")/../libexec/mafft-oxide" && pwd)"'
awk -v lib="$LIBEXEC" -v fb="$FALLBACK" '
  $0 == "\tprefix=_LIBDIR" { print "\tprefix=" lib; print fb; n++; next }
  { print }
  END { if (n != 1) { print "expected exactly one prefix=_LIBDIR line, saw " n > "/dev/stderr"; exit 1 } }
' "$SRC_DIR/core/mafft.tmpl" > "$PREFIX/bin/mafft-oxide"
chmod 755 "$PREFIX/bin/mafft-oxide"

# That substitution is the one piece of this recipe that edits executable logic
# rather than copying it, so check both lines landed and that nothing else
# changed.
grep -qxF "	prefix=$LIBEXEC" "$PREFIX/bin/mafft-oxide" \
  || { echo "the _LIBDIR substitution did not bake the literal libexec path" >&2; exit 1; }
grep -qF '[ -d "$prefix" ] || prefix="$(cd "$(dirname "$0")/../libexec/mafft-oxide" && pwd)"' \
  "$PREFIX/bin/mafft-oxide" \
  || { echo "the self-locating fallback line is missing" >&2; exit 1; }
test "$(grep -cF '[ -d "$prefix" ] || prefix=' "$PREFIX/bin/mafft-oxide")" = 1 \
  || { echo "the fallback line appears more than once" >&2; exit 1; }
# Drop the one line this recipe adds, normalise the two `prefix=` lines that
# upstream already has, and the rest must be byte-identical to the template.
d=$(diff <(sed 's:^\tprefix=.*:PREFIXLINE:' "$SRC_DIR/core/mafft.tmpl") \
         <(grep -vF '[ -d "$prefix" ] || prefix=' "$PREFIX/bin/mafft-oxide" \
           | sed 's:^\tprefix=.*:PREFIXLINE:') | wc -l)
test "$d" = 0 \
  || { echo "generated driver differs from mafft.tmpl beyond the prefix lines" >&2; exit 1; }

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
