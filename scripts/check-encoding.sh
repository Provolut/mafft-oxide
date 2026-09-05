#!/bin/bash
# Refuse a source file that is not valid UTF-8, and offer to convert it.
#
#     scripts/check-encoding.sh          # report; exit 1 if anything is not UTF-8
#     scripts/check-encoding.sh --fix    # convert from EUC-JP, then report again
#
# Upstream MAFFT ships eight files whose Japanese comments are EUC-JP
# (constants.c, io.c, mltaln9.c, Lalign11.c, Lalignmm.c, MSalignmm.c, plus one
# each under MPI/ and extensions/). That costs more than
# unreadable comments: GNU grep in a UTF-8 locale treats such a file as binary
# and reports *nothing* for it -- no match, no count, no warning -- so "this
# function is called nowhere" silently excludes a fifteenth of `core/`,
# including all of I/O, the memsave aligner, and mltaln9.c -- where `fillimp`
# lives, and where the comment "be careful not to apply the weight twice" sits
# a few lines from it. That has already produced two wrong conclusions in this
# repo's own documentation.
#
# One accepted gap: an untracked file is not checked. It becomes visible to
# `git ls-files` the moment it is staged, which is before any push can carry
# it, so the gate still sees it in time.
#
# Converting is safe and was checked rather than assumed: the affected bytes are
# all inside comments (no string literal contains one), and building `core/`
# before and after the conversion gives 34 of 34 binaries byte-identical.
#
# Run with --fix after importing a new upstream release, which will bring the
# EUC-JP files back.
set -u
here=$(cd "$(dirname "$0")/.." && pwd)
cd "$here" || exit 2
fix=0
[ "${1:-}" = "--fix" ] && fix=1

# Every tracked file, minus the built binaries this repo also tracks. The
# text/binary split is `grep -I` under LC_ALL=C, which decides on NUL bytes
# alone -- in a UTF-8 locale it would call an EUC-JP source binary and skip
# exactly the files this script exists to find. Deliberately not an extension
# whitelist: upstream imports arrive with .cpp, .hpp and .txt too, and a guard
# that only knows the extensions we have seen stops guarding the moment one
# arrives that we have not.
mapfile -t -d '' files < <(git ls-files -z)

bad=()
for f in "${files[@]}"; do
  [ -f "$f" ] || continue
  LC_ALL=C grep -qI . "$f" 2>/dev/null || continue   # binary; not our business
  iconv -f UTF-8 -t UTF-8 "$f" > /dev/null 2>&1 || bad+=("$f")
done

if [ ${#bad[@]} -eq 0 ]; then
  echo "check-encoding: ${#files[@]} tracked sources, all valid UTF-8"
  exit 0
fi

if [ $fix -eq 0 ]; then
  echo "check-encoding: not valid UTF-8:" 1>&2
  printf '  %s\n' "${bad[@]}" 1>&2
  echo "check-encoding: grep in a UTF-8 locale silently skips these files." 1>&2
  echo "check-encoding: run 'scripts/check-encoding.sh --fix' to convert them." 1>&2
  exit 1
fi

for f in "${bad[@]}"; do
  if iconv -f EUC-JP -t UTF-8 "$f" -o "$f.utf8tmp" 2>/dev/null; then
    mv "$f.utf8tmp" "$f"
    echo "check-encoding: converted $f from EUC-JP"
  else
    rm -f "$f.utf8tmp"
    echo "check-encoding: $f is neither UTF-8 nor EUC-JP; convert it by hand" 1>&2
    exit 1
  fi
done

# Re-check, so --fix cannot report success on a file it failed to convert.
exec "$0"
