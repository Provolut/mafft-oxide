#!/bin/bash
# Differential test for the small utilities: run the C `countlen` / `version`
# and their ports on the same inputs and require byte-identical stdout, stderr
# and exit status.
#
#   verify_utils.sh <path-to-core-dir> [corpus-dir]
#   verify_utils.sh --regenerate <path-to-core-dir>
#
# With no corpus directory, the committed fixtures in tests/golden/utils/in are
# used, plus a set of generated edge cases. Point it at a directory of FASTA
# files to widen it -- that is the whole reason this exists next to
# tests/utils.rs, which pins the same behaviour without needing a C at all.
#
# `--regenerate` rewrites tests/golden/utils/*.expected from the C. Use it when
# the C changes or a case is added, never to accommodate a diff in the port:
# a port whose output moved is the finding, not the fixture.

set -u
here=$(cd "$(dirname "$0")" && pwd)

if [ "${1:-}" = "--regenerate" ]; then
  REGEN=1; shift
else
  REGEN=0
fi

CORE="${1:?usage: verify_utils.sh [--regenerate] <core-dir> [corpus-dir]}"
CORPUS="${2:-}"
RSDIR="$here/target/release"
GOLDEN="$here/tests/golden/utils"

# Absolute, because the `dndpre` cases below run each side in a scratch
# directory of its own -- `dndpre` writes `hat2` into the working directory, so
# comparing the two means changing into them, and a relative `../core` stops
# resolving the moment that happens.
CORE=$(cd "$CORE" 2>/dev/null && pwd) || { echo "verify_utils: no such core dir" 1>&2; exit 2; }

for b in "$CORE/countlen" "$CORE/version" "$CORE/replaceu" "$CORE/f2cl" "$CORE/dndpre" \
         "$RSDIR/countlen-rs" "$RSDIR/version-rs" "$RSDIR/replaceu-rs" "$RSDIR/f2cl-rs" \
         "$RSDIR/dndpre-rs"; do
  [ -x "$b" ] || { echo "verify_utils: $b is not executable" 1>&2; exit 2; }
done

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# Edge cases the committed corpus cannot carry as files: an empty input, a
# record with no trailing newline, and an input past the 100000-site sampling
# cap, which is where the C's N tally starts carrying the last sampled record
# forward (see issue #8 -- the port reproduces it deliberately).
gen_edge_cases() {
  local d="$1"
  printf '' > "$d/edge_empty.fa"
  printf '>a\nMKVLTA' > "$d/edge_nonewline.fa"
  printf '>a\n\n>b\nMKV\n' > "$d/edge_blankrecord.fa"
  printf '>a\n----\n>b\nMKVL\n' > "$d/edge_allgap.fa"
  printf '>a\nACGTNNNNACGT\n>b\nACGTACGTACGT\n' > "$d/edge_dna_n.fa"
  awk 'BEGIN{for(i=0;i<400;i++){printf(">s%d\n",i); for(j=0;j<400;j++) printf("%s","A"); printf("\n")}}' \
    > "$d/edge_over_cap.fa"
}

inputs="$work/inputs"
mkdir -p "$inputs"
cp "$GOLDEN"/in/* "$inputs/" 2>/dev/null
gen_edge_cases "$inputs"
if [ -n "$CORPUS" ]; then
  n=0
  for f in "$CORPUS"/*; do
    [ -f "$f" ] || continue
    cp "$f" "$inputs/corpus_$(basename "$f")"
    n=$((n+1))
  done
  echo "verify_utils: added $n file(s) from $CORPUS"
fi

if [ "$REGEN" = "1" ]; then
  : > "$GOLDEN/countlen.expected"
  for f in $(ls "$GOLDEN/in" | sort); do
    printf '%s\t%s\n' "$f" "$("$CORE/countlen" -i "$GOLDEN/in/$f" </dev/null 2>/dev/null)" \
      >> "$GOLDEN/countlen.expected"
  done
  "$CORE/version" > "$GOLDEN/version.expected"
  echo "verify_utils: regenerated $GOLDEN/*.expected from $CORE"
  exit 0
fi

pass=0; fail=0; indet=0

# Inputs whose *reference* output is undefined, so there is nothing to be
# bit-exact against. `countATGCandN` returns early for an empty sequence body
# and sets `*total` without setting `*countN`, so `getnumlen_nogap_countn`'s
# uninitialised `nN` reaches the `nfreq` division -- valgrind names the stack
# allocation, and the value it finds is not a contract (issue #8). The port
# initialises to zero and so prints `-nan` where this C happens to print `inf`.
#
# Listed here rather than dropped from the corpus: the run still executes them
# and still reports what each side said, it just does not score the comparison.
# Anything not on this list is a real failure.
INDETERMINATE="edge_blankrecord.fa"

is_indeterminate() {
  case " $INDETERMINATE " in *" $1 "*) return 0 ;; *) return 1 ;; esac
}

# Compare stdout and exit status only. For `f2cl`, whose stderr carries heap
# addresses from a `reporterr` left in the shipped binary -- different on every
# run of the C against itself, ASLR -- stderr cannot be a contract, so it is
# not compared. `rescale = 1`, printed from `constants()`, is on the same
# stream and is skipped for the same reason: the FASTA path never needs the
# scoring constants, and the script sends this stream to /dev/null.
cmp_stdout() {
  local label="$1" cbin="$2" rbin="$3"; shift 3
  local co cs ro rs
  co=$("$cbin" "$@" </dev/null 2>/dev/null); cs=$?
  ro=$("$rbin" "$@" </dev/null 2>/dev/null); rs=$?
  if [ "$co" = "$ro" ] && [ "$cs" = "$rs" ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL $label (status C=$cs RS=$rs)"
    [ "$co" = "$ro" ] || { echo "  stdout C : $(echo "$co" | head -c 200)"; echo "  stdout RS: $(echo "$ro" | head -c 200)"; }
  fi
}

# Compare one argument vector on both binaries: stdout, stderr and status.
cmp_run() {
  local label="$1" cbin="$2" rbin="$3"; shift 3
  local co ce cs ro re rs
  co=$("$cbin" "$@" </dev/null 2>"$work/ce"); cs=$?
  ro=$("$rbin" "$@" </dev/null 2>"$work/re"); rs=$?
  ce=$(cat "$work/ce"); re=$(cat "$work/re")
  if [ "$co" = "$ro" ] && [ "$ce" = "$re" ] && [ "$cs" = "$rs" ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL $label"
    [ "$co" = "$ro" ] || { echo "  stdout C : $co"; echo "  stdout RS: $ro"; }
    [ "$ce" = "$re" ] || { echo "  stderr C : $ce"; echo "  stderr RS: $re"; }
    [ "$cs" = "$rs" ] || echo "  status C: $cs  RS: $rs"
  fi
}

for f in "$inputs"/*; do
  base=$(basename "$f")
  if is_indeterminate "$base"; then
    indet=$((indet+1))
    echo "INDETERMINATE $base -- C: $("$CORE/countlen" -i "$f" </dev/null 2>/dev/null) | RS: $("$RSDIR/countlen-rs" -i "$f" </dev/null 2>/dev/null)"
    continue
  fi
  cmp_run "countlen -i $base" "$CORE/countlen" "$RSDIR/countlen-rs" -i "$f"
done

# The stdin path is a separate code path in both, and the script uses it twice.
for f in "$inputs"/*; do
  is_indeterminate "$(basename "$f")" && continue
  co=$("$CORE/countlen" < "$f" 2>/dev/null); cs=$?
  ro=$("$RSDIR/countlen-rs" < "$f" 2>/dev/null); rs=$?
  if [ "$co" = "$ro" ] && [ "$cs" = "$rs" ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL countlen < $(basename "$f")"
    echo "  C : $co (status $cs)"
    echo "  RS: $ro (status $rs)"
  fi
done

# Argument handling, where the C's loop is easy to get subtly wrong.
cmp_run "countlen -Z"            "$CORE/countlen" "$RSDIR/countlen-rs" -Z
cmp_run "countlen -Zi file"      "$CORE/countlen" "$RSDIR/countlen-rs" -Zi "$inputs/edge_dna_n.fa"
cmp_run "countlen positional"    "$CORE/countlen" "$RSDIR/countlen-rs" "$inputs/edge_dna_n.fa"
cmp_run "countlen -i (no arg)"   "$CORE/countlen" "$RSDIR/countlen-rs" -i
cmp_run "countlen -i f extra"    "$CORE/countlen" "$RSDIR/countlen-rs" -i "$inputs/edge_dna_n.fa" extra
cmp_run "countlen --"            "$CORE/countlen" "$RSDIR/countlen-rs" --
cmp_run "countlen -i /nope"      "$CORE/countlen" "$RSDIR/countlen-rs" -i /nope
cmp_run "version"                "$CORE/version"  "$RSDIR/version-rs"

# replaceu: the alphabet fold and the `_os_<i>_oe_<name>` marker restoreu reads.
# Its own indeterminate set is the same one, for the same reason -- it shares
# the reader, though not `getnumlen_nogap_countn`, so only the argument paths
# below can reach an undefined value.
for f in "$inputs"/*; do
  base=$(basename "$f")
  is_indeterminate "$base" && continue
  cmp_run "replaceu -i $base" "$CORE/replaceu" "$RSDIR/replaceu-rs" -i "$f"
done

for f in "$inputs"/*; do
  is_indeterminate "$(basename "$f")" && continue
  co=$("$CORE/replaceu" < "$f" 2>/dev/null); cs=$?
  ro=$("$RSDIR/replaceu-rs" < "$f" 2>/dev/null); rs=$?
  if [ "$co" = "$ro" ] && [ "$cs" = "$rs" ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL replaceu < $(basename "$f")"
    echo "  C : $(echo "$co" | head -c 200) (status $cs)"
    echo "  RS: $(echo "$ro" | head -c 200) (status $rs)"
  fi
done

# The seed offset shifts every `_os_` index, and `-D`/`-P` pick the fold that
# is applied -- both change the bytes downstream, so both are compared.
dnaf="$inputs/edge_dna_n.fa"
protf="$inputs/sample.fa"
cmp_run "replaceu -o 7"          "$CORE/replaceu" "$RSDIR/replaceu-rs" -o 7 -i "$protf"
cmp_run "replaceu -o -3"         "$CORE/replaceu" "$RSDIR/replaceu-rs" -o -3 -i "$protf"
cmp_run "replaceu -o junk"       "$CORE/replaceu" "$RSDIR/replaceu-rs" -o junk -i "$protf"
cmp_run "replaceu -D on protein" "$CORE/replaceu" "$RSDIR/replaceu-rs" -D -i "$protf"
cmp_run "replaceu -P on dna"     "$CORE/replaceu" "$RSDIR/replaceu-rs" -P -i "$dnaf"
cmp_run "replaceu -Do clustered" "$CORE/replaceu" "$RSDIR/replaceu-rs" -Do 2 -i "$protf"
cmp_run "replaceu -Z"            "$CORE/replaceu" "$RSDIR/replaceu-rs" -Z
cmp_run "replaceu positional"    "$CORE/replaceu" "$RSDIR/replaceu-rs" "$protf"
cmp_run "replaceu -i /nope"      "$CORE/replaceu" "$RSDIR/replaceu-rs" -i /nope

# f2cl: only the FASTA path is ported (see the binary's own header), so only
# the FASTA path is compared. `-n -1 -f -l 60` is what `--auto` passes.
for f in "$inputs"/*; do
  base=$(basename "$f")
  is_indeterminate "$base" && continue
  co=$("$CORE/f2cl" -n -1 -f -l 60 < "$f" 2>/dev/null); cs=$?
  ro=$("$RSDIR/f2cl-rs" -n -1 -f -l 60 < "$f" 2>/dev/null); rs=$?
  if [ "$co" = "$ro" ] && [ "$cs" = "$rs" ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL f2cl -f -l 60 < $base (status C=$cs RS=$rs)"
  fi
done

# The unwrapped branch is a different code path in the C, and differs from the
# wrapped one on an empty record: no sequence line at all against a blank one.
for f in "$inputs"/*; do
  base=$(basename "$f")
  is_indeterminate "$base" && continue
  co=$("$CORE/f2cl" -n -1 -f < "$f" 2>/dev/null); cs=$?
  ro=$("$RSDIR/f2cl-rs" -n -1 -f < "$f" 2>/dev/null); rs=$?
  if [ "$co" = "$ro" ] && [ "$cs" = "$rs" ]; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL f2cl -f (unwrapped) < $base (status C=$cs RS=$rs)"
  fi
done

cmp_stdout "f2cl -l 1"  "$CORE/f2cl" "$RSDIR/f2cl-rs" -n -1 -f -l 1 -i "$inputs/sample.fa"
cmp_stdout "f2cl -l 7"  "$CORE/f2cl" "$RSDIR/f2cl-rs" -n -1 -f -l 7 -i "$inputs/sample.fa"
cmp_stdout "f2cl -l 0"  "$CORE/f2cl" "$RSDIR/f2cl-rs" -n -1 -f -l 0 -i "$inputs/sample.fa"
cmp_stdout "f2cl -i /nope" "$CORE/f2cl" "$RSDIR/f2cl-rs" -n -1 -f -i /nope
cmp_stdout "f2cl -Z"    "$CORE/f2cl" "$RSDIR/f2cl-rs" -Z

# dndpre: needs an *aligned* input and writes `hat2` into the working
# directory, so each side runs in its own and the files are compared. Uses the
# committed alignments, which is the only part of this corpus that is aligned.
for f in "$GOLDEN"/aligned/*; do
  [ -f "$f" ] || continue
  cd1=$(mktemp -d); cd2=$(mktemp -d)
  ( cd "$cd1" && "$CORE/dndpre" -y hat2 -b 62 -M 2 -C 0 -i "$f" >/dev/null 2>/dev/null )
  ( cd "$cd2" && "$RSDIR/dndpre-rs" -y hat2 -b 62 -M 2 -C 0 -i "$f" >/dev/null 2>/dev/null )
  if cmp -s "$cd1/hat2" "$cd2/hat2"; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "FAIL dndpre $(basename "$f")"
    diff "$cd1/hat2" "$cd2/hat2" | head -6
  fi
  rm -rf "$cd1" "$cd2"
done

cmp_stdout "dndpre -y zzz"  "$CORE/dndpre" "$RSDIR/dndpre-rs" -y zzz -i "$GOLDEN/aligned/sample.fa"
cmp_stdout "dndpre 1 record" "$CORE/dndpre" "$RSDIR/dndpre-rs" -y hat2 -i "$GOLDEN/in/one.fa"

# A run that compared nothing is vacuous, not green.
if [ "$pass" -eq 0 ]; then
  echo "verify_utils: pass=0 -- nothing was compared, which is a harness failure" 1>&2
  exit 2
fi

echo "pass=$pass fail=$fail indeterminate=$indet"
[ "$fail" -eq 0 ] || exit 1
