#!/bin/bash
# Differential test: run the C `disttbfast` and `disttbfast-rs` on the same
# inputs with the same arguments and require byte-identical stdout and `order`.
#
#   verify.sh <path-to-C-disttbfast> <path-to-rust-disttbfast-rs> [corpus-dir]
#
# With no corpus directory, a corpus is generated from ../test/sample plus a
# set of edge cases (see gen_corpus below).

set -u
CBIN="${1:?usage: verify.sh <c-disttbfast> <rust-disttbfast-rs> [corpus-dir]}"
RBIN="${2:?usage: verify.sh <c-disttbfast> <rust-disttbfast-rs> [corpus-dir]}"
CORPUS="${3:-}"

here=$(cd "$(dirname "$0")" && pwd)
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

# The argument vectors the mafft wrapper produces for the ported strategies.
ARGSETS=(
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100"
  "-q 0 -E 1 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100"
  "-q 0 -E 3 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100"
  "-q 0 -E 2 -V -2.00 -s 0.0 -W 6 -O -C 0-0 -b 62 -g -0.05 -f -2.00 -Q 100.0 -h -0.1 -X 0.0 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.5 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100 -u"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100 -n"
  # FFT anchoring (FFT-NS-1 / FFT-NS-2, i.e. plain mafft and fftns)
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 1 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -2.00 -s 0.0 -W 6 -O -C 0-0 -b 62 -g -0.05 -f -2.00 -Q 100.0 -h -0.1 -F -X 0.0 -x 100"
  "-q 0 -E 3 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.5 -x 100 -u"
  # -G forces the FFT path even when the sequences are short for it
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -G -X 0.1 -x 100"
  # non-default FFT window and threshold
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -w 40 -X 0.1 -x 100"
  # --anysymbol: a no-op for protein, but it rewrites the `n` row and column of
  # the nucleotide matrix (nscore())
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -: -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  # the other built-in BLOSUM tables (--bl); -b is ignored for nucleotides
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 30 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 45 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 80 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 100 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  # JTT (--jtt) and transmembrane (--tm); a negative PAM distance skips the
  # recentring step (makeaverage0)
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -j 200 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -j 100 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -m 200 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -j -200 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  # --fmodel (-a): residue frequencies taken from the input rather than the
  # table's, which changes the matrix for all three model families
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -a -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -j 200 -a -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  # --legacygapcost (-L): drops the gap-frequency weighting, but only on the
  # branch of A__align that Falign takes, so both -F and non-FFT matter
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -L -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -L -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100"
  # -y (distance output) and -T (no alignment), including the combination that
  # stops before `order` is written
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -y h -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -y h -T -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  # 10-mer nucleotide tuples (--tuplesize 10); protein refuses it, as in the C
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 10 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  # --allowshift / --unalignlevel: the shift (warp) DP and the per-node
  # dynamic matrix
  "-q 0 -E 2 -V -1.53 -s 0.8 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 2.00 -h 0 -F -X 0.1 -x 100"
  # -t (--treeout): the Newick file, including the -n name rewriting
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -t -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -t -T -n -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  # -^ (--nodeout): the tree plus per-sequence density and per-node membership
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -^ -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
  # -q N (sparse pickup): keep a subset, write the rest to `notused`, and stop
  "-q 8 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
)

gen_corpus() {
  local d="$1"
  mkdir -p "$d"
  cp "$here/../test/sample" "$d/sample.fa"
  cp "$here/../test/samplerna" "$d/samplerna.fa"
  python3 "$here/gencorpus.py" "$here/../test/sample" "$d" full
}

if [ -z "$CORPUS" ]; then
  CORPUS="$work/corpus"
  gen_corpus "$CORPUS"
fi

pass=0
fail=0
failed=""
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  for i in "${!ARGSETS[@]}"; do
    args="${ARGSETS[$i]}"
    rm -rf "$work/c" "$work/r"
    mkdir -p "$work/c" "$work/r"
    ( cd "$work/c" && $CBIN $args < "$f" > out.fa 2> err.log )
    cst=$?
    ( cd "$work/r" && $RBIN $args < "$f" > out.fa 2> err.log )
    rst=$?
    label="$(basename "$f") [argset $i]"
    if [ $cst -ne $rst ]; then
      fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
      continue
    fi
    if ! cmp -s "$work/c/out.fa" "$work/r/out.fa"; then
      fail=$((fail+1)); failed="$failed\n  $label (stdout differs)"
      continue
    fi
    bad=""
    for out in order hat2 infile.tree notused; do
      if [ -f "$work/c/$out" ] || [ -f "$work/r/$out" ]; then
        cmp -s "$work/c/$out" "$work/r/$out" || bad="$bad $out"
      fi
    done
    if [ -n "$bad" ]; then
      fail=$((fail+1)); failed="$failed\n  $label (differs:$bad)"
      continue
    fi
    pass=$((pass+1))
  done
done

# --treein (-U): a user guide tree read from `_guidetree`, which means no
# distance matrix is built at all. The tree must match the input's sequence
# count, so it is generated per input rather than being an ARGSETS entry.
UARGS="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100 -U"
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  n=$(grep -c '^>' "$f")
  [ "$n" -ge 2 ] || continue
  for shape in ladder balanced; do
    for extra in "" "-t"; do
      rm -rf "$work/c" "$work/r"
      mkdir -p "$work/c" "$work/r"
      python3 "$here/gen_treein.py" "$n" "$shape" > "$work/c/_guidetree"
      cp "$work/c/_guidetree" "$work/r/_guidetree"
      ( cd "$work/c" && $CBIN $UARGS $extra < "$f" > out.fa 2> err.log )
      cst=$?
      ( cd "$work/r" && $RBIN $UARGS $extra < "$f" > out.fa 2> err.log )
      rst=$?
      label="$(basename "$f") [-U $shape $extra]"
      if [ $cst -ne $rst ]; then
        fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
        continue
      fi
      bad=""
      for o in out.fa order infile.tree; do
        if [ -f "$work/c/$o" ] || [ -f "$work/r/$o" ]; then
          cmp -s "$work/c/$o" "$work/r/$o" || bad="$bad $o"
        fi
      done
      if [ -n "$bad" ]; then
        fail=$((fail+1)); failed="$failed\n  $label (differs:$bad)"
      else
        pass=$((pass+1))
      fi
    done
  done
done

echo "pass=$pass fail=$fail"
if [ $fail -ne 0 ]; then
  printf 'FAILED:%b\n' "$failed"
  exit 1
fi
