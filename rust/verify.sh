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
    if [ -f "$work/c/order" ] || [ -f "$work/r/order" ]; then
      if ! cmp -s "$work/c/order" "$work/r/order"; then
        fail=$((fail+1)); failed="$failed\n  $label (order differs)"
        continue
      fi
    fi
    pass=$((pass+1))
  done
done

echo "pass=$pass fail=$fail"
if [ $fail -ne 0 ]; then
  printf 'FAILED:%b\n' "$failed"
  exit 1
fi
