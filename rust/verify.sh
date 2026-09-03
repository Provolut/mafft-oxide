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

# Resolve the binaries to absolute paths and check they are runnable. Every
# comparison below runs from a scratch directory, so a relative path would not
# resolve -- and if *both* binaries fail to launch the run looks like a pass,
# which is how this harness once reported thousands of identical results while
# testing nothing at all.
abspath() {
  case "$1" in
    /*) printf '%s\n' "$1" ;;
    *)  printf '%s/%s\n' "$(pwd)" "$1" ;;
  esac
}
CBIN=$(abspath "$CBIN")
RBIN=$(abspath "$RBIN")
for b in "$CBIN" "$RBIN"; do
  if [ ! -x "$b" ]; then
    echo "verify: $b is not executable" 1>&2
    exit 2
  fi
done

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
  # dynamic matrix. With -O the C refuses ("headgp and tailgp must be 1")
  # after writing order, so only the three without it reach the warp DP.
  "-q 0 -E 2 -V -1.53 -s 0.8 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 2.00 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.8 -W 6 -C 0-0 -b 62 -g 0 -f -1.53 -Q 2.00 -h 0 -F -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.8 -W 6 -C 0-0 -b 62 -g 0 -f -1.53 -Q 2.00 -h 0 -X 0.1 -x 100"
  "-q 0 -E 2 -V -1.53 -s 0.8 -W 6 -C 0-0 -b 62 -g -0.1 -f -2.0 -Q 0.5 -h 0 -F -X 0.5 -x 100"
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

# An empty or missing corpus would make every loop below run zero times and
# report a clean pass. Refuse instead -- this harness has been vacuous once
# already (see the abspath check above).
n_corpus=$(find "$CORPUS" -maxdepth 1 -type f | wc -l)
if [ "$n_corpus" -lt 10 ]; then
  echo "verify: corpus at $CORPUS has only $n_corpus files; refusing to run" 1>&2
  exit 2
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

# -r (--oneiteration): leave-one-out refinement after the progressive pass.
#
# The C segfaults on part of this: `disttbfast()` frees `cpmxhist` and NULLs
# it before calling `dooneiteration`, which still passes `cpmxhist + l` to
# `A__align`, so any node past the first that reaches `A__align` writes
# through a small non-null offset from NULL. The port refuses exactly there,
# so a C crash paired with a Rust exit 1 is the expected outcome and counts as
# a pass; every other combination fails.
RARGS="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100 -r"
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  rm -rf "$work/c" "$work/r"
  mkdir -p "$work/c" "$work/r"
  ( cd "$work/c" && $CBIN $RARGS < "$f" > out.fa 2> err.log )
  cst=$?
  ( cd "$work/r" && $RBIN $RARGS < "$f" > out.fa 2> err.log )
  rst=$?
  label="$(basename "$f") [-r]"
  if [ $cst -ge 128 ] && [ $rst -eq 1 ]; then
    pass=$((pass+1))
    continue
  fi
  if [ $cst -ne $rst ]; then
    fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
    continue
  fi
  bad=""
  for o in out.fa order hat2; do
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

# --aamatrix (-b -1): the scoring matrix is read from `_aamtx` rather than
# built in. The variants cover a permuted residue order, a `frequency` line
# and `norescale` -- and, since `constants()` starts with
# `if( nblosum < 0 ) dorp = 'p';`, the nucleotide half of the corpus is
# expected to be refused by both binaries as protein with illegal characters.
MARGS="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b -1 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  for variant in plain perm freq norescale permfreq nocomment; do
    rm -rf "$work/c" "$work/r"
    mkdir -p "$work/c" "$work/r"
    python3 "$here/gen_aamtx.py" "$variant" > "$work/c/_aamtx" || continue
    cp "$work/c/_aamtx" "$work/r/_aamtx"
    ( cd "$work/c" && $CBIN $MARGS < "$f" > out.fa 2> err.log )
    cst=$?
    ( cd "$work/r" && $RBIN $MARGS < "$f" > out.fa 2> err.log )
    rst=$?
    label="$(basename "$f") [--aamatrix $variant]"
    if [ $cst -ne $rst ]; then
      fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
      continue
    fi
    bad=""
    for o in out.fa order hat2; do
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

# --treein's generated trees: `pileup` and `shuffle` build a comb instead of
# clustering (shuffle permutes the input first with glibc's rand, which the
# port reproduces), `nodepair` is a keyword the C accepts and then rejects,
# and `auto` / `test` sample pairwise scores first.
#
# `auto` and `test` are refused by the port whenever a picked sequence is
# gap-free and of maximal length: `varpairscore` gives each row exactly
# `nlenmax` bytes and `gappick0` writes a terminator past the end, which the
# next malloc overwrites, after which the C aligns whatever follows the row.
# A C crash, or a Rust refusal naming that, is the expected outcome there.
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  n=$(grep -c '^>' "$f")
  [ "$n" -ge 4 ] || continue
  np=4
  [ "$n" -lt 4 ] && np=$n
  for kw in "pileup" "shuffle" "shuffle 42" "shuffle 12345" "nodepair" "auto 0 $np" "test 1 $np" \
            "compact 10GB" "compact 1gb" "compact 5m" "compact 100kb" "compact 3" \
            "compact 100" "compact 1kb" "very compact" "stepadd" "youngestlinkage"; do
    for extra in "" "-t"; do
      rm -rf "$work/c" "$work/r"
      mkdir -p "$work/c" "$work/r"
      printf '%s\n' "$kw" > "$work/c/_guidetree"
      cp "$work/c/_guidetree" "$work/r/_guidetree"
      GARGS="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100 -U $extra"
      ( cd "$work/c" && $CBIN $GARGS < "$f" > out.fa 2> err.log )
      cst=$?
      ( cd "$work/r" && $RBIN $GARGS < "$f" > out.fa 2> err.log )
      rst=$?
      label="$(basename "$f") [--treein $kw $extra]"
      # The C aborting, or the port refusing the undefined varpairscore case,
      # is expected and counts as a pass -- but only for auto/test.
      case "$kw" in
        auto*|test*)
          if [ $cst -ge 128 ]; then pass=$((pass+1)); continue; fi
          if [ $rst -eq 1 ] && grep -qa "writes its terminator outside" "$work/r/err.log"; then
            pass=$((pass+1)); continue
          fi
          ;;
      esac
      if [ $cst -ne $rst ]; then
        fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
        continue
      fi
      bad=""
      for o in out.fa order _guidetree infile.tree; do
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

# --text / --textmatrix (-b -2): the alphabet becomes every byte, so the
# corpus does not apply -- gen_text.py makes inputs out of characters the
# ordinary reader would refuse, plus the `_aamtx` that scores them.
# `minusone` sets a letter's diagonal to -1, which is how `constants()` spells
# "no score given", and the C then refuses the run; both binaries must agree.
for variant in full sparse nomatrix unscored minusone; do
  for seed in 1 2 3 4 5 6 7 8; do
    for extra in "" "-a" "-F"; do
      rm -rf "$work/c" "$work/r"
      mkdir -p "$work/c" "$work/r"
      python3 "$here/gen_text.py" "$seed" "$variant" "$work/c" || continue
      cp "$work/c/in.fa" "$work/r/in.fa"
      [ -f "$work/c/_aamtx" ] && cp "$work/c/_aamtx" "$work/r/_aamtx"
      TARGS="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b -2 -P -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100 $extra"
      ( cd "$work/c" && $CBIN $TARGS < in.fa > out.fa 2> err.log )
      cst=$?
      ( cd "$work/r" && $RBIN $TARGS < in.fa > out.fa 2> err.log )
      rst=$?
      label="text/$variant/$seed [$extra]"
      if [ $cst -ne $rst ]; then
        fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
        continue
      fi
      bad=""
      for o in out.fa order hat2; do
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

# --merge (-H): sub-alignments that must come out monophyletic. The groups
# have to be aligned already, so the C is run once to produce an alignment and
# gen_merge.py splits it -- `interleave` picks groups the unconstrained UPGMA
# would split, which is what exercises the veto-and-rescan path.
HARGS="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  n=$(grep -c '^>' "$f")
  [ "$n" -ge 6 ] || continue
  rm -rf "$work/pre"; mkdir -p "$work/pre"
  ( cd "$work/pre" && $CBIN $HARGS < "$f" > aligned.fa 2> /dev/null ) || continue
  [ -s "$work/pre/aligned.fa" ] || continue
  for spec in "2" "3" "2 interleave" "3 interleave" "4 interleave" "2 interleave preservegaps"; do
    python3 "$here/gen_merge.py" "$work/pre/aligned.fa" $spec > "$work/pre/st.txt" || continue
    [ -s "$work/pre/st.txt" ] || continue
    rm -rf "$work/c" "$work/r"
    mkdir -p "$work/c" "$work/r"
    cp "$work/pre/st.txt" "$work/c/_subalignmentstable"
    cp "$work/pre/st.txt" "$work/r/_subalignmentstable"
    ( cd "$work/c" && $CBIN $HARGS -H 0 < "$work/pre/aligned.fa" > out.fa 2> err.log )
    cst=$?
    ( cd "$work/r" && $RBIN $HARGS -H 0 < "$work/pre/aligned.fa" > out.fa 2> err.log )
    rst=$?
    label="$(basename "$f") [-H $spec]"
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

# --add / --addprofile (-I n, -K): the first njob-n sequences are an existing
# alignment whose columns must survive. The input has to be built from a real
# alignment, so the C is run once first to produce one.
AARGS="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -F -X 0.1 -x 100"
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  n=$(grep -c '^>' "$f")
  [ "$n" -ge 4 ] || continue
  rm -rf "$work/pre"; mkdir -p "$work/pre"
  ( cd "$work/pre" && $CBIN $AARGS < "$f" > aligned.fa 2> /dev/null ) || continue
  [ -s "$work/pre/aligned.fa" ] || continue
  for nadd in 1 $((n/2)); do
    [ "$nadd" -ge 1 ] || continue
    for mode in add profile; do
      python3 "$here/gen_add.py" "$work/pre/aligned.fa" "$nadd" "$mode" > "$work/pre/in.fa" || continue
      kflag=""
      [ "$mode" = "add" ] && kflag="-K"
      # -Y (--keeplength) rides on --add: it cuts the columns the alignment
      # inserted into the reference back out, and -Z / -z say what that cost.
      for extra in "" "-Y" "-Y -Z" "-Y -z" "-Y -Z -n"; do
        rm -rf "$work/c" "$work/r"
        mkdir -p "$work/c" "$work/r"
        ( cd "$work/c" && $CBIN $AARGS -I $nadd $kflag $extra < "$work/pre/in.fa" > out.fa 2> err.log )
        cst=$?
        ( cd "$work/r" && $RBIN $AARGS -I $nadd $kflag $extra < "$work/pre/in.fa" > out.fa 2> err.log )
        rst=$?
        label="$(basename "$f") [-I $nadd $mode $extra]"
        if [ $cst -ne $rst ]; then
          fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
          continue
        fi
        bad=""
        for o in out.fa order _deletelist _deletemap; do
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
done

# -l (--anchors): anchors read from `_externalanchors`, which pin an interval
# of one sequence to an interval of another and cut the DP into blocks around
# them. gen_anchors.py writes both consistent anchors and deliberately
# conflicting ones, since the conflicting ones are what drive
# `checkanchors_strongestfirst`.
#
# The C segfaults on part of this too. `readexternalanchors` range-checks only
# the two *start* positions, and even those with `>` rather than `>=`; the end
# positions are never checked. An end position past the sequence's residue
# count -- which is what a user gets by counting on their own FASTA when MAFFT
# drops characters from it, `*` for instance -- reads `recountpositions`' map
# above the entries it filled, so the anchor comes back ending before it
# starts. `Falign_givenanchors` lays that inversion into its cut list and then
# copies the block with a negative `strncpy` length. As with -r, a C crash
# paired with a Rust exit 1 is the expected outcome and counts as a pass.
LARGS_BASE="-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100"
LARGSETS=(
  "$LARGS_BASE -l"
  "$LARGS_BASE -l -u"
  "-q 0 -E 1 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100 -l"
  "-q 0 -E 2 -V -2.00 -s 0.0 -W 6 -O -C 0-0 -b 62 -g -0.05 -f -2.00 -Q 100.0 -h -0.1 -X 0.0 -x 100 -l"
)
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  n=$(grep -c '^>' "$f")
  [ "$n" -ge 2 ] || continue
  # seed + count, and whether the intervals may conflict
  for spec in "1 4" "2 10" "3 6 overlap" "4 20 overlap"; do
    for i in "${!LARGSETS[@]}"; do
      rm -rf "$work/c" "$work/r"
      mkdir -p "$work/c" "$work/r"
      python3 "$here/gen_anchors.py" "$f" $spec > "$work/c/_externalanchors" || continue
      [ -s "$work/c/_externalanchors" ] || continue
      cp "$work/c/_externalanchors" "$work/r/_externalanchors"
      ( cd "$work/c" && $CBIN ${LARGSETS[$i]} < "$f" > out.fa 2> err.log )
      cst=$?
      ( cd "$work/r" && $RBIN ${LARGSETS[$i]} < "$f" > out.fa 2> err.log )
      rst=$?
      label="$(basename "$f") [anchors $spec, argset $i]"
      if [ $cst -ge 128 ] && [ $rst -eq 1 ]; then
        pass=$((pass+1))
        continue
      fi
      if [ $cst -ne $rst ]; then
        fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
        continue
      fi
      bad=""
      for o in out.fa order hat2; do
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

# -M (--memsave): the linear-space aligner in MSalignmm.c. The C also reaches
# it without being asked, whenever a row passes 30000 columns (disttbfast.c
# :2836), and `alg` is a global that never switches back -- so this is not a
# rarely-taken branch. It is not equivalent to A__align either: MSalignmm.c
# compiles with USE_PENALTY_EX 0 where Salignmm.c uses 1, which is why forcing
# it changes the alignment on 17 of the 86 corpus inputs.
#
# The FFT branch under -M is `Falign_udpari_long`, a different anchoring pass,
# and is not ported -- so these argument vectors deliberately leave -F off.
MARGSETS=(
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100 -M -B"
  "-q 0 -E 1 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100 -M -B"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100 -M -B"
  "-q 0 -E 2 -V -2.00 -s 0.0 -W 6 -O -C 0-0 -b 62 -g -0.05 -f -2.00 -Q 100.0 -h -0.1 -X 0.0 -x 100 -M -B"
  "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 100 -M -B -u"
)
for f in "$CORPUS"/*; do
  [ -f "$f" ] || continue
  for i in "${!MARGSETS[@]}"; do
    rm -rf "$work/c" "$work/r"
    mkdir -p "$work/c" "$work/r"
    ( cd "$work/c" && $CBIN ${MARGSETS[$i]} < "$f" > out.fa 2> err.log )
    cst=$?
    ( cd "$work/r" && $RBIN ${MARGSETS[$i]} < "$f" > out.fa 2> err.log )
    rst=$?
    label="$(basename "$f") [memsave argset $i]"
    if [ $cst -ne $rst ]; then
      fail=$((fail+1)); failed="$failed\n  $label (exit $cst vs $rst)"
      continue
    fi
    bad=""
    for o in out.fa order hat2; do
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

echo "pass=$pass fail=$fail"
if [ $fail -ne 0 ]; then
  printf 'FAILED:%b\n' "$failed"
  exit 1
fi
