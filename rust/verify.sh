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
)

gen_corpus() {
  local d="$1"
  mkdir -p "$d"
  local sample="$here/../test/sample"
  cp "$sample" "$d/sample.fa"

  python3 - "$sample" "$d" <<'PY'
import random, sys
src, out = sys.argv[1], sys.argv[2]

recs, name, seq = [], None, []
for line in open(src):
    line = line.rstrip('\n')
    if line.startswith('>'):
        if name is not None:
            recs.append((name, ''.join(seq)))
        name, seq = line[1:], []
    else:
        seq.append(line)
recs.append((name, ''.join(seq)))

def write(path, rs):
    with open(path, 'w') as f:
        for n, s in rs:
            f.write('>%s\n%s\n' % (n, s))

rng = random.Random(20240607)

# Random subsets of the real sample, sizes 2..len(recs).
for i in range(40):
    k = rng.randint(2, len(recs))
    write('%s/subset_%02d.fa' % (out, i), rng.sample(recs, k))

# Truncated sequences: exercises short profiles and the `too short for a
# 6-mer` sentinel in seq_grp().
for L in (1, 2, 3, 5, 6, 7, 12, 40):
    rs = [(n, s[:L]) for n, s in recs[:8]]
    write('%s/trunc_%02d.fa' % (out, L), rs)

# Two sequences only, and duplicated sequences (ties everywhere).
write('%s/pair.fa' % out, recs[:2])
write('%s/dup.fa' % out, [(('d%d' % i), recs[0][1]) for i in range(6)])
write('%s/dup_mixed.fa' % out,
      [('a', recs[0][1]), ('b', recs[0][1]), ('c', recs[1][1]),
       ('d', recs[1][1]), ('e', recs[2][1])])

# Ragged lengths, including one very short sequence among long ones.
write('%s/ragged.fa' % out,
      [(n, s[:rng.randint(3, len(s))]) for n, s in recs[:12]])

# Ambiguity codes and pre-existing gaps in the input.
alpha = 'ARNDCQEGHILKMFPSTWYVBZX'
rs = []
for i, (n, s) in enumerate(recs[:10]):
    t = list(s[:120])
    for _ in range(12):
        t[rng.randrange(len(t))] = rng.choice(alpha)
    rs.append((n, ''.join(t)))
write('%s/ambig.fa' % out, rs)

rs = []
for i, (n, s) in enumerate(recs[:10]):
    t = list(s[:150])
    for _ in range(20):
        t[rng.randrange(len(t))] = '-'
    rs.append((n, ''.join(t)))
write('%s/gappy.fa' % out, rs)

# Random sequences: no phylogenetic signal, so lots of near-ties in the tree.
for i, (nseq, L) in enumerate(((5, 30), (12, 60), (25, 45), (60, 25), (30, 200))):
    rs = [('r%d' % j, ''.join(rng.choice(alpha[:20]) for _ in range(L)))
          for j in range(nseq)]
    write('%s/random_%02d.fa' % (out, i), rs)

# Long names (the C truncates at B-2 = 254 bytes) and odd FASTA formatting.
with open('%s/longname.fa' % out, 'w') as f:
    for i, (n, s) in enumerate(recs[:4]):
        f.write('>%s\n' % (('N' * 300) + str(i)))
        for k in range(0, len(s), 37):
            f.write(s[k:k+37] + '\n')

# Lower case, and '*' which the reader turns into '-'.
write('%s/lower.fa' % out, [(n, s.lower()) for n, s in recs[:8]])
write('%s/stars.fa' % out,
      [(n, s[:100].replace('A', '*')) for n, s in recs[:8]])

# Tandem repeats with varying copy number, shuffled domains and inversions.
# These are the inputs on which the FFT scan proposes mutually inconsistent
# anchors, so blockAlign2 has to prune -- the only place its DP quirks
# (a penalised `maxj` carried into the next cell, a `result1[]` element left
# over from the previous call) can change the answer.
for idx, (motiflen, ncopies, nseq, mut) in enumerate(
        [(15, (6, 20), 16, 0.10), (25, (4, 12), 20, 0.15), (40, (3, 8), 12, 0.10),
         (50, (2, 6), 10, 0.05), (20, (8, 16), 30, 0.20)]):
    motif = ''.join(rng.choice(alpha[:20]) for _ in range(motiflen))
    rs = []
    for i in range(nseq):
        k = rng.randint(*ncopies)
        s = ''
        for _ in range(k):
            t = list(motif)
            for _ in range(max(1, int(motiflen * mut))):
                t[rng.randrange(motiflen)] = rng.choice(alpha[:20])
            s += ''.join(t)
        rs.append(('r%02d' % i, s))
    write('%s/tandem_%d.fa' % (out, idx), rs)

domains = [''.join(rng.choice(alpha[:20]) for _ in range(70)) for _ in range(4)]
rs = []
for i in range(14):
    order = list(range(4))
    rng.shuffle(order)
    t = list(''.join(domains[d] for d in order))
    for _ in range(len(t) // 15):
        t[rng.randrange(len(t))] = rng.choice(alpha[:20])
    rs.append(('dom%02d' % i, ''.join(t)))
write('%s/shuffled_domains.fa' % out, rs)

base = ''.join(rng.choice(alpha[:20]) for _ in range(300))
rs = []
for i in range(12):
    a, b = sorted(rng.sample(range(300), 2))
    t = list(base[:a] + base[a:b][::-1] + base[b:])
    for _ in range(20):
        t[rng.randrange(len(t))] = rng.choice(alpha[:20])
    rs.append(('inv%02d' % i, ''.join(t)))
write('%s/inversions.fa' % out, rs)
PY
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
