#!/usr/bin/env python3
"""Generate the differential-test corpus shared by verify.sh and verify_tbfast.sh.

    gencorpus.py <test/sample> <outdir> [full|small]

Every case exists to reach a branch that ordinary input does not: sequences too
short for a 6-mer, exact duplicates (which tie everywhere in the guide tree),
`X` (stripped by the score-only variants), tandem repeats and inversions (the
only inputs on which the FFT scan proposes mutually inconsistent anchors, so
`blockAlign2` has to prune), and both sequence types, since `constants()`, the
k-tuple distance and the FFT vector all branch on protein vs nucleotide.

`small` emits a reduced set for the all-pairs harness, which is far slower.
"""
import random
import sys

src, out = sys.argv[1], sys.argv[2]
mode = sys.argv[3] if len(sys.argv) > 3 else 'full'

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


alpha = 'ARNDCQEGHILKMFPSTWYVBZX'

if mode == 'full':
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
else:
    rng = random.Random(31337)
    for i in range(12):
        k = rng.randint(2, len(recs))
        write('%s/subset_%02d.fa' % (out, i), rng.sample(recs, k))
    for L in (2, 6, 12, 40):
        write('%s/trunc_%02d.fa' % (out, L), [(n, s[:L]) for n, s in recs[:6]])
    write('%s/pair.fa' % out, recs[:2])
    write('%s/dup.fa' % out, [(('d%d' % i), recs[0][1]) for i in range(5)])
    # X is stripped for the score-only variants, so it exercises `thereisx`.
    rs = []
    for n, s in recs[:8]:
        t = list(s[:150])
        for _ in range(10):
            t[rng.randrange(len(t))] = 'X'
        rs.append((n, ''.join(t)))
    write('%s/withx.fa' % out, rs)
    rs = []
    for n, s in recs[:8]:
        t = list(s[:150])
        for _ in range(10):
            t[rng.randrange(len(t))] = rng.choice(alpha)
        rs.append((n, ''.join(t)))
    write('%s/ambig.fa' % out, rs)
    motif = ''.join(rng.choice(alpha[:20]) for _ in range(30))
    rs = []
    for i in range(10):
        s = ''
        for _ in range(rng.randint(2, 6)):
            t = list(motif)
            for _ in range(3):
                t[rng.randrange(len(t))] = rng.choice(alpha[:20])
            s += ''.join(t)
        rs.append(('r%02d' % i, s))
    write('%s/tandem.fa' % out, rs)
    for i, (nseq, L) in enumerate(((4, 40), (15, 70))):
        write('%s/random_%d.fa' % (out, i),
              [('q%d' % j, ''.join(rng.choice(alpha[:20]) for _ in range(L)))
               for j in range(nseq)])

# ------------------------------------------------------------- nucleotide
# The reader decides the type from the input itself (>75% atgcun), so these
# files need no flag to take the DNA branch of constants(), the base-4 k-tuple
# distance and the single-indicator FFT vector.
rng = random.Random(20260903)


def mutate(s, sub=0.12, indel=0.03):
    t = []
    for ch in s:
        r = rng.random()
        if r < indel * 0.5:
            continue
        if r < indel:
            t.append(rng.choice('acgt'))
        t.append(rng.choice('acgt') if rng.random() < sub else ch)
    return ''.join(t)


fams = ((8, 300), (20, 150), (4, 600), (30, 80)) if mode == 'full' else ((8, 300), (20, 150))
for k, (n, L) in enumerate(fams):
    anc = ''.join(rng.choice('acgt') for _ in range(L))
    write('%s/d_fam_%d.fa' % (out, k), [('s%02d' % i, mutate(anc)) for i in range(n)])
write('%s/d_rand.fa' % out,
      [('r%d' % i, ''.join(rng.choice('acgt') for _ in range(120))) for i in range(10)])
# Upper case on input: the reader lower-cases nucleotides and upper-cases protein.
anc = ''.join(rng.choice('acgt') for _ in range(200))
write('%s/d_upper.fa' % out, [('u%d' % i, mutate(anc).upper()) for i in range(6)])
# 'u' rather than 't', scored by a copied row of the matrix.
anc = ''.join(rng.choice('acgu') for _ in range(200))
write('%s/d_rna.fa' % out,
      [('n%d' % i, mutate(anc.replace('u', 't')).replace('t', 'u')) for i in range(6)])
# IUPAC ambiguity codes, which exercise ambiguousscore().
anc = ''.join(rng.choice('acgt') for _ in range(250))
rs = []
for i in range(8):
    t = list(mutate(anc))
    for _ in range(25):
        t[rng.randrange(len(t))] = rng.choice('rykmswbdhvn')
    rs.append(('a%d' % i, ''.join(t)))
write('%s/d_ambig.fa' % out, rs)
motif = ''.join(rng.choice('acgt') for _ in range(40))
rs = []
for i in range(10):
    rs.append(('t%02d' % i, ''.join(mutate(motif, 0.08, 0.0)
                                    for _ in range(rng.randint(2, 6)))))
write('%s/d_tandem.fa' % out, rs)
anc = ''.join(rng.choice('acgt') for _ in range(150))
write('%s/d_gapped.fa' % out,
      [('g%d' % i, mutate(anc)[:60] + '-' * rng.randint(1, 10) + mutate(anc)[60:])
       for i in range(6)])
for L in (2, 5, 6, 12):
    anc = ''.join(rng.choice('acgt') for _ in range(L))
    write('%s/d_trunc_%02d.fa' % (out, L),
          [('s%d' % i, mutate(anc, 0.2, 0.0)) for i in range(5)])
anc = ''.join(rng.choice('acgt') for _ in range(200))
write('%s/d_pair.fa' % out, [('p0', anc), ('p1', mutate(anc))])
write('%s/d_dup.fa' % out, [('d%d' % i, anc) for i in range(5)])
write('%s/d_ragged.fa' % out,
      [('q%d' % i, mutate(anc)[:rng.randint(20, 200)]) for i in range(8)])
