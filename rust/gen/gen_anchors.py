#!/usr/bin/env python3
"""Generate `_externalanchors` for a FASTA input.

    gen_anchors.py <input.fa> <seed> <nanchors> [overlap]

Anchors pin an interval of one sequence to an interval of another. Positions
are 1-origin over the *ungapped* sequence, as MAFFT wants them. With `overlap`
the intervals are allowed to conflict, which is what exercises
`checkanchors_strongestfirst`.
"""
import random
import sys

path, seed, nanchors = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
overlap = len(sys.argv) > 4 and sys.argv[4] == "overlap"
random.seed(seed)

seqs, cur = [], None
for line in open(path):
    if line.startswith('>'):
        if cur is not None:
            seqs.append(cur)
        cur = ''
    else:
        cur += line.strip()
seqs.append(cur)
lens = [len(s.replace('-', '')) for s in seqs]
n = len(seqs)
if n < 2:
    sys.exit(0)

print("# generated anchors")
made = 0
tries = 0
while made < nanchors and tries < nanchors * 40:
    tries += 1
    i = random.randrange(n - 1)
    j = random.randrange(i + 1, n)
    li, lj = lens[i], lens[j]
    w = min(li, lj) // 6
    if w < 2:
        continue
    if overlap:
        si = random.randrange(1, max(2, li - w))
        sj = random.randrange(1, max(2, lj - w))
    else:
        # Keep the anchors ordered and disjoint, one per sixth.
        k = made % 3
        si = 1 + k * (li // 3)
        sj = 1 + k * (lj // 3)
        if si + w > li or sj + w > lj:
            continue
    ei, ej = si + w - 1, sj + w - 1
    if ei > li or ej > lj:
        continue
    print("%d %d %d %d %d %d %d" % (i + 1, j + 1, si, ei, sj, ej, random.randrange(1, 100)))
    made += 1
