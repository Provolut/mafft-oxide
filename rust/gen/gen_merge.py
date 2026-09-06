#!/usr/bin/env python3
"""Build a `--merge` case: align an input, then split it into sub-alignments.

    gen_merge.py <aligned.fa> <ngroups> [preservegaps] > _subalignmentstable

Groups are contiguous blocks of the aligned input, so they are already aligned
among themselves, which `--merge` requires.
"""
import sys

aligned, ngroups = sys.argv[1], int(sys.argv[2])
flags = sys.argv[3:]
preserve = "preservegaps" in flags
interleave = "interleave" in flags
n = sum(1 for line in open(aligned) if line.startswith('>'))
lines = []
if interleave:
    # Groups that the unconstrained UPGMA will almost certainly split, so the
    # builder has to reject merges and retry.
    for g in range(ngroups):
        ids = [i for i in range(1, n + 1) if (i - 1) % ngroups == g]
        if len(ids) < 2:
            continue
        if preserve:
            ids = [-i for i in ids]
        lines.append(" ".join(str(i) for i in ids))
else:
    per = n // ngroups
    p = 1
    for g in range(ngroups):
        last = n if g == ngroups - 1 else p + per - 1
        ids = list(range(p, last + 1))
        if len(ids) < 2:
            p = last + 1
            continue
        if preserve:
            ids = [-i for i in ids]
        lines.append(" ".join(str(i) for i in ids))
        p = last + 1
print("\n".join(lines))
