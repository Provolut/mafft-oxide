#!/usr/bin/env python3
"""Generate `_aamtx` files for the --aamatrix (-b -1) tests.

    gen_aamtx.py <variant> > _aamtx

The variants exercise every branch of `loadaamtx()`: the residue order in the
header need not be the canonical one (the matrix is permuted into ours and only
the lower triangle is read, so a value is looked up at (max, min) of the two
mapped indices); a `frequency` line overrides the background frequencies; a
`norescale` line stops `constants()` recentring the matrix and makes it scale
by a flat 600 instead.
"""
import random
import sys

AA = "ARNDCQEGHILKMFPSTWYV"
variant = sys.argv[1] if len(sys.argv) > 1 else "plain"

random.seed(7)
base = {}
for a in AA:
    for b in AA:
        if (b, a) in base:
            base[(a, b)] = base[(b, a)]
        else:
            base[(a, b)] = 8 if a == b else random.randint(-4, 4)
freq = [round(random.uniform(0.02, 0.09), 4) for _ in range(20)]

order = list(AA)
if "perm" in variant:
    random.shuffle(order)

out = []
if variant != "nocomment":
    out.append("# generated test matrix")
out.append("   " + "  ".join(order))
for i, a in enumerate(order):
    out.append(a + " " + " ".join("%3d" % base[(a, order[j])] for j in range(i + 1)))
if "freq" in variant:
    out.append("frequency " + " ".join(str(x) for x in freq))
if "norescale" in variant:
    out.append("norescale")
print("\n".join(out))
