#!/usr/bin/env python3
"""Generate a `--text` (-b -2) case: arbitrary-byte sequences plus `_aamtx`.

    gen_text.py <seed> <variant> <outdir>

variants: full (a score for every pair of letters used), sparse (only the
diagonal), nomatrix (no `_aamtx` at all -- the C warns and uses its default
1.0/-1.0 matrix), unscored (a letter is left out of the matrix, which still
leaves it the default 1.0 on the diagonal), minusone (a letter's diagonal is
explicitly set to -1, which is how `constants()` spells "no score given" -- the
C then refuses the run).
"""
import os
import random
import sys

seed, variant, outdir = int(sys.argv[1]), sys.argv[2], sys.argv[3]
random.seed(seed)
os.makedirs(outdir, exist_ok=True)

# Bytes outside the amino alphabet, so the ordinary reader would reject them.
alpha = "abcdefghij0123456789!$%&@+/:;?[]^_{}~"
alpha = alpha[: 8 + seed % 20]
nseq = 4 + seed % 5
recs = []
for i in range(nseq):
    ln = 30 + random.randrange(40)
    recs.append((f">t{i}", "".join(random.choice(alpha) for _ in range(ln))))
with open(os.path.join(outdir, "in.fa"), "w") as fh:
    fh.write("".join(n + "\n" + s + "\n" for n, s in recs))

used = sorted(set("".join(s for _, s in recs)))
odd = None
if variant == "unscored":
    used = used[:-1]  # leave the last letter out of the matrix entirely
elif variant == "minusone":
    odd = used[-1]
lines = ["# generated text matrix"]
for i, a in enumerate(used):
    for b in used[i:]:
        if a == odd or b == odd:
            # Only the diagonal, and explicitly -1: "no score given".
            if a == b:
                lines.append("0x%02x 0x%02x -1" % (ord(a), ord(b)))
            continue
        if variant in ("sparse", "unscored") and a != b:
            continue
        lines.append("0x%02x 0x%02x %d" % (ord(a), ord(b), 10 if a == b else random.randrange(-4, 3)))
path = os.path.join(outdir, "_aamtx")
if variant == "nomatrix":
    if os.path.exists(path):
        os.remove(path)
else:
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
