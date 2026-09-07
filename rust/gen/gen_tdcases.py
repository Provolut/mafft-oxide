#!/usr/bin/env python3
"""Emit the guide trees `verify_td.sh` feeds to both branch-weighting families.

One case per line, tab-separated, in the order `gen_tdharness.c` and
`examples/tdref.rs` read them:

    nseq  then nseq-1 x ( members0 members1 len0 len1 )

A member list is comma-separated; both sides append the C's -1 terminator
themselves. The trees are merge sequences, exactly what `topol` holds: at every
step two existing clusters join, and each side's list is that cluster's members.

Branch lengths deliberately include zeros and negatives. They are not noise:
`checkMinusLength` clamps them to MINLEN in place, and a zero-length arm is what
drives `syntheticLength`'s `nanflag` and `calcW`'s two early returns, which are
the only paths in this family that are not plain arithmetic.
With --small, only the trees up to twelve sequences are emitted. That is what
`tests/golden/td/` holds: the fixtures replay in `cargo test`, and the per-case
output is quadratic in the number of sequences, so the large trees would add
megabytes to the repository to re-check arithmetic the small ones already cover.
The full corpus, up to 36 sequences, is what `verify_td.sh` runs against a live
C.
"""
import random
import sys


def merge_tree(rng, n, shape):
    """One merge sequence over `n` leaves, as (topol, len) rows."""
    clusters = [[i] for i in range(n)]
    topol, lens = [], []
    for _ in range(n - 1):
        if shape == "ladder":
            # Always grow the same cluster: maximum depth, so the recursions in
            # syntheticLength and the two _rec walks go as deep as they can.
            a = 0
            b = 1 if len(clusters) > 1 else 0
        elif shape == "balanced":
            sizes = sorted(range(len(clusters)), key=lambda k: len(clusters[k]))
            a, b = sizes[0], sizes[1]
        else:
            a, b = rng.sample(range(len(clusters)), 2)
        ca, cb = clusters[a], clusters[b]
        topol.append((sorted(ca), sorted(cb)))
        merged = sorted(ca + cb)
        clusters = [c for k, c in enumerate(clusters) if k not in (a, b)]
        clusters.append(merged)
        lens.append(None)  # filled by the caller
    return topol


def case(rng, n, shape, lengths):
    topol = merge_tree(rng, n, shape)
    f = [str(n)]
    for k, (m0, m1) in enumerate(topol):
        l0, l1 = lengths(rng, k)
        f += [
            ",".join(str(x) for x in m0),
            ",".join(str(x) for x in m1),
            "%.17g" % l0,
            "%.17g" % l1,
        ]
    return "\t".join(f)


def main():
    small = "--small" in sys.argv
    rng = random.Random(20260905)
    ordinary = lambda r, k: (r.uniform(0.001, 0.4), r.uniform(0.001, 0.4))
    with_zeros = lambda r, k: (
        0.0 if r.random() < 0.35 else r.uniform(0.001, 0.4),
        0.0 if r.random() < 0.35 else r.uniform(0.001, 0.4),
    )
    negatives = lambda r, k: (r.uniform(-0.5, 0.3), r.uniform(-0.5, 0.3))
    tiny = lambda r, k: (r.uniform(0.0, 0.002), r.uniform(0.0, 0.002))

    out = []
    # treeCnv indexes stopol[locnseq-3], so three sequences is the floor.
    sizes = (3, 4, 5, 6, 8, 12) if small else (3, 4, 5, 6, 8, 12, 20, 36)
    for n in sizes:
        for shape in ("random", "ladder", "balanced"):
            for lengths in (ordinary, with_zeros, negatives, tiny):
                out.append(case(rng, n, shape, lengths))
    for line in out:
        print(line)


if __name__ == "__main__":
    sys.exit(main())
