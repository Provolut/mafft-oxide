#!/usr/bin/env python3
"""Emit the refinement runs `verify_tdi.sh` feeds to both `TreeDependentIteration`s.

One case per line, tab-separated, in the order `gen_tdiharness.c` and
`examples/tdiref.rs` read them:

    nseq niter cut constraint fastathreshold
    <nseq aligned rows>
    nseq-1 x ( members0 members1 len0 len1 )
    nseq-1 x ( skip0 skip1 )
    nlh  then nlh x ( i j start1 end1 start2 end2 opt )

The rows are an *alignment* -- equal length, with gaps -- because that is what
the progressive stage hands the refinement loop. Constraints are given by `opt`
only: both sides run `calcimportance_half` to turn those into the importances
`fillimp` reads, so that step is under test too.
"""
import random
import sys

AA = "ACDEFGHIKLMNPQRSTVWY"


def alignment(rng, n, length, gap_rate, sim):
    base = [rng.choice(AA) for _ in range(length)]
    out = []
    for _ in range(n):
        r = list(base)
        for p in range(length):
            if rng.random() > sim:
                r[p] = rng.choice(AA)
            if rng.random() < gap_rate:
                r[p] = "-"
        out.append("".join(r))
    return out


def merge_tree(rng, n, shape):
    clusters = [[i] for i in range(n)]
    topol = []
    for _ in range(n - 1):
        if shape == "ladder":
            a, b = 0, 1
        else:
            a, b = rng.sample(range(len(clusters)), 2)
        ca, cb = clusters[a], clusters[b]
        topol.append((sorted(ca), sorted(cb)))
        merged = sorted(ca + cb)
        clusters = [c for k, c in enumerate(clusters) if k not in (a, b)]
        clusters.append(merged)
    return topol


def constraints(rng, n, length, density):
    """Segments in the shape a `hat3` table holds them: pair, then a diagonal."""
    out = []
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() > density:
                continue
            for _ in range(rng.randint(1, 2)):
                s1 = rng.randrange(0, max(1, length - 12))
                ln = rng.randint(4, max(5, min(30, length - s1 - 1)))
                s2 = min(max(0, s1 + rng.randint(-3, 3)), length - ln - 1)
                out.append((i, j, s1, s1 + ln, s2, s2 + ln, rng.uniform(10.0, 150.0)))
    return out


def case(rng, n, length, niter, cut, constraint, shape="random", skip_rate=0.0,
         gap_rate=0.08, sim=0.8, density=0.5, fthr=2.7):
    rows = alignment(rng, n, length, gap_rate, sim)
    topol = merge_tree(rng, n, shape)
    f = [str(n), str(niter), "%.17g" % cut, str(constraint), "%.17g" % fthr]
    f += rows
    for (m0, m1) in topol:
        f += [",".join(map(str, m0)), ",".join(map(str, m1)),
              "%.17g" % rng.uniform(0.005, 0.4), "%.17g" % rng.uniform(0.005, 0.4)]
    for _ in range(n - 1):
        f += ["1" if rng.random() < skip_rate else "0",
              "1" if rng.random() < skip_rate else "0"]
    lh = constraints(rng, n, length, density) if constraint else []
    f += [str(len(lh))]
    for (i, j, s1, e1, s2, e2, opt) in lh:
        f += [str(i), str(j), str(s1), str(e1), str(s2), str(e2), "%.17g" % opt]
    return "\t".join(f)


def main():
    small = "--small" in sys.argv
    rng = random.Random(70707)
    out = []
    sizes = [(4, 40), (5, 60), (6, 50), (8, 70)] if small else [
        (4, 40), (5, 60), (6, 50), (8, 70), (10, 90), (14, 80), (20, 60)
    ]
    for (n, length) in sizes:
        for constraint in (2, 0):
            # One round and several, since the walk direction alternates and
            # the oscillation check only arms from the second round on.
            for niter in (1, 2, 5):
                out.append(case(rng, n, length, niter, 0.0, constraint))
            # A branch the caller has excluded, and a non-zero acceptance
            # margin, which is the only thing `cut` changes.
            out.append(case(rng, n, length, 3, 0.0, constraint, skip_rate=0.3))
            out.append(case(rng, n, length, 3, 5.0, constraint))
            # Ladder trees, dense gaps, and near-identical rows: the last one
            # drives the `identical` branch and, through it, convergence.
            out.append(case(rng, n, length, 3, 0.0, constraint, shape="ladder"))
            out.append(case(rng, n, length, 3, 0.0, constraint, gap_rate=0.18, sim=0.68))
            out.append(case(rng, n, length, 4, 0.0, constraint, sim=0.995, gap_rate=0.01))
            if constraint:
                out.append(case(rng, n, length, 3, 0.0, constraint, density=0.0))
                out.append(case(rng, n, length, 3, 0.0, constraint, density=1.0))
    for line in out:
        print(line)


if __name__ == "__main__":
    sys.exit(main())
