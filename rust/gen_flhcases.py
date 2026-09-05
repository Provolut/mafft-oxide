#!/usr/bin/env python3
"""Emit the argument vectors `verify_flh.sh` feeds to both `Falign_localhom`s.

One case per line, tab-separated, in the order `gen_flhharness.c` and
`examples/flhref.rs` read them:

    n1 <rows> <effs> <orinums> n2 <rows> <effs> <orinums>
    outgap legacygapcost fastathreshold fftWinSize fftThreshold
    nlh <nlh x (i j start1 end1 start2 end2 opt importance)>

Both groups are rows of one alignment, so every row -- in both groups -- has
the same length; that is what `tditeration` passes and the C never sees
anything else here. The default `fftWinSize` / `fftThreshold` / `fastathreshold`
are the values a traced `mafft --auto` run reports (20 / 50 / 2.7).
"""
import random
import sys

AA = "ACDEFGHIKLMNPQRSTVWY"


def rows(rng, n, length, gap_rate, sim):
    """`n` rows of one alignment: a common backbone, mutated, with gap columns."""
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


def chains(rng, n1, n2, length, density):
    """Constraint segments, in the shape `msshrinklocalhom_fast_half` leaves."""
    out = []
    for i in range(n1):
        for j in range(n2):
            if rng.random() > density:
                continue
            for _ in range(rng.randint(1, 3)):
                s1 = rng.randrange(0, max(1, length - 10))
                ln = rng.randint(3, max(4, min(40, length - s1 - 1)))
                s2 = min(max(0, s1 + rng.randint(-4, 4)), length - ln - 1)
                out.append(
                    (i, j, s1, s1 + ln, s2, s2 + ln,
                     rng.uniform(10.0, 200.0), rng.uniform(0.05, 1.5))
                )
    return out


def case(rng, n1, n2, length, gap_rate, sim, density, outgap, legacy,
         fthr=2.7, fftwin=20, fftthr=50, swap=False):
    g1 = rows(rng, n1, length, gap_rate, sim)
    g2 = rows(rng, n2, length, gap_rate, sim)
    e1 = [rng.uniform(0.05, 1.0) for _ in range(n1)]
    e2 = [rng.uniform(0.05, 1.0) for _ in range(n2)]
    # `fillimp` reads a constraint's two halves the other way round when the
    # first group's member has the higher original index, so both orders are
    # generated.
    if swap:
        o1 = list(range(n2, n2 + n1))
        o2 = list(range(n2))
    else:
        o1 = list(range(n1))
        o2 = list(range(n1, n1 + n2))
    lh = chains(rng, n1, n2, length, density)
    f = ["%d" % n1] + g1 + ["%.17g" % v for v in e1] + ["%d" % v for v in o1]
    f += ["%d" % n2] + g2 + ["%.17g" % v for v in e2] + ["%d" % v for v in o2]
    f += ["1" if outgap else "0", "1" if legacy else "0", "%.17g" % fthr]
    f += ["%d" % fftwin, "%d" % fftthr, "%d" % len(lh)]
    for (i, j, s1, e1_, s2, e2_, opt, imp) in lh:
        f += ["%d" % i, "%d" % j, "%d" % s1, "%d" % e1_, "%d" % s2, "%d" % e2_,
              "%.17g" % opt, "%.17g" % imp]
    return "\t".join(f)


def main():
    rng = random.Random(20260905)
    out = []
    # The shapes tditeration produces: one member against the rest, and even
    # splits, over a range of alignment widths.
    for (n1, n2) in [(1, 1), (1, 3), (1, 8), (2, 2), (3, 5), (5, 5), (2, 9)]:
        for length in (40, 120, 240):
            for outgap in (True, False):
                out.append(case(rng, n1, n2, length, 0.06, 0.85, 0.7, outgap, False))
                out.append(case(rng, n1, n2, length, 0.20, 0.60, 0.4, outgap, False,
                                swap=True))
    # Legacy gap cost, and the anchoring parameters off their defaults: a
    # smaller window and a lower threshold both make `alignable_reagion` cut
    # more blocks, which is the branch this function exists to drive.
    for length in (60, 180):
        out.append(case(rng, 2, 4, length, 0.10, 0.75, 0.6, True, True))
        out.append(case(rng, 2, 4, length, 0.10, 0.75, 0.6, True, False,
                        fftwin=8, fftthr=20))
        out.append(case(rng, 3, 3, length, 0.10, 0.75, 0.6, False, False,
                        fftwin=4, fftthr=10))
        out.append(case(rng, 1, 5, length, 0.10, 0.75, 0.6, True, False,
                        fftwin=1000, fftthr=50))
    # Degenerate corners: no constraints at all, a single column, no gaps,
    # every row identical, and a constraint chain on every pair.
    out.append(case(rng, 2, 3, 80, 0.05, 0.8, 0.0, True, False))
    out.append(case(rng, 1, 1, 1, 0.0, 1.0, 0.0, True, False))
    out.append(case(rng, 2, 2, 50, 0.0, 1.0, 1.0, True, False))
    out.append(case(rng, 3, 4, 90, 0.35, 0.95, 1.0, False, False))
    out.append(case(rng, 4, 4, 150, 0.15, 0.7, 1.0, True, False, fthr=0.5))
    out.append(case(rng, 4, 4, 150, 0.15, 0.7, 1.0, True, False, fthr=10.0))
    for line in out:
        print(line)


if __name__ == "__main__":
    sys.exit(main())
