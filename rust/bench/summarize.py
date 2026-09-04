#!/usr/bin/env python3
"""Aggregate `profile.tsv` into the numbers GOAL.md's definition of done uses.

    ./summarize.py [profile.tsv]

Two aggregates, because they disagree and only one of them is honest.

The *unweighted mean of ratios* treats a 50 ms case and a 1.5 s case as equally
important. On this database that flatters the port badly: the three small
inputs run faster than the C (startup dominates them, and the C's startup is
slower), which pulls the mean down and hides that the two heaviest cases are
the two worst.

The *work-weighted ratio* -- the sum of the port's cost over the sum of the
C's -- is what a caller actually experiences across the workload, and is the
figure GOAL.md means by "on average". It is reported first for that reason.

With a second argument it also compares against a frozen baseline, which is how
GOAL.md's "<= 90% of today" is checked rather than eyeballed:

    ./summarize.py profile.tsv baseline-2026-09-04.tsv
"""
import sys


def totals(path):
    """(wall clock ms summed, mean peak RSS KB) over the rows that measured."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(header)}
        ms = rss = n = 0.0
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header) or f[idx["note"]] != "ok":
                continue
            ms += float(f[idx["r_ms"]])
            rss += float(f[idx["r_rss_kb"]])
            n += 1
    return (ms, rss / n if n else 0.0)


def compare_to_baseline(path, base):
    ms, rss = totals(path)
    bms, brss = totals(base)
    print("\nagainst %s (GOAL.md target: <= 90%% of either)" % base)
    print("-" * 70)
    hit = False
    for label, now, was, unit in (
        ("wall clock", ms, bms, "ms"),
        ("mean peak RSS", rss, brss, "KB"),
    ):
        pct = 100.0 * now / was if was else 0.0
        hit = hit or pct <= 90.0
        print(
            "%-14s %8.0f %s   was %8.0f %s   %6.1f%% of today  %s"
            % (label, now, unit, was, unit, pct, "MET" if pct <= 90.0 else "")
        )
    print("\n%s" % ("TARGET MET" if hit else "target not met"))
    return 0 if hit else 1


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "profile.tsv"
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header) or f[idx["note"]] != "ok":
                continue
            rows.append(f)

    if not rows:
        print("no usable rows in %s" % path, file=sys.stderr)
        return 1

    # Instructions and milliseconds accumulate across a workload, so summing
    # them is meaningful and the weighted ratio is the figure that counts.
    # Peak RSS does not: two runs peaking at 20 MB do not need 40 MB, they need
    # 20 MB twice. For memory the per-run mean and the worst case are what a
    # caller provisions for, so the weighted column is suppressed rather than
    # printed and quietly misread.
    metrics = [
        ("instructions", "c_ir", "r_ir", True),
        ("wall clock", "c_ms", "r_ms", True),
        ("peak RSS", "c_rss_kb", "r_rss_kb", False),
    ]

    print("%d cases\n" % len(rows))
    print("%-14s %14s %14s %10s %10s" % ("metric", "C total", "port total", "weighted", "mean"))
    print("-" * 66)
    for label, ckey, rkey, summable in metrics:
        ctot = sum(float(r[idx[ckey]]) for r in rows)
        rtot = sum(float(r[idx[rkey]]) for r in rows)
        weighted = rtot / ctot if ctot else 0.0
        ratios = [
            float(r[idx[rkey]]) / float(r[idx[ckey]])
            for r in rows
            if float(r[idx[ckey]])
        ]
        mean = sum(ratios) / len(ratios) if ratios else 0.0
        if summable:
            print(
                "%-14s %14.0f %14.0f %9.3fx %9.3fx"
                % (label, ctot, rtot, weighted, mean)
            )
        else:
            print(
                "%-14s %14s %14s %10s %9.3fx   (worst %.3fx)"
                % (label, "-", "-", "n/a", mean, max(ratios) if ratios else 0.0)
            )

    print("\nworst cases by weighted contribution:")
    for label, ckey, rkey, _ in metrics:
        worst = sorted(
            rows,
            key=lambda r: (float(r[idx[rkey]]) - float(r[idx[ckey]])),
            reverse=True,
        )[:3]
        names = ", ".join(
            "%s (%.2fx)" % (r[idx["tag"]], float(r[idx[rkey]]) / float(r[idx[ckey]]))
            for r in worst
            if float(r[idx[ckey]])
        )
        print("  %-13s %s" % (label + ":", names))

    if len(sys.argv) > 2:
        return compare_to_baseline(path, sys.argv[2])
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        # Piped through `head`; not an error worth a traceback.
        sys.exit(0)
