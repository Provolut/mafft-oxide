#!/usr/bin/env python3
"""Generate argument lines for the `Gaplen` differential harness.

    gen_glcases.py [seed] > cases.txt

Each line is a tab-separated argument vector for both `gen_glharness.c` (built
against `core/`) and `cargo run --example glref`: two groups of aligned rows
with their weights. `verify_gaplen.sh` runs both and diffs them.
"""
import random
import sys

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 4242)
AA = "ACDEFGHIKLMNPQRSTVWY"


def rows(n, width, gapprob, maxgap):
    out = []
    for _ in range(n):
        r = []
        while len(r) < width:
            if random.random() < gapprob:
                r += ['-'] * random.randint(1, maxgap)
            else:
                r += [random.choice(AA)]
        out.append(''.join(r[:width]))
    return out


cases = []
for _ in range(200):
    n1, n2 = random.randint(1, 8), random.randint(1, 8)
    g = random.choice([0.0, 0.05, 0.2, 0.35, 0.5, 0.75, 0.95])
    mg = random.choice([1, 2, 5, 12])
    r1 = rows(n1, random.randint(1, 40), g, mg)
    r2 = rows(n2, random.randint(1, 40), g, mg)
    e1 = [round(random.random(), 6) + 0.001 for _ in range(n1)]
    e2 = [round(random.random(), 6) + 0.001 for _ in range(n2)]
    cases.append([str(n1)] + r1 + [str(x) for x in e1] + [str(n2)] + r2 + [str(x) for x in e2])

# All-gap rows, a single column, and a group whose rows share no gap pattern:
# the cases where a column list is empty or every entry has relend != 0.
cases.append(['1', '-', '1.0', '1', '-', '1.0'])
cases.append(['2', '----', '----', '0.5', '0.5', '1', 'ACGT', '1.0'])
cases.append(['3', 'A---B', 'AB---', '---AB', '0.3', '0.3', '0.4',
              '2', '-----', 'ABCDE', '0.5', '0.5'])

for c in cases:
    print('\t'.join(c))
