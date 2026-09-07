#!/usr/bin/env python3
"""Split an alignment into an aligned head and `nadd` sequences to add back.

    gen_add.py <aligned.fa> <nadd> add|profile

`add` ungaps the tail (--add); `profile` leaves it aligned (--addprofile).
"""
import sys

path, nadd, mode = sys.argv[1], int(sys.argv[2]), sys.argv[3]
recs, name, seq = [], None, []
for line in open(path):
    line = line.rstrip('\n')
    if line.startswith('>'):
        if name is not None:
            recs.append((name, ''.join(seq)))
        name, seq = line, []
    else:
        seq.append(line)
recs.append((name, ''.join(seq)))

head, tail = recs[:len(recs) - nadd], recs[len(recs) - nadd:]
out = []
for n, s in head:
    out.append('%s\n%s' % (n, s))
for n, s in tail:
    out.append('%s\n%s' % (n, s if mode == 'profile' else s.replace('-', '')))
print('\n'.join(out))
