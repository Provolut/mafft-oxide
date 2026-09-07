import random, sys
n, shape = int(sys.argv[1]), sys.argv[2]
rng = random.Random(1234 if shape == 'ladder' else 5678)
out = []
if shape == 'ladder':
    for j in range(2, n + 1):
        out.append('%d %d %.5f %.5f' % (1, j, rng.random() * 0.3, rng.random() * 0.3))
else:
    slots = list(range(1, n + 1))
    while len(slots) > 1:
        nxt = []
        for i in range(0, len(slots) - 1, 2):
            a, b = slots[i], slots[i + 1]
            out.append('%d %d %.5f %.5f'
                       % (min(a, b), max(a, b), rng.random() * 0.2, rng.random() * 0.2))
            nxt.append(min(a, b))
        if len(slots) % 2:
            nxt.append(slots[-1])
        slots = nxt
print('\n'.join(out))
