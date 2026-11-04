#!/usr/bin/env python3
"""Cold (brick-wall) vs hot (all-monomer) start: Q_un per time window at the required z.
Equilibrium is reached where the two starts agree and stop moving.
usage: diag_compare.py COLDDIR HOTDIR N [WINDOW_BLOCKS]"""
import sys, glob, numpy as np
cold, hot, N = sys.argv[1], sys.argv[2], int(sys.argv[3]); W = int(sys.argv[4]) if len(sys.argv) > 4 else 100
REQ = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1]
def series(d, z):
    out = []
    for f in sorted(glob.glob(f'{d}/N{N}_s*_r*.txt')):
        if f.endswith('_hist.txt') or abs(float(open(f).readline().split()[-1]) / z - 1) > 1e-6: continue
        Q = np.array([float(l.split()[2]) for l in open(f) if l[0] == 'B']); k = len(Q) // W
        out.append(N / 2 - Q[:k * W].reshape(k, W).mean(1))
    k = min(map(len, out)) if out else 0
    return np.array([o[:k] for o in out])
per = None
for z in REQ:
    c, h = series(cold, z), series(hot, z)
    k = min(c.shape[1] if c.size else 0, h.shape[1] if h.size else 0)
    if not k: continue
    print(f"z={z:g}  window | cold seeds | hot seeds")
    for j in range(k):
        print(f"   {j:3d}  {' '.join(f'{x:7.3f}' for x in c[:, j])} | {' '.join(f'{x:7.3f}' for x in h[:, j])}")
