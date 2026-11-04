#!/usr/bin/env python3
"""Exact Q_edge / Q_un from the transfer computation (aztec.py), with the same edge weights as aztec_pt.c.
usage: exact_qun.py N z1 z2 ...   -> lines 'N z Q_un s n_mono'"""
import sys, numpy as np
import aztec as A

def qun(N, z):
    lnZ, mono, left, right = A.Aztec(N, z).run()
    Q = n = 0.0
    for k in range(2 * N + 1):
        for i, pm in enumerate(mono[k]):
            x, y = A.cell_coords(N, k, i)
            kk, b = x + y - (N - 1), x - y; uu = kk - N
            w = (1 if uu < 0 else 0) if abs(uu) > abs(b) else (0 if abs(b) > abs(uu) else (0.5 if uu < 0 else 0))
            Q += (1 if kk % 2 == 0 else -1) * w * pm; n += pm
    return N / 2 - Q, 2 * Q / N, n

if __name__ == '__main__':
    N = int(sys.argv[1])
    for z in map(float, sys.argv[2:]):
        q, s, n = qun(N, z); print(f"{N} {z:g} {q:.6f} {s:.6f} {n:.6f}", flush=True)
