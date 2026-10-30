#!/usr/bin/env python3
"""Equilibrium estimates from long runs: discard the first DISCARD sweeps, pool all seeds/dirs.
Per required z: mean Q_un, error (seed-to-seed scatter), integrated autocorrelation time of the
block series (from binning: tau = BL/2 * (err_binned/err_naive)^2), and the sweeps per seed needed
for a target error with 4 seeds.
usage: diag_analyze.py N DISCARD_SWEEPS TARGET_ERR DIR1 [DIR2 ...]"""
import sys, glob, numpy as np
N, disc, target = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]); dirs = sys.argv[4:]
REQ = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1]
def per_block(d):
    hdr = open(glob.glob(f'{d}/N{N}_s*.ladder')[0]).readline()
    return None
print(f"{'z':>8} {'Q_un':>8} {'err':>6} {'tau(sweeps)':>12} {'std(Q)':>7} {'sweeps/seed for err':>20}  per-seed means")
for z in REQ:
    series = []
    for d in dirs:
        for f in sorted(glob.glob(f'{d}/N{N}_s*_r*.txt')):
            if f.endswith('_hist.txt') or abs(float(open(f).readline().split()[-1]) / z - 1) > 1e-6: continue
            Q = N / 2 - np.array([float(l.split()[2]) for l in open(f) if l[0] == 'B'])
            series.append(Q)
    if not series: continue
    BL = 32000000 / 800  # sweeps per block in the diagnostic runs
    k0 = int(disc / BL); S = [s[k0:] for s in series]
    means = np.array([s.mean() for s in S]); err = means.std(ddof=1) / np.sqrt(len(means))
    taus = []
    for s in S:
        naive = s.std(ddof=1) / np.sqrt(len(s)); best = naive; b = 2
        while len(s) // b >= 8:
            y = s[: len(s) // b * b].reshape(-1, b).mean(1); best = max(best, y.std(ddof=1) / np.sqrt(len(y))); b *= 2
        taus.append(BL / 2 * (best / naive) ** 2)
    tau = np.mean(taus); L = len(S[0]) * BL
    var_mean = err ** 2 * len(means)  # variance of a single-seed mean over length L
    need = L * var_mean / (4 * target ** 2)  # sweeps per seed so that 4 seeds give `target`
    print(f"{z:8.0e} {means.mean():8.3f} {err:6.3f} {tau:12,.0f} {np.mean([s.std() for s in S]):7.3f} {need:20,.0f}  {np.round(means, 3)}")
