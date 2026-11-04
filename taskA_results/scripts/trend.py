#!/usr/bin/env python3
"""Thermalisation test that does not rely on binning errors: for each (N, z), split every seed's
measurement into quarters; trend = mean(Q_edge, last quarter) - mean(first quarter), with the error
taken from the seed-to-seed scatter of the per-seed trends.  Also reports seed-to-seed scatter of the
last-half means.   usage: trend.py OUTDIR [N]"""
import sys, glob, re, os, numpy as np
out = sys.argv[1]; Nsel = sys.argv[2] if len(sys.argv) > 2 else r'\d+'
g = {}
for f in glob.glob(os.path.join(out, 'N*_s*_r*.txt')):
    m = re.search(rf'N({Nsel})_s(\d+)_r\d+\.txt$', f)
    if not m: continue
    z = float(open(f).readline().split()[-1]); Q = np.array([float(l.split()[2]) for l in open(f) if l[0] == 'B'])
    if len(Q) >= 8: g.setdefault((int(m[1]), z), []).append(Q[:len(Q) // 4 * 4].reshape(4, -1).mean(1))
print(f"{'N':>3} {'z':>9} {'Q_un(2nd half)':>14} {'err':>6} {'q1->q4 trend':>12} {'err':>6} {'t':>6}  quarters(Q_un, seed avg)")
for (N, z), qs in sorted(g.items()):
    q = np.array(qs); k = len(q); tr = q[:, 3] - q[:, 0]; h2 = q[:, 2:].mean(1)
    et = tr.std(ddof=1) / np.sqrt(k) if k > 1 else np.nan; eh = h2.std(ddof=1) / np.sqrt(k) if k > 1 else np.nan
    flag = '  <-- trend' if abs(tr.mean()) > 2 * et else ''
    print(f"{N:3d} {z:9.3g} {N/2-h2.mean():14.3f} {eh:6.3f} {-tr.mean():12.3f} {et:6.3f} {-tr.mean()/et:6.1f}  {np.round(N/2-q.mean(0),3)}{flag}")
