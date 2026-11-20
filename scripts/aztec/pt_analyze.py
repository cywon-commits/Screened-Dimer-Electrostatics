#!/usr/bin/env python3
"""Analysis of aztec_pt / aztec_mc output.

usage:
  python3 pt_analyze.py OUTDIR            # OUTDIR holds files  N<N>_s<seed>_r<k>.txt  and  ..._hist.txt
Produces summary.csv with one row per (N, z):
  N, z, Q_edge, err, Q_un=N/2-Q_edge, s=2Q_edge/N, Psi, n_mono, drift_sigma, seeds
Errors: blocks from every seed pooled; error = std of per-seed means / sqrt(#seeds) when >=3 seeds,
otherwise binning error of the pooled blocks (largest of bin sizes 1,2,4,...).
drift_sigma = (mean of first half of blocks - second half) / combined error  -> |drift_sigma| > 3 means NOT thermalised.
Psi: per-cell phi = max over {monomer, 4 bond directions} of the state frequencies;
     Psi = < (phi - phi_b)/(1 - phi_b) >, phi_b = max(n_bulk, (1 - n_bulk)/4),
     n_bulk from the table below (bulk monomer density of the infinite plane), interpolated in log z.
"""
import sys, glob, re, os, numpy as np

# bulk monomer density per site n(z) of the infinite square lattice (strip transfer matrix, W<=22)
NB_Z = np.array([1e-6, 1e-3, 3e-3, 1e-2, 3e-2, 5e-2, 7e-2, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0])
NB_N = np.array([0.0, 2.0e-4, 7.0e-4, 2.67e-3, 9.18e-3, 1.634e-2, 2.375e-2, 3.543e-2, 5.522e-2, 7.529e-2, 0.11534, 0.19272, 0.36188])
def n_bulk(z):
    if z <= NB_Z[0]: return 0.0
    # n ~ z^(8/7) at small z
    if z < 1e-3: return NB_N[1] * (z / 1e-3) ** (8 / 7)
    return float(np.interp(np.log(z), np.log(NB_Z[1:]), NB_N[1:]))

def blocks(fn):
    z = None; Q = []; n = []
    for l in open(fn):
        if l.startswith('#'):
            m = re.search(r'z (\S+)', l); z = float(m[1]) if m else z
        elif l.startswith('B'):
            t = l.split(); Q.append(float(t[2])); n.append(float(t[4]))
    return z, np.array(Q), np.array(n)

def binning_err(x):
    best = x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan
    b = 2
    while len(x) // b >= 4:
        y = x[: len(x) // b * b].reshape(-1, b).mean(1)
        best = max(best, y.std(ddof=1) / np.sqrt(len(y))); b *= 2
    return best

def psi(histfile, z):
    a = np.loadtxt(histfile, comments='#'); h = a[:, 2:]
    P = h / np.maximum(h.sum(1, keepdims=True), 1); phi = P.max(1)
    nb = n_bulk(z); pb = max(nb, (1 - nb) / 4)
    return float(np.mean((phi - pb) / (1 - pb)))

def main(outdir):
    files = glob.glob(os.path.join(outdir, 'N*_s*_r*.txt'))
    groups = {}
    for f in files:
        if f.endswith('_hist.txt'): continue
        m = re.search(r'N(\d+)_s(\d+)_r(\d+)\.txt$', f)
        N, seed, r = int(m[1]), int(m[2]), int(m[3])
        z, Q, n = blocks(f)
        groups.setdefault((N, z), []).append((seed, Q, n, f[:-4] + '_hist.txt'))
    rows = []
    for (N, z), runs in sorted(groups.items()):
        means = [Q.mean() for _, Q, _, _ in runs if len(Q)]
        allQ = np.concatenate([Q for _, Q, _, _ in runs])
        if len(allQ) < 4: continue
        err = np.std(means, ddof=1) / np.sqrt(len(means)) if len(means) >= 3 else binning_err(allQ)
        h = len(allQ) // 2
        e1, e2 = binning_err(allQ[:h]), binning_err(allQ[h:])
        drift = (allQ[:h].mean() - allQ[h:].mean()) / np.hypot(e1, e2) if e1 > 0 else np.nan
        ps = np.mean([psi(hf, z) for *_, hf in runs if os.path.exists(hf)])
        nm = np.mean(np.concatenate([n for _, _, n, _ in runs]))
        Qm = allQ.mean()
        rows.append((N, z, Qm, err, N / 2 - Qm, 2 * Qm / N, ps, nm, drift, len(runs)))
    with open(os.path.join(outdir, 'summary.csv'), 'w') as f:
        f.write('N,z,Q_edge,err,Q_un,s,Psi,n_mono,drift_sigma,seeds\n')
        for r in rows: f.write(','.join(f'{v:.6g}' for v in r) + '\n')
    print(f"{'N':>4} {'z':>9} {'Q_un':>9} {'err':>7} {'s':>7} {'Psi':>7} {'drift':>6} seeds")
    for r in rows:
        flag = '  <-- not thermalised' if abs(r[8]) > 3 else ''
        print(f"{r[0]:4d} {r[1]:9.2e} {r[4]:9.3f} {r[3]:7.3f} {r[5]:7.4f} {r[6]:7.4f} {r[8]:6.2f} {r[9]}{flag}")

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
