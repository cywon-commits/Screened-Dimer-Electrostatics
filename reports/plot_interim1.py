#!/usr/bin/env python3
"""Figures for the first interim report (N=32).  usage: python3 reports/plot_interim1.py"""
import glob, re, csv, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6, 'figure.facecolor': '#fcfcfb',
                     'axes.facecolor': '#fcfcfb'})

def quarters(out, N, z):
    q = []
    for f in glob.glob(f'{out}/N{N}_s*_r*.txt'):
        if f.endswith('_hist.txt') or abs(float(open(f).readline().split()[-1]) - z) > 1e-12 * max(z, 1): continue
        Q = np.array([float(l.split()[2]) for l in open(f) if l[0] == 'B']); q.append(N / 2 - Q.reshape(4, -1).mean(1))
    return np.array(q)

# Fig 1: thermalisation trend at N=32
fig, ax = plt.subplots(figsize=(6.4, 4))
t = 2 + np.array([1, 3, 5, 7])  # quarter centres, 1e6 sweeps (2e6 thermalisation)
for c, z in zip(C, [3e-5, 1e-4, 3e-4, 1e-3]):
    q = quarters('ana32', 32, z); d = q - q[:, 3:4]; m = d.mean(0); e = d.std(0, ddof=1) / np.sqrt(len(d))
    ax.errorbar(t, m, yerr=e, color=c, lw=2, marker='o', ms=6, capsize=3, label=f'z = {z:g}')
    ax.annotate(f'{z:g}', (t[0], m[0]), xytext=(-8, 0), textcoords='offset points', ha='right', va='center', color=INK, fontsize=9)
ax.axhline(0, color=MUTED, lw=0.8)
ax.set_xlabel('sweeps (10⁶), quarter centres of the measurement'); ax.set_ylabel('Q_un(quarter) − Q_un(last quarter)')
ax.set_xlim(1.5, 9.5); ax.legend(frameon=False, loc='upper right')
ax.set_title('N=32: Q_un still drifting at z ≤ 6×10⁻⁴ (mean of 4 seeds, ±1σ seed scatter)', fontsize=10, color=INK, loc='left')
fig.tight_layout(); fig.savefig('reports/fig1_N32_trend.png', dpi=150)

# Fig 2: Psi(s) collapse, N=10, 20, 32
fig, ax = plt.subplots(figsize=(6.4, 4))
for c, (out, N) in zip(C, [('val', 10), ('val20', 20), ('ana32', 32)]):
    r = [x for x in csv.DictReader(open(f'{out}/summary.csv')) if int(x['N']) == N]
    s = np.array([float(x['s']) for x in r]); p = np.array([float(x['Psi']) for x in r]); z = np.array([float(x['z']) for x in r])
    ok = z >= 8.5e-4 if N == 32 else np.ones_like(z, bool)
    ax.plot(s[ok], p[ok], color=c, lw=2, marker='o', ms=5, label=f'N = {N}')
    if (~ok).any():
        ax.plot(s[~ok], p[~ok], color=c, lw=1, ls=':', marker='o', ms=5, mfc='#fcfcfb', label=f'N = {N}, not thermalised (z < 8.5×10⁻⁴)')
ax.set_xlabel('s = 2 Q_edge / N'); ax.set_ylabel('Ψ'); ax.legend(frameon=False)
ax.set_title('Frozen order vs screened fraction', fontsize=10, color=INK, loc='left')
fig.tight_layout(); fig.savefig('reports/fig2_psi_vs_s.png', dpi=150)
