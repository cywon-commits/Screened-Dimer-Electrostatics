#!/usr/bin/env python3
"""Figure for the diagnostic report: cold vs hot start Q_un in 4e6-sweep windows, N=32.
usage: python3 reports/plot_diag.py"""
import glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
C = {'cold': '#2a78d6', 'hot': '#eb6834'}; INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.6, 'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb'})
W = 100; BL = 4e4
fig, axs = plt.subplots(1, 3, figsize=(10, 3.4), sharex=True)
for ax, z in zip(axs, [1e-5, 3e-5, 1e-4]):
    for kind, d in (('cold', 'diag32c'), ('hot', 'diag32h')):
        first = True
        for f in sorted(glob.glob(f'{d}/N32_s*_r*.txt')):
            if f.endswith('_hist.txt') or abs(float(open(f).readline().split()[-1]) / z - 1) > 1e-6: continue
            Q = 16 - np.array([float(l.split()[2]) for l in open(f) if l[0] == 'B']); k = len(Q) // W
            y = Q[:k * W].reshape(k, W).mean(1); t = (np.arange(k) + 0.5) * W * BL / 1e6
            ax.plot(t, y, color=C[kind], lw=2, marker='o', ms=4, label=f'{kind} start' if first else None); first = False
    ax.set_title(f'z = {z:g}', fontsize=10, color=INK, loc='left'); ax.set_xlabel('sweeps (10⁶)')
axs[0].set_ylabel('Q_un (4×10⁶-sweep window)'); axs[0].legend(frameon=False)
fig.suptitle('N=32: cold (brick-wall) vs hot (all-monomer) start, 2 seeds each', x=0.01, ha='left', fontsize=10, color=INK)
fig.tight_layout(); fig.savefig('reports/fig3_diag_cold_hot.png', dpi=150)
