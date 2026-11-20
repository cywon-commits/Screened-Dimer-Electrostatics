"""Fig. 2: vacancy isotherms of the two uniform geometries (exact transfer matrix vs closed forms)."""
import sys, numpy as np, matplotlib.pyplot as plt
sys.path.insert(0, "../isotherms")
from style import *
import tm_diag_b as T
from loz_cap import lnZ
from ckp import L
fig, ax = plt.subplots(1, 2, figsize=(TEXTW, 2.5), gridspec_kw=dict(wspace=0.28))
ts = np.linspace(0.0, 1, 200); ax[0].plot(ts, (2 / np.pi) * np.arctan(ts), "k-", label=r"$(2/\pi)\arctan z^{2/W}$")
T.RIN = 0.0
for W, m, c in ((6, "v", SKY), (10, "s", BLUE), (14, "o", VERM)):
    tt = [t for t in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0) if t ** (W / 2) > 1e-3]   # smaller z: the sparse eigensolver does not converge
    ax[0].plot(tt, [T.densities(W, t ** (W / 2))[2][0] for t in tt], m, mfc="none", color=c, label=f"$W={W}$")
ax[0].set_xlabel(r"$z^{2/W}$"); ax[0].set_ylabel(r"edge vacancy density $\rho$"); ax[0].set_xlim(0, 1.03); ax[0].set_ylim(0, 0.55); ax[0].legend(loc="upper left")
s = lambda r: (2 * L(np.pi * r / 2) + L(np.pi * (1 - r))) / np.pi
xs = np.linspace(0, 1, 300); r = (2 / np.pi) * np.arccos(2 ** (xs - 1)); ax[1].plot(xs, s(r) - r * xs * np.log(2), "k-", label="closed form")
M = 12; xp = (0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95); f = {N: np.array([lnZ(M, N, 2.0 ** (-x * N / 2)) / (M * N) for x in xp]) for N in (32, 64)}
ax[1].plot(xp, f[32], "s", color=SKY, mfc="none", label="$N=32$"); ax[1].plot(xp, f[64], "o", color=BLUE, mfc="none", label="$N=64$")
ax[1].plot(xp, 2 * f[64] - f[32], "x", color=VERM, label=r"$N\to\infty$ (linear in $1/N$)")
ax[1].set_xlabel(r"$x=N^*/N$"); ax[1].set_ylabel("free energy per site"); ax[1].set_xlim(0, 1.03); ax[1].set_ylim(bottom=0); ax[1].legend(loc="upper right")
panel(ax[0], "(a)", x=-0.13); panel(ax[1], "(b)", x=-0.13)
save(fig, "fig2_isotherms")
