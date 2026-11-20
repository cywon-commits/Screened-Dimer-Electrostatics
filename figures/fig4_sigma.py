"""Fig. 4: uncompensated fraction sigma_inf(x): exact solution vs boundary-Pfaffian enumeration."""
import sys, pickle, numpy as np, matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.optimize import brentq
sys.path.insert(0, "../src")
from style import *
from elliptic_sol import Sol
d = pickle.load(open("../data/bx_exact.pkl", "rb"))
rs = np.concatenate([np.geomspace(2e-3, 0.3, 25), np.linspace(0.32, 0.985, 45)]); S = [Sol(r) for r in rs]
X = np.array([s.x() for s in S]); Sg = np.array([s.sigma() for s in S]); o = np.argsort(X); X, Sg = X[o], Sg[o]
fig, ax = plt.subplots(1, 2, figsize=(TEXTW, 2.5), gridspec_kw=dict(wspace=0.28))
A = 8 * np.log(2) / np.pi ** 2; B = 4 * np.pi ** 1.5 * np.e / (gamma(0.25) ** 2 * np.log(2)); xx = np.linspace(0.005, 0.75, 200)
ax[0].plot(np.r_[0, X, 1], np.r_[0, Sg, 1], "k-", label="exact solution"); ax[0].plot(xx, A * xx * np.log(B / xx), "--", color=GREY, label=r"small-$x$ law")
xs = sorted({x for (_, x) in d})
for N, m, c in ((64, "s", SKY), (192, "o", VERM)):
    ax[0].plot(xs, [d[(N, x)] for x in xs], m, mfc="none", color=c, label=f"$N={N}$")
ax[0].set_xlabel(r"$x=N^*/N$"); ax[0].set_ylabel(r"uncompensated fraction $\sigma$"); ax[0].set_xlim(0, 1); ax[0].set_ylim(0, 1.02); ax[0].legend(loc="lower right")
for x in xs:
    Ns = np.array(sorted(N for (N, y) in d if y == x)); se = Sol(brentq(lambda r: Sol(r).x() - x, 1e-3, 0.995, xtol=1e-13)).sigma()
    ax[1].plot(1 / Ns, [se - d[(N, x)] for N in Ns], "-", marker=MARK[x], color=xcolor(x), mfc="none", lw=0.6, label=f"${x}$")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel(r"$1/N$"); ax[1].set_ylabel(r"$\sigma_\infty-\sigma_N$"); ax[1].legend(title="$x$", ncol=2, loc="lower right", title_fontsize=7)
panel(ax[0], "(a)", x=-0.13); panel(ax[1], "(b)", x=-0.13)
save(fig, "fig4_sigma")
for x in xs: print(x, [(N, d[(N, x)]) for N in sorted(N for (N, y) in d if y == x)])
