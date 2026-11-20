"""Fig. 6: frozen-corner boundaries and frozen length from the exact solution."""
import sys, numpy as np, matplotlib.pyplot as plt
from scipy.optimize import brentq
sys.path.insert(0, "../corner_islands")
from style import *
from exact_sol import solve
from exact_more_lib import arctic
from elliptic_sol import Sol
fig = plt.figure(figsize=(TEXTW, 2.7))
ax = fig.add_axes([0.02, 0.03, 0.40, 0.90]); fig.text(0.005, 0.965, "(a)", fontsize=9, va="top")
rot = lambda a, b: (np.asarray(a) - np.asarray(b), np.asarray(a) + np.asarray(b))
for x in (0.9, 0.7, 0.5, 0.3, 0.1):
    u, v = arctic(1 + np.exp(brentq(lambda l: solve(1 + np.exp(l))["x"] - x, -12, 11, xtol=1e-10)))
    U = np.r_[u, v[::-1]]; V = np.r_[v, u[::-1]]
    for k, (a, b) in enumerate(((U, V), (2 - U, V), (U, 2 - V), (2 - U, 2 - V))): ax.plot(*rot(a, b), "-", color=xcolor(x), label=f"${x}$" if k == 0 else None)
th = np.linspace(0, 2 * np.pi, 400); ax.plot(np.sqrt(2) * np.cos(th), 2 + np.sqrt(2) * np.sin(th), ":", color="k", lw=0.8, label="$1$ (arctic circle)")
ax.plot(*rot([0, 2, 2, 0, 0], [0, 0, 2, 2, 0]), "k-", lw=0.8)
for t, (px, py), ha, va in (("S", (0, -0.08), "center", "top"), ("N", (0, 4.08), "center", "bottom"), ("W", (-2.08, 2), "right", "center"), ("E", (2.08, 2), "left", "center")): ax.text(px, py, t, ha=ha, va=va)
ax.legend(title="$x$", loc="center", title_fontsize=7, handlelength=1.4); ax.set_xlim(-2.4, 2.4); ax.set_ylim(-0.4, 4.4); ax.set_aspect("equal"); ax.axis("off")
ax = fig.add_axes([0.56, 0.17, 0.42, 0.76]); fig.text(0.47, 0.965, "(b)", fontsize=9, va="top")
rs = np.concatenate([np.geomspace(2e-3, 0.3, 25), np.linspace(0.32, 0.999, 60)]); S = [Sol(r) for r in rs]
X = np.array([s.x() for s in S]); D = np.array([s.dstar() for s in S]); o = np.argsort(X)
ax.plot(np.r_[0, X[o], 1], np.r_[0, D[o], 1], "k-", label="exact solution"); xx = np.linspace(0, 1, 50); ax.plot(xx, 4 * xx * np.log(2) / np.pi, "--", color=GREY, label=r"isolated corner, $(4\ln 2/\pi)\,x$")
ax.set_xlabel(r"$x=N^*/N$"); ax.set_ylabel(r"frozen length $d^*/(N/2)$"); ax.set_xlim(0, 1); ax.set_ylim(0, 1.02); ax.legend(loc="upper left")
save(fig, "fig6_frozen")
print("check d*:", [(round(x, 3), round(dd, 4)) for x, dd in zip(X[o][::12], D[o][::12])])
