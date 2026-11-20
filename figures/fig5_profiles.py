"""Fig. 5: vacancy density along an edge: exact solution vs boundary-Pfaffian enumeration."""
import sys, pickle, numpy as np, matplotlib.pyplot as plt
from scipy.optimize import brentq
sys.path.insert(0, "../corner_islands")
from style import *
from exact_sol import solve
sm = pickle.load(open("../data/smallx.pkl", "rb")); pr = pickle.load(open("../data/profiles.pkl", "rb"))
fig, ax = plt.subplots(figsize=(COLW, 2.6))
for x in (0.2, 0.3, 0.5, 0.7):
    s = solve(1 + np.exp(brentq(lambda l: solve(1 + np.exp(l))["x"] - x, -9, 11, xtol=1e-10))); u = s["u"]; rho = (1 - s["h"]) / 2
    ax.plot(np.r_[0, u, (2 - u)[::-1], 2] / 2, np.r_[0, rho, rho[::-1], 0], "-", color=xcolor(x))
    n56 = sm[(56, x)][1] if (56, x) in sm else pr[(56, x)][1][:56]; n40 = sm[(40, x)][1] if (40, x) in sm else pr[(40, x)][1][:40]
    dl = (np.arange(56) + .5) / 56; p40 = np.interp(dl, (np.arange(40) + .5) / 40, n40)
    ax.plot(dl[::2], ((56 * n56 - 40 * p40) / 16)[::2], MARK[x], color=xcolor(x), mfc="none", ms=3, mew=0.6, label=f"${x}$")
ax.set_xlabel(r"position along the edge, $s/N$"); ax.set_ylabel(r"vacancy density $\rho$"); ax.set_xlim(0, 1); ax.set_ylim(0, 0.5)
ax.legend(title="$x$", loc="lower center", bbox_to_anchor=(0.5, 0.02), title_fontsize=7, ncol=2, columnspacing=1.0, handletextpad=0.3)
save(fig, "fig5_profiles")
