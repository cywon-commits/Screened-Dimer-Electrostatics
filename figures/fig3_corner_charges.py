"""Fig. 3: corner charges, their reflected checkerboard, and the hodograph domain."""
import sys, numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Wedge
from scipy.optimize import brentq
sys.path.insert(0, "../corner_islands")
from style import *
from exact_sol import solve
from exact_more_lib import arctic
from elliptic_sol import Sol
X0 = 0.8
lj = brentq(lambda l: solve(1 + np.exp(l))["x"] - X0, -12, 11, xtol=1e-10); u, v = arctic(1 + np.exp(lj))
rM = brentq(lambda r: Sol(r).x() - X0, 0.05, 0.95, xtol=1e-10)
U = np.concatenate([[0], u, v[::-1], [0]]); V = np.concatenate([[0], v, u[::-1], [0]])     # frozen corner at S, square [0,2]^2
rot = lambda a, b: (np.asarray(a) - np.asarray(b), np.asarray(a) + np.asarray(b))           # S bottom, E right, N top, W left
C_EDGE, C_AXIS, C_DIAG, C_ARC = "k", GREEN, PURPLE, ORANGE
def diamond(ax, a0=0, b0=0, lw=0.8, fill=True, ec="k", zo=3):
    ax.plot(*rot([a0, a0 + 2, a0 + 2, a0, a0], [b0, b0, b0 + 2, b0 + 2, b0]), "-", color=ec, lw=lw, zorder=zo)
    if fill:
        for (a, b) in ((U, V), (2 - U, V), (U, 2 - V), (2 - U, 2 - V)):
            ax.add_patch(Polygon(np.c_[rot(a0 + a, b0 + b)], fc=C_FROZEN, ec="none", zorder=1))
def charge(ax, a, b, s, ms=5.5):
    x, y = rot(a, b); c = C_POS if s > 0 else C_NEG
    ax.plot([x], [y], "o", ms=ms, mfc=c, mec="w", mew=0.5, zorder=5, clip_on=False)
    ax.text(x, y - 0.02, "$+$" if s > 0 else "$-$", color="w", ha="center", va="center", fontsize=6, zorder=6)
fig = plt.figure(figsize=(TEXTW, 2.45)); Y = 0.965
# (a) one diamond
ax = fig.add_axes([0.03, 0.02, 0.31, 0.92]); fig.text(0.005, Y, "(a)", fontsize=9, va="top")
diamond(ax)
ds = u[np.argmin(v)]; ax.plot(*rot([ds, 1], [0, 0]), "-", color=C_EDGE, lw=1.8, zorder=4); ax.plot(*rot([1, 1], [0, 1]), "-", color=C_AXIS, lw=1.2, zorder=4)
dd = u[np.argmin(abs(u - v))]; ax.plot(*rot([dd, 1], [dd, 1]), "-", color=C_DIAG, lw=1.2, zorder=4); ax.plot(*rot(u, v), "-", color=C_ARC, lw=1.4, zorder=4)
for (a, b, s) in ((0, 0, 1), (2, 2, 1), (2, 0, -1), (0, 2, -1)): charge(ax, a, b, s)
ax.text(0, -0.22, r"S  $(+q)$", ha="center", va="top"); ax.text(0, 4.22, r"N  $(+q)$", ha="center", va="bottom")
ax.text(2.18, 2.0, "E\n$(-q)$", ha="left", va="center", linespacing=1.1)
ax.text(-2.18, 2.0, "W\n$(-q)$", ha="right", va="center", linespacing=1.1); ax.text(-2.9, -0.55, r"$q=x\ln 2$", ha="left", va="bottom", fontsize=7)
ax.text(1.1, 0.9, "M", ha="left", va="top"); ax.text(0, 2.12, "O", ha="center", va="bottom")
ax.text(0, 2.75, "liquid", ha="center", va="center", fontsize=7); ax.text(0, 3.55, "frozen", ha="center", va="center", fontsize=7)
ax.set_xlim(-3.0, 3.0); ax.set_ylim(-0.75, 4.75); ax.set_aspect("equal"); ax.axis("off")
# (b) reflected checkerboard
ax = fig.add_axes([0.37, 0.02, 0.31, 0.92]); fig.text(0.355, Y, "(b)", fontsize=9, va="top")
R = (-1, 0, 1)
for i in R:
    for j in R: diamond(ax, 2 * i, 2 * j, lw=0.4, ec=GREY)
diamond(ax, 0, 0, lw=1.0, fill=False, zo=4)
for i in range(-1, 3):
    for j in range(-1, 3): charge(ax, 2 * i, 2 * j, (-1) ** (i + j), ms=4.5)
ax.set_xlim(-6.4, 6.4); ax.set_ylim(-4.6, 8.6); ax.set_aspect("equal"); ax.axis("off")
# (c) hodograph sector
ax = fig.add_axes([0.72, 0.10, 0.27, 0.80]); fig.text(0.705, Y, "(c)", fontsize=9, va="top")
ax.add_patch(Wedge((0, 0), 1, -45, 0, fc="#F4F4F4", ec="none"))
th = np.linspace(-np.pi / 4, 0, 100)
ax.plot([0, rM], [0, 0], "-", color=C_AXIS, lw=1.2); ax.plot([rM, 1], [0, 0], "-", color=C_EDGE, lw=1.6)
ax.plot(np.cos(th), np.sin(th), "-", color=C_ARC, lw=1.4); ax.plot([0, np.cos(np.pi / 4)], [0, -np.sin(np.pi / 4)], "-", color=C_DIAG, lw=1.2)
ax.plot([0, rM], [0, 0], "o", ms=3, mfc="w", mec="k", mew=0.6, zorder=5)
ax.text(-0.03, 0.03, "O", ha="right", va="bottom"); ax.text(rM, 0.04, "M", ha="center", va="bottom")
ax.text(rM / 2, 0.04, "$u=1$", ha="center", va="bottom", fontsize=7, color=C_AXIS); ax.text((1 + rM) / 2 + 0.03, 0.04, "$v=0$", ha="center", va="bottom", fontsize=7)
ax.text(0.30, -0.40, "$u=v$", ha="right", va="top", fontsize=7, color=C_DIAG)
ax.text(0.96, -0.42, "frozen\nboundary", ha="left", va="center", fontsize=7, color=C_ARC, linespacing=0.9)
ax.text(0.58, -0.24, r"$\zeta$", ha="center", va="center", fontsize=9); ax.text(1.0, 0.04, "1", ha="center", va="bottom", fontsize=7); ax.text(rM, -0.04, r"$\rho_M$", ha="center", va="top", fontsize=7)
ax.set_xlim(-0.15, 1.42); ax.set_ylim(-0.85, 0.22); ax.set_aspect("equal"); ax.axis("off")
save(fig, "fig3_corner_charges")
print("rho_M =", rM)
