"""Fig. 1: vacancies as sublattice charges, the two boundary types, and an Aztec diamond with boundary vacancies."""
import numpy as np, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from scipy.optimize import linear_sum_assignment
from style import *
rng = np.random.default_rng(7)
def matching(cells, removed=()):
    cells = [c for c in cells if c not in removed]; S = set(cells)
    Bk = [c for c in cells if (c[0] + c[1]) % 2 == 0]; Wh = [c for c in cells if (c[0] + c[1]) % 2 == 1]
    assert len(Bk) == len(Wh), (len(Bk), len(Wh))
    iw = {c: i for i, c in enumerate(Wh)}; cost = np.full((len(Bk), len(Wh)), 1e6)
    for i, (x, y) in enumerate(Bk):
        for d in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            n = (x + d[0], y + d[1])
            if n in S: cost[i, iw[n]] = rng.random()
    r, c = linear_sum_assignment(cost); assert cost[r, c].max() < 10
    return [(Bk[i], Wh[j]) for i, j in zip(r, c)]
def draw_tiling(ax, cells, dom, removed, lw=0.5):
    for (x, y) in cells:
        if (x + y) % 2 == 0: ax.add_patch(Rectangle((x, y), 1, 1, fc="#F2F2F2", ec="none", zorder=0))
    for (a, b) in dom:
        x0, y0 = min(a[0], b[0]), min(a[1], b[1]); w, h = (2, 1) if a[1] == b[1] else (1, 2)
        ax.add_patch(FancyBboxPatch((x0 + .08, y0 + .08), w - .16, h - .16, boxstyle="round,pad=0,rounding_size=0.12", fc="none", ec="k", lw=lw, zorder=2))
    for (x, y) in removed:
        ax.add_patch(Rectangle((x + .12, y + .12), .76, .76, fc=C_POS if (x + y) % 2 == 0 else C_NEG, ec="none", zorder=3))
    ax.set_aspect("equal"); ax.axis("off")

def sq(ax, x, y, lw=0.3, dy=0.0):
    ax.add_patch(Rectangle((x, y + dy), 1, 1, fc="#F2F2F2" if (x + y) % 2 == 0 else "w", ec="#BBBBBB", lw=lw, zorder=0))
def vac(ax, x, y, sign=True, dy=0.0):
    c = C_POS if (x + y) % 2 == 0 else C_NEG
    ax.add_patch(Rectangle((x + .12, y + dy + .12), .76, .76, fc=c, ec="none", zorder=3))
    if sign: ax.text(x + .5, y + dy + .47, "$+$" if (x + y) % 2 == 0 else "$-$", color="w", ha="center", va="center", fontsize=7, zorder=4)
fig = plt.figure(figsize=(TEXTW, 2.6)); Y = 0.955
# (a) bulk patch with a +/- vacancy pair
ax = fig.add_axes([0.035, 0.06, 0.27, 0.84]); fig.text(0.005, Y, "(a)", fontsize=9, va="top")
cells = [(x, y) for x in range(8) for y in range(8)]; rem = [(2, 4), (5, 2)]
draw_tiling(ax, cells, matching(cells, rem), []); [vac(ax, *c) for c in rem]; ax.set_xlim(-.2, 8.2); ax.set_ylim(-.2, 8.2)
# (b) boundary types
ax = fig.add_axes([0.37, 0.06, 0.24, 0.84]); fig.text(0.335, Y, "(b)", fontsize=9, va="top")
for x in range(8):
    for y in range(3): sq(ax, x, y, dy=8.6)
    vac(ax, x, 0, dy=8.6)
ax.plot([0, 8], [8.6, 8.6], "k-", lw=1.0); ax.text(4, 8.2, "straight edge (neutral)", ha="center", va="top", fontsize=7)
c0 = 5; xs, ys = [], []
for x in range(8):
    for y in range(6):
        if x + y >= c0: sq(ax, x, y)
        if x + y == c0: vac(ax, x, y)
for x in range(0, 6): xs += [x, x, x + 1]; ys += [c0 - x + 1, c0 - x, c0 - x]
ax.plot(xs[1:], ys[1:], "k-", lw=1.0); ax.text(4, -0.4, "staircase edge (charged)", ha="center", va="top", fontsize=7)
ax.set_xlim(-.3, 8.3); ax.set_ylim(-1.6, 11.8); ax.set_aspect("equal"); ax.axis("off")
# (c) Aztec diamond with boundary vacancies
ax = fig.add_axes([0.66, 0.0, 0.34, 0.96]); fig.text(0.655, Y, "(c)", fontsize=9, va="top")
N = 10; cells = [(x, y) for y in range(2 * N) for x in range(2 * N) if abs(x - N + .5) + abs(y - N + .5) <= N]
def edge_cell(e, p): return [(N - 1 - p, p), (N + p, p), (N + p, 2 * N - 1 - p), (N - 1 - p, 2 * N - 1 - p)][e]
rem = [edge_cell(e, p) for e in range(4) for p in (3, 5, 6)]
draw_tiling(ax, cells, matching(cells, rem), [], lw=0.4); [vac(ax, *c, sign=False) for c in rem]
for t, (x, y), ha, va in (("S", (N, -0.4), "center", "top"), ("N", (N, 2 * N + .4), "center", "bottom"), ("W", (-0.4, N), "right", "center"), ("E", (2 * N + .4, N), "left", "center")):
    ax.text(x, y, t, ha=ha, va=va, fontsize=8)
ax.set_xlim(-1.8, 2 * N + 1.8); ax.set_ylim(-1.8, 2 * N + 1.8)
save(fig, "fig1_model")
