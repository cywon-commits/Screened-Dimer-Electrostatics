"""Fig. 7: full model against boundary model.  (a) edge free energy dg(rho); (b) D = Q_un(full) - Q_un(boundary) vs 1/N
with the predicted limits.  Inputs (made by scripts/full_vs_boundary/): dg.csv, D_table.csv, D_inf.csv
(D_inf.csv: PYTHONPATH=src:scripts/full_vs_boundary python3 scripts/full_vs_boundary/predict.py 12,16 0.4 0.8,0.7,0.6,0.5,0.4152,0.3)."""
import csv, numpy as np, matplotlib.pyplot as plt
from style import *
P = "../data/full_vs_boundary/"
dg = list(csv.DictReader(open(P + "dg.csv"))); rho = np.array([float(r["rho"]) for r in dg])
col = lambda k: np.array([float(r[k]) if r[k] else np.nan for r in dg])
best = np.full(len(rho), np.nan)
for k in ("dg_W8", "dg_W10", "dg_W12", "dg_W14", "dg_W16"):       # the largest available width wins
    v = col(k); best[~np.isnan(v)] = v[~np.isnan(v)]
fig, ax = plt.subplots(1, 2, figsize=(TEXTW, 2.5), gridspec_kw=dict(wspace=0.30))
a = ax[0]; t = np.tan(np.pi * rho / 2)
a.plot(rho, best, "k-", label="one-sided strip, largest $W$"); a.plot(rho, col("dg_W8"), ":", color=GREY, label="$W=8$")
a.plot(rho, 0.545 * t * t, "--", color=VERM, label=r"$0.545\,\tan^2(\pi\rho/2)$")
a.set_xlabel(r"edge vacancy density $\rho$"); a.set_ylabel(r"edge free energy $\delta g$"); a.set_xlim(0, 0.42); a.set_ylim(0, 0.33); a.legend(loc="upper left")
D = [(int(r["N"]), float(r["x"]), float(r["D"]), float(r["err"]), r["source"]) for r in csv.DictReader(open(P + "D_table.csv"))]
extra = {0.4152: [(64, -1.467, 0.038, "MC")]}   # Q_un = 18.475 +- 0.038 (taskA_results/reports/interim4_N64pl.md) minus 19.942 (boundary, exact)
pred = {float(r["x"]): float(r["D_inf_W16"]) for r in csv.DictReader(open(P + "D_inf.csv"))}
MK = dict(MARK); MK[0.4152] = "o"; a = ax[1]
for x in (0.7, 0.5, 0.4152, 0.3):
    pts = [(N, d, e, "exact" if s == "exact" else "MC") for N, xx, d, e, s in D if abs(xx - x) < 2e-4] + extra.get(x, [])
    sel = {}
    for N, d, e, s in pts:
        if N not in sel or e < sel[N][1]: sel[N] = (d, e, s)
    Ns = np.array(sorted(sel)); c = xcolor(x)
    a.plot(1 / Ns, [sel[N][0] for N in Ns], "-", color=c, lw=0.6)
    for N in Ns:
        d, e, s = sel[N]; a.errorbar([1 / N], [d], [e], fmt=MK[x], color=c, mfc="w" if s == "exact" else c, ms=3.5, mew=0.7, capsize=1.5, lw=0.7)
    a.plot([], [], MK[x] + "-", color=c, ms=3.5, lw=0.6, label=f"${x}$")
    a.plot([0], [pred[x]], "*", color=c, ms=7, mec="k", mew=0.4, clip_on=False, zorder=5)
a.plot([], [], "*", color="w", mec="k", mew=0.4, ms=7, label="predicted limit")
a.set_xlim(0, 0.09); a.set_ylim(-2.75, -0.2); a.set_xlabel("$1/N$"); a.set_ylabel(r"$Q_{\rm un}^{\rm full}-Q_{\rm un}^{\rm bdry}$")
a.legend(title="$x$", title_fontsize=7, loc="lower right", ncol=2, columnspacing=1.0)
panel(ax[0], "(a)", x=-0.14); panel(ax[1], "(b)", x=-0.14)
save(fig, "fig7_full_model")
