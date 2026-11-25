# Manuscript versions

## v0.3 (2026-11-25) leading term of the edge free energy
- dg(rho) = (4/3) rho^2 + O(rho^3): Eq. (dglead) in Sec. VI C, derivation in Appendix F (wall picture, COE, Pfaffian
  evaluation, closed form kappa_n). Fig. 7(a) dashed curve is now this leading term; the fitted 0.545 tan^2 is removed.
- Status in the text: the closed form for kappa_n is identified from exact values and checked to n = 40, not proved;
  the order of limits (width before length) is assumed; no closed form for the full function.

## v0.2 (2026-11-25) full model
- Sec. VI C rewritten: the difference between the full and the boundary model is an edge free energy dg(rho);
  D_inf(x) = (1/ln 2) dI/dx, Eqs. (dg), (Dinf). Evidence level stated in the text: locality is an argument; predicted
  constants agree with Monte Carlo at x = 0.7, 0.5, 0.4152 (two sizes each in the asymptotic range); x = 0.3 is not a test.
- Table III replaced by the fixed-x series (request B + task A), with D and the predicted D_inf.
- New Fig. 7 (`figures/fig7_full_model.py`) and Appendix F (one-sided strips).
- Source of the N = 64, z = 1e-4 point fixed: 18.475 +- 0.038 (task A interim report 4); v0.1 quoted "about 0.02".
- Sec. II A, Discussion (range of validity, open problems), Appendix E (full model) updated accordingly. 13 pages.
- Still open from v0.1: author list, arXiv:2506.15769, 13 `TODO-` references, definition of Q_edge against the code.
- Not blind: the predictions were compared with data that already existed. Predictions recorded before any data:
  D_inf(0.6) = -0.586, D_inf(0.8) = -0.234 (`data/full_vs_boundary/D_inf.csv`).

## v0.1 (2026-11-21) first full draft
- All sections I-VIII and Appendices A-E written; compiles with REVTeX 4.2 (`pre`, `reprint`), 12 pages, no warnings.
- No new calculations except the sampler behind Fig. 1(c) (`src/sample_aztec.py`, checked at N=24, x=0.5).
- Figures 1-6 from `../figures`. Fig. 7 (full model at fixed x) is pending Monte Carlo data (`request/REQUEST_B_ko.md`).
- Convention fixed in this draft: `ln Z / N^2 = (1/2) int s(grad H) du dv + (x ln2/4)[H(S)+H(N)-H(W)-H(E)]`, with `s` the entropy
  per cell, (u,v) along the edges with half edge = 1, and |dH/du| = 1 - 2 rho on an edge. (The research notes wrote the same
  functional in Cartesian area units.)
- Symbols: `r` = modulus of zeta on the real side (notes: rho), `rho` = vacancy density, `omega`, `eta` = Burgers variables (notes: z, w).

### Open items before review
- Author list; whether to cite the companion papers of the series.
- Reference for arXiv:2506.15769 (bibliographic data not confirmed; left as a comment in Sec. I).
- 13 bibliography keys with prefix `TODO-` are entered from memory (see `refs_check.md`).
- Statistical errors of the Monte Carlo values in Table III are quoted as "about 0.02 or less"; replace by the run-by-run errors.
- Definition of Q_edge in Sec. VI C should be checked against `scripts/aztec/aztec_pt.c` (region weights `wtop`).
- Six of the seven figures are two-column; consider reducing Figs. 2 and 4 to one column.
