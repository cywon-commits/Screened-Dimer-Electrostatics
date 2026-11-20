# Manuscript versions

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
