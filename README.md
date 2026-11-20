# Screened-Dimer-Electrostatics

Monomers at finite fugacity z in square-lattice dimer (domino) models, treated as a screened 2D Coulomb plasma
(charge = sublattice colour, coupling Gamma = 1/2, sine-Gordon beta^2 = 1/8), and their interplay with boundaries:
straight (neutral) edges, staircase (charged) edges, and the Aztec diamond, whose arctic region is melted by monomers.

Follow-on to the forced-defect electrostatics papers (Domino-Defects).

## Main results so far
- Bulk: n ~ z^(8/7); Stillinger-Lovett sum rule holds; charge mass m1^2/n -> 3.041(9) (sine-Gordon 7 sin(pi/7) = 3.0372;
  Debye-Hueckel pi excluded) from strip transfer matrices.
- Boundaries: straight edge -> same-sign image, boundary dimension Delta_b = q^2/2 (checked q = 1, 2);
  staircase edge -> charged plate, adsorbed monomer charge (1 - n)/2 per edge site (Gauss law), decaying as e^(-m1 d).
- Aztec diamond: dilute screening s ~ 5 z^2 2^N / N with rate ln 2 = exact dent asymptotics (2^N/(pi N) for a pair at the
  two tangent points); frozen order depends mainly on the screened fraction s; uncompensated edge charge
  Q_un - log2(1/z) looks z-independent for N* < N << xi (open question -> request/).

## Layout
| folder | contents |
|---|---|
| `src/` | core modules used by the figures and tests: exact solution (`elliptic_sol.py`, `exact_sol.py`), frozen-corner boundary (`frozen_boundary.py`), boundary Pfaffian in arbitrary precision (`bpf_exact.py`, needs `python-flint`), staircase strip and lozenge capacitor transfer matrices (`tm_diag_b.py`, `loz_cap.py`), Aztec sampler (`sample_aztec.py`), dimer entropy (`ckp.py`) |
| `scripts/` | working scripts by topic (`corner_islands/`, `isotherms/`, `boundary/`, `boundary_pfaffian/`, `bulk_strip/`, `mc_bulk/`, `aztec/`). Run from the repository root, e.g. `PYTHONPATH=src python3 scripts/corner_islands/bx_fit.py`; data paths point to `data/`. Only the scripts listed below were re-run after the files were moved. |
| `data/` | raw outputs behind the figures and tables |
| `figures/` | one script per paper figure (see `figures/README.md`); `figures/working/` holds older working plots |
| `manuscript/` | REVTeX source (`main.tex`, `refs.bib`); `make` builds `main.pdf` |
| `tests/` | `pytest -q tests` (about 3 s): tabulated values, series vs elliptic form, small-x law, Aztec count, boundary Pfaffian vs direct enumeration, strip isotherm |
| `request/`, `taskA_results/` | large-scale Monte Carlo requests and returned results |

Re-run after the move: all `figures/fig*.py`, `tests/`, `scripts/corner_islands/bx_fit.py`.

## Build
    gcc -O3 -march=native -o scripts/aztec/aztec_pt scripts/aztec/aztec_pt.c -lm
    gcc -O3 -march=native -o scripts/aztec/aztec_mc scripts/aztec/aztec_mc.c -lm
    gcc -O3 -march=native -o scripts/mc_bulk/md_mc scripts/mc_bulk/md_mc.c -lm
Python 3 with numpy, scipy, matplotlib.

## Large-scale request
See `request/README.md`. Validation first (`scripts/aztec/validate.sh`), then `scripts/aztec/run_pt.sh`, then `scripts/aztec/pt_analyze.py`.
