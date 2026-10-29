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
| `bulk_strip/` | cylinder transfer matrices (`tm2.py`: n, m1, m2 by momentum sector) |
| `mc_bulk/` | grand-canonical bulk MC (`md_mc.c`) and correlator analysis |
| `boundary/` | open strip (`tm_open.py`, straight edges) and diagonal strip (`tm_diag.py`, staircase edges) |
| `aztec/` | exact Aztec transfer (`aztec.py`, N <= 20), fixed pairs (`pairs.py`), MC (`aztec_mc.c`), parallel tempering (`aztec_pt.c`), analysis |
| `request/` | large-scale computation request (task A: N = 32-128, z = 1e-5 - 0.1) |
| `data/` | raw outputs behind the figures and tables |
| `figures/` | Aztec maps, screening curves, scaling plots |

## Build
    gcc -O3 -march=native -o aztec/aztec_pt aztec/aztec_pt.c -lm
    gcc -O3 -march=native -o aztec/aztec_mc aztec/aztec_mc.c -lm
    gcc -O3 -march=native -o mc_bulk/md_mc mc_bulk/md_mc.c -lm
Python 3 with numpy, scipy, matplotlib.

## Large-scale request
See `request/README.md`. Validation first (`aztec/validate.sh`), then `aztec/run_pt.sh`, then `aztec/pt_analyze.py`.
