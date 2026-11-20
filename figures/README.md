# Paper figures

Each figure is made by one script, at final size (PRE two-column: 3.375 in or 7.0 in), with the common style in `style.py`.
Run from this folder; each script writes `<name>.pdf` (vector, for submission) and `<name>.png` (preview).

| Figure | Command | Input |
|---|---|---|
| 1 model, boundary types, Aztec diamond | `python3 fig1_model.py` | none (schematic) |
| 2 isotherms of uniform geometries | `python3 fig2_isotherms.py` | computed (`../isotherms`) |
| 3 corner charges, reflection, hodograph sector | `python3 fig3_corner_charges.py` | exact solution (`../corner_islands`) |
| 4 uncompensated fraction | `python3 fig4_sigma.py` | `../data/bx_exact.pkl` |
| 5 edge vacancy profiles | `python3 fig5_profiles.py` | `../data/smallx.pkl`, `../data/profiles.pkl` |
| 6 frozen corners | `python3 fig6_frozen.py` | exact solution |
| 7 full model at fixed x | (pending Monte Carlo data, see `../request/REQUEST_B_ko.md`) | |

The `az_*.png` files are older working plots, kept for the record; they are not paper figures.
