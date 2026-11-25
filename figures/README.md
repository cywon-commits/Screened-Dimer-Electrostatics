# Paper figures

Each figure is made by one script, at final size (PRE two-column: 3.375 in or 7.0 in), with the common style in `style.py`.
Run from this folder; each script writes `<name>.pdf` (vector, for submission) and `<name>.png` (preview).

| Figure | Command | Input |
|---|---|---|
| 1 model, boundary types, Aztec diamond | `python3 fig1_model.py` | `../data/fig1_sample.pkl` (made by `python3 ../src/sample_aztec.py sample 11` run from `figures/`; `... check 3` compares the sampler with exact edge densities) |
| 2 isotherms of uniform geometries | `python3 fig2_isotherms.py` | computed (`../src`) |
| 3 corner charges, reflection, hodograph sector | `python3 fig3_corner_charges.py` | exact solution (`../src`) |
| 4 uncompensated fraction | `python3 fig4_sigma.py` | `../data/bx_exact.pkl` |
| 5 edge vacancy profiles | `python3 fig5_profiles.py` | `../data/smallx.pkl`, `../data/profiles.pkl` |
| 6 frozen corners | `python3 fig6_frozen.py` | exact solution |
| 7 full model against boundary model | `python3 fig7_full_model.py` | `../data/full_vs_boundary/dg.csv`, `D_table.csv`, `D_inf.csv` (made by `../scripts/full_vs_boundary/`) |

The files in `working/` are older working plots, kept for the record; they are not paper figures.
