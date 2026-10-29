# Monomer-dimer Aztec diamond: code package for the large-scale request

Square-lattice dimers with monomer fugacity z (weight z^{#monomers}) on the Aztec diamond of order N.
Monomers carry charge +1 (even sublattice) / -1 (odd); each staircase edge carries 1/2 per edge site.
Main observable: uncompensated edge charge Q_un = N/2 - Q_edge (Q_edge = signed monomer charge in the top
triangle cut by the two diagonals), screened fraction s = 2 Q_edge / N, frozen order Psi.

## Build
    gcc -O3 -march=native -o aztec_pt aztec_pt.c -lm
    gcc -O3 -march=native -o aztec_mc aztec_mc.c -lm
Python 3 with numpy + scipy for the analysis and exact-reference scripts.

## Step 0 - validation (mandatory)
    ./validate.sh
Must reproduce the exact N=10 values printed by the script within 2 sigma, and the N=20 references in the
request document (Psi within 0.002, Q_un within 2 sigma).

## Step 1 - production (task A)
    ./run_pt.sh N SEED R ZMIN ZMAX THERM MEAS OUTDIR [NBLK] [EVERY]
e.g. ./run_pt.sh 64 1 42 1e-5 0.1 6000000 24000000 out 200 200
Run 4 independent seeds per N.
aztec_pt checkpoints to OUTDIR/N<N>_s<SEED>.ckpt every AZTEC_CKPT_SEC seconds (default 300); if a run is
killed, rerun the identical command and it resumes (output is bitwise identical to an uninterrupted run).  Increase R until every neighbour swap acceptance printed at the end
(file OUTDIR/N<N>_s<SEED>.ladder) is >= 0.2.

## Step 2 - analysis
    python3 pt_analyze.py OUTDIR      -> OUTDIR/summary.csv
Flags runs whose first/second-half drift exceeds 3 sigma ("not thermalised").

## Files
- aztec_pt.c   replica exchange in z; local moves: break/join/slide of monomers, plaquette flips
- aztec_mc.c   single-z version of the same local dynamics (cross-check)
- pt_analyze.py  Q_un, s, Psi, errors (seed-to-seed scatter), drift test
- aztec.py     exact transfer computation for N <= 20 (diagonal frame); pairs.py fixed monomer pairs
- tm2.py, tm_open.py, tm_diag.py, md_mc.c   bulk / strip references (n(z), charge mass m1, boundary results)

## Output files of aztec_pt
- <prefix>_r<k>.txt       '# N z' header, then block lines 'B Q <Q_edge> n <monomer density>'
- <prefix>_r<k>_hist.txt  per cell: x y count(monomer) count(+x) count(+y) count(-x) count(-y)
