#!/bin/bash
# Step 0 before production: must reproduce the exact transfer-matrix values (aztec.py) within 2 sigma.
set -e
gcc -O3 -march=native -o aztec_pt aztec_pt.c -lm
gcc -O3 -march=native -o aztec_mc aztec_mc.c -lm
mkdir -p val
./aztec_pt 10 50000 6000000 60 11 50 val/N10_s11 0.001 0.0015 0.002 0.003 0.005 0.007 0.01 0.015 0.02 0.03 0.05 0.07 0.1
python3 pt_analyze.py val
echo "exact (N=10): z=0.001 Q_un 4.9974 | 0.003 4.9775 | 0.01 4.8024 | 0.03 4.2269 | 0.1 3.0638 ;  Psi(0.01)=0.4280 Psi(0.1)=0.2316"
python3 -c "
import aztec as A, numpy as np
for z in (0.003,0.01,0.1):
    lnZ,m,l,r=A.Aztec(16,z).run()
    print('exact N=16 z',z,'lnZ',round(lnZ,6))
"
