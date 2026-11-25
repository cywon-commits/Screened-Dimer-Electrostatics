import numpy as np, sys
from qb import Qb
sig={0.3:0.520891,0.5:0.720715,0.7:0.867527,0.4152:0.643450,0.2:0.393292}
for N in (64,96,128):
    for x in (0.2,0.3,0.5,0.7):
        q=Qb(N,2.0**(-x*N/2)); print(f'N={N} x={x}: Q_bdry={q:.5f}  Q-(N/2)sigma={q-N/2*sig[x]:+.5f}',flush=True)
