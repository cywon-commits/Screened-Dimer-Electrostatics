"""Adsorption isotherm of the staircase strip: prediction rho = (2/pi) arctan(z^(2/W)),
lnLam = W*ent(rho) + 2 rho ln z, ent(rho) = (2/pi)[L(pi(1-rho)/2)+L(pi rho/2)]  (CKP at tilt (h_u,h_v)=(0,1-2rho))."""
import numpy as np, sys, tm_diag_b as T
from ckp import L
def ent(r): return (2/np.pi)*(L(np.pi*(1-r)/2)+L(np.pi*r/2))
for rin,label in ((0.0,'edge-only'),(1.0,'full')):
    T.RIN=rin
    for t in (0.2,0.5,0.8):
        row=[]
        for W in (6,8,10,12,14):
            z=t**(W/2); Lam,_,ne,no=T.densities(W,z)
            row.append((W,ne[0],no[W-1],np.log(Lam)/W))
        pred=(2/np.pi)*np.arctan(t); fpred=ent(pred)+2*pred*np.log(t)/2   # lnLam/W = ent + 2 rho ln z / W = ent + rho ln t
        print(label,'t=z^(2/W)=',t,' pred rho %.4f  pred lnLam/W %.4f'%(pred,fpred))
        for W,a,b,f in row: print('   W %2d  rho_edge %.4f %.4f   lnLam/W %.4f'%(W,a,b,f))
        Ws=np.array([r[0] for r in row[-3:]]); rs=np.array([r[1] for r in row[-3:]]); c=np.polyfit(1/Ws,rs,1)
        print('   1/W extrapolation of rho: %.4f'%c[1],flush=True)
