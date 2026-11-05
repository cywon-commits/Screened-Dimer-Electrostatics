import numpy as np, sys
from scipy.optimize import least_squares
from img_test import prof
for W in (16,20,22):
  for z in (0.002,0.005):
    n=prof(W,np.full(W,z)); r=n/n[W//2]; d=np.arange(W)
    def res(p):
        a,de,do,Ae,Ao=p
        out=[]
        for x in range(1,W//2):   # skip outermost site
            s=np.sin(np.pi*(x+ (de if x%2==0 else do))/(W+2*(de if x%2==0 else do)))
            out.append(np.log(r[x])-np.log((Ae if x%2==0 else Ao))-a*np.log(s))
        return out
    f=least_squares(res,[0.25,0.5,0.5,1,1],bounds=([-2,-0.9,-0.9,0.1,0.1],[2,5,5,10,10]))
    print(f'W {W} z {z}: alpha = {f.x[0]:.3f}  (Dirichlet/same-sign +0.25, Neumann/opposite -0.25)  rms {np.sqrt(np.mean(np.array(f.fun)**2)):.4f}',flush=True)
