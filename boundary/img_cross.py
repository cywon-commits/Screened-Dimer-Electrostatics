import numpy as np
from scipy.optimize import least_squares, brentq
from img_test import prof
z=0.002
def alpha(W,zb):
    zs=np.full(W,z); zs[0]=zs[-1]=zb; n=prof(W,zs); r=n/n[W//2]
    def res(p):
        a,de,do,Ae,Ao=p; out=[]
        for x in range(1,W//2):
            dd=de if x%2==0 else do
            out.append(np.log(r[x])-np.log(Ae if x%2==0 else Ao)-a*np.log(np.sin(np.pi*(x+dd)/(W+2*dd))))
        return out
    return least_squares(res,[0,0.5,0.5,1,1],bounds=([-2,-0.9,-0.9,0.01,0.01],[2,5,5,10,10])).x[0]
for W in (10,12,14,16,18,20):
    zs=brentq(lambda lz: alpha(W,np.exp(lz)), np.log(0.005), np.log(0.5), xtol=1e-3)
    print(f'W {W}: crossover zb* = {np.exp(zs):.4f}   zb*^2 * W = {np.exp(2*zs)*W:.4f}',flush=True)
