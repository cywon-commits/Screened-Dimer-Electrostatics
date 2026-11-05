import numpy as np
from scipy.optimize import least_squares
from img_test import prof
W=20; z=0.002
for zb in (0.002,0.02,0.1,0.3,1.0,3.0,10.0):
    zs=np.full(W,z); zs[0]=zs[-1]=zb
    n=prof(W,zs); r=n/n[W//2]
    def res(p,skip):
        a,de,do,Ae,Ao=p; out=[]
        for x in range(skip,W//2):
            dd=de if x%2==0 else do
            out.append(np.log(r[x])-np.log(Ae if x%2==0 else Ao)-a*np.log(np.sin(np.pi*(x+dd)/(W+2*dd))))
        return out
    al=[least_squares(lambda p: res(p,s),[0,0.5,0.5,1,1],bounds=([-2,-0.9,-0.9,0.01,0.01],[2,5,5,10,10])).x[0] for s in (1,3)]
    print(f'zb {zb:6.3f}: n_edge/n_mid {r[0]:.3f}  alpha(skip1) {al[0]:+.3f}  alpha(skip3) {al[1]:+.3f}   env ' + ' '.join(f'{(r[2*k]+r[2*k+1])/2:.3f}' for k in range(5)),flush=True)
