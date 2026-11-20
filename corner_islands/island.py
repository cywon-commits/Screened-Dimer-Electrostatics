import numpy as np
from ckp import ent
# 1. stiffness: ent(s,t) ~ ent0 - a (s^2+t^2)
h=0.02; a=-(ent(h,0)-2*ent(0,0)+ent(-h,0))/h**2/2
print('stiffness a = %.5f   pi/16 = %.5f'%(a,np.pi/16))
# 2. lattice constant M0 = sum' (-1)^(a+b) ln|n|  (Gaussian-damped, eps->0)
n=np.arange(-600,601); A,B=np.meshgrid(n,n,indexing='ij'); r2=(A*A+B*B).astype(float); sg=(-1.0)**(A+B)
r2[600,600]=1.0; lnr=0.5*np.log(r2); lnr[600,600]=0
for eps in (4e-4,2e-4,1e-4,5e-5):
    w=np.exp(-eps*r2); w[600,600]=0
    print('eps',eps,'M0 =',(sg*lnr*w).sum(),' sum charges =',(sg*w).sum())
