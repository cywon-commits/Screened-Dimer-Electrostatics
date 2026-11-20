import numpy as np
from exact_sol import cheb
def arctic(JM, m=600):
    c0,c1,a=cheb(JM); n=np.arange(len(a))
    T=np.linspace(1e-4,1-1e-4,m); J=1-8*T**2/(1+T**2)**2
    Bp=(c0-J)/c1; ep=Bp-np.sqrt(Bp*Bp-1)
    Th=-2**-0.25*(1-J)**-0.25*(a[None,:]*ep[:,None]**n[None,:]).sum(1)
    Psi=1+T+np.sqrt(T)*Th
    dP=np.gradient(Psi,T); return dP, Psi-T*dP      # (u,v)
