import numpy as np, vp_outer as VO, sys
from scipy.optimize import minimize
L2=np.log(2)
def run(x,N,M,eps,rho0=None):
    O=VO.Outer(M,eps); du=2/M; half=M//2
    unpack=lambda r: np.concatenate([r,r[::-1]])
    def f(r):
        rho=unpack(r); S,dS=O.S(rho)
        q=np.clip(rho,1e-12,1-1e-12); import os; mode=os.environ.get('MIX','bin')
        if mode=='bin': s=-(q*np.log(q)+(1-q)*np.log(1-q)); ds=np.log((1-q)/q)
        else: s=-0.5*q*np.log(q); ds=-0.5*(np.log(q)+1)
        mix=(4.0/N)*s.sum()*du/2 if N else 0.0
        phi=S-2*x*L2*rho.sum()*du/2+mix
        g=dS-2*x*L2*du/2+((4.0/N)*ds*du/2 if N else 0); g=g[:half]+g[half:][::-1]
        return -phi,-g
    r0=rho0 if rho0 is not None else np.full(half,0.2)
    res=minimize(f,r0,jac=True,method='L-BFGS-B',bounds=[(1e-9,1-1e-9)]*half,options=dict(maxiter=500,ftol=1e-14,gtol=1e-10))
    rho=unpack(res.x); return 1-2*rho.sum()*du/2, res.x
M=int(sys.argv[1])
for x in [float(v) for v in sys.argv[2].split(',')]:
    out={}
    for eps in (0.02,0.01):
        s0,r0=run(x,0,M,eps)
        for N in (400,):
            sN,_=run(x,N,M,eps,r0); out[(eps,N)]=(N/2)*(sN-s0)
    c={N:2*out[(0.01,N)]-out[(0.02,N)] for N in (400,)}
    print('x',x,'c(N=400) %.3f'%c[400],' raw',{k:round(v,3) for k,v in out.items()},flush=True)
