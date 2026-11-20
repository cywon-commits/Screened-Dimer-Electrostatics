"""TF with position-dependent exact asymptotic kernels on all four edges.
SW dent at p (fraction from bottom corner), symmetric density rho(q) on every edge (q from that edge's own corner):
 H(p)-ln2/2 - x ln2/2 + int rho(q) [ ln|p-q| + X - ln|D(p,q)| - ln|D(1-p,q)| - ln2 ] dq = 0 on support
 D(p,q) = 1/2 - (1-p)(1-q)  (crossing-regime asymptotics of the dent kernel), X = far same-sign term.
"""
import numpy as np
from scipy.optimize import minimize
def H(a): a=np.clip(a,1e-15,1-1e-15); return -(a*np.log(a)+(1-a)*np.log(1-a))
def solve(x, M=300, X=-np.log(2), Dmin=0.02):
    h=1.0/M; p=h*(np.arange(M)+0.5)
    f=H(p)-np.log(2)/2-x*np.log(2)/2
    P=p[:,None]; Q=p[None,:]
    L=np.log(np.abs(P-Q)+np.eye(M)); np.fill_diagonal(L,np.log(h/2)-1)
    D1=np.maximum(0.5-(1-P)*(1-Q),Dmin); D2=np.maximum(0.5-P*(1-Q),Dmin)
    G=(L+X-np.log(D1)-np.log(D2)-np.log(2))
    # symmetric kernel needed for a variational form: symmetrise (rho is symmetric about 1/2 anyway)
    G=0.5*(G+G.T)*h*h
    obj=lambda r: -(f@r*h+0.5*r@G@r); grad=lambda r: -(f*h+G@r)
    r0=np.clip(f,0,None)
    res=minimize(obj,r0,jac=grad,bounds=[(0,1)]*M,method='L-BFGS-B',options=dict(maxiter=8000,ftol=1e-15,gtol=1e-12))
    return p,res.x,1-2*res.x.sum()*h
