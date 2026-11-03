"""Parameter-free Thomas-Fermi (equilibrium-measure) solution for boundary-adsorbed monomers.
Per SW-edge dent at a = 1/2 + delta:   gain u(delta) = H(a) - ln2/2   (H natural binary entropy, exact dent asymptotics)
cost of a dent: x ln2 / 2 ;  interaction with another dent (same edge, its partner on SE included):
   ln|delta-delta'| + K(delta-delta')/2,   K from the exact two-pair determinant (N -> infinity).
Maximise  G[rho] = int (u - x ln2/2) rho + 1/2 int int rho rho' [ln|d| + K(d)/2],  0 <= rho <= 1.
Then 1 - sigma = 2 int rho.
"""
import numpy as np, pickle
from scipy.optimize import minimize
def Kinf():
    r={}
    for f in ('az/twopair_200_400.pkl','az/twopair_800.pkl'): r.update(pickle.load(open(f,'rb')))
    tab={}
    for N,rows in r.items():
        for d,ld,a,b,c,e in rows: tab.setdefault(round(d/N,3),{})[N]=(ld-a-b)-2*np.log(d/N)
    xs=[];ys=[]
    for dl,v in sorted(tab.items()):
        if 400 in v and 800 in v: xs.append(dl); ys.append(2*v[800]-v[400])   # Richardson in 1/N
    xs=np.array([0.0]+xs); ys=np.array([3*np.log(2)]+ys)
    return lambda d: np.interp(np.abs(d),xs,ys)
def H(a): a=np.clip(a,1e-15,1-1e-15); return -(a*np.log(a)+(1-a)*np.log(1-a))
def solve(x, M=400, K=None, quad=False):
    h=1.0/M; d=-0.5+h*(np.arange(M)+0.5)
    u=(np.log(2)/2-2*d**2) if quad else (H(0.5+d)-np.log(2)/2)
    f=u-x*np.log(2)/2
    D=d[:,None]-d[None,:]
    # cell-averaged log kernel: int over cell of ln|D| / h  ~ ln|D| off-diagonal, (ln(h/2)-1) on diagonal
    L=np.log(np.abs(D)+np.eye(M)); np.fill_diagonal(L,np.log(h/2)-1)
    Km=K(D) if K else 2.0*np.ones_like(D)
    G=(L+Km/2)*h*h
    obj=lambda r: -(f@r*h+0.5*r@G@r)
    grad=lambda r: -(f*h+G@r)
    r0=np.clip(f,0,None)
    res=minimize(obj,r0,jac=grad,bounds=[(0,1)]*M,method='L-BFGS-B',options=dict(maxiter=5000,ftol=1e-14,gtol=1e-12))
    rho=res.x; return d,rho,1-2*rho.sum()*h
