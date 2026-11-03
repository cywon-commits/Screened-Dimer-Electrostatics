import numpy as np, tf_solve as S
from scipy.special import k0
from scipy.optimize import minimize
K=S.Kinf()
def solve(x, Nxi, M=300, screen_cross=False):
    """constant-cross TF (10단계) with same-edge log kernel replaced by the screened 2D kernel:
       ln|d| -> -K0(Nxi*|d|) + ln(2/Nxi) - gamma   (Nxi = N/xi)."""
    h=1.0/M; d=-0.5+h*(np.arange(M)+0.5)
    f=S.H(0.5+d)-np.log(2)/2-x*np.log(2)/2
    D=np.abs(d[:,None]-d[None,:])
    g=np.euler_gamma
    if Nxi>0:
        Ds=np.where(D>0,D,h/4)
        L=-k0(Nxi*Ds)+np.log(2/Nxi)-g
        np.fill_diagonal(L,np.log(h/2)-1)
    else:
        L=np.log(D+np.eye(len(d))); np.fill_diagonal(L,np.log(h/2)-1)
    cross=(K(d[:,None]-d[None,:])+np.log(2))/2
    if screen_cross and Nxi>0: cross=cross*np.exp(-Nxi*0.7)   # crude: cross partners ~0.7 edge lengths away
    G=(L+cross)*h*h
    obj=lambda r: -(f@r*h+0.5*r@G@r); grad=lambda r: -(f*h+G@r)
    res=minimize(obj,np.clip(f,0,None),jac=grad,bounds=[(0,1)]*len(d),method='L-BFGS-B',options=dict(maxiter=8000,ftol=1e-15,gtol=1e-12))
    return 1-2*res.x.sum()*h
