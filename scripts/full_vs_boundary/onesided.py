"""Boundary free energy dg(rho) of a staircase edge in the limit z -> 0 at fixed edge density.

One-sided strip: vacancies of the left-edge colour (even diagonals) are allowed on every cell; vacancies of the other
colour only on the outermost cells of the right edge.  No neutral pair can form in the interior, so the limit z -> 0 is
reached up to terms of order t^W (a left-colour vacancy crossing the whole strip), t = tan(pi rho/2).
Reference: boundary model (both colours on outermost cells only).
Compared at equal vacancy number per period (Legendre transform), which removes the O(1/W) shift of the state:
    dg_W(rho) = phi_one(rho) - phi_bdry(rho),   phi(rho) = ln Lambda(z) - 2 rho ln z  at  z d lnLambda/dz = 2 rho."""
import sys, numpy as np, pickle, os, time
import tm_diag_b as T
from scipy.sparse.linalg import LinearOperator, eigs
from scipy.interpolate import CubicSpline
def set_model(m):
    if m=='bdry':
        T.zE=lambda i,z,W: z if i==0 else 0.0
        T.zO=lambda j,z,W: z if j==W-1 else 0.0
    elif m=='one':
        T.zE=lambda i,z,W: z
        T.zO=lambda j,z,W: z if j==W-1 else 0.0
    elif m=='full':
        T.zE=lambda i,z,W: z
        T.zO=lambda j,z,W: z
def lnlam(W,z,m):
    set_model(m); n=1<<W
    R=LinearOperator((n,n),matvec=lambda v: T.apply_O(T.apply_E(np.real(v).reshape(-1),z,W),z,W),dtype=float)
    v=np.ones(n)/np.sqrt(n)
    for _ in range(40):
        w=R.matvec(v); v=w/np.linalg.norm(w)
    vals,vecs=eigs(R,k=2,which='LR',v0=v,tol=1e-13,ncv=24)
    return float(np.log(vals.real.max()))
def curve(W,m,rho_lo=0.02,rho_hi=0.47,npts=46):
    """F(ln z) on a grid; returns (rho, phi) with rho from the spline derivative"""
    lt=np.linspace(np.log(np.tan(np.pi*rho_lo/2)),np.log(np.tan(np.pi*rho_hi/2)),npts)
    lz=lt*W/2; F=np.array([lnlam(W,float(np.exp(a)),m) for a in lz])
    cs=CubicSpline(lz,F); n=cs(lz,1)             # total vacancies per period = 2 rho
    return lz,F,n/2,F-n*lz
if __name__=='__main__':
    f='data/full_vs_boundary/onesided.pkl'; out=pickle.load(open(f,'rb')) if os.path.exists(f) else {}
    for W in [int(a) for a in sys.argv[1].split(',')]:
        for m in ('bdry','one'):
            if (W,m) in out: continue
            t0=time.time(); out[(W,m)]=curve(W,m); pickle.dump(out,open(f,'wb')); print(W,m,'%.0fs'%(time.time()-t0),flush=True)
