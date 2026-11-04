import numpy as np, vp, pickle, sys
from scipy.optimize import minimize
L2=np.log(2)
class Outer:
    def __init__(self,M,eps):
        self.V=vp.VP(M); self.V.eps=eps; self.M=M; self.H=None
        B0=self.V.boundary(np.zeros(M)); self.mask=np.isnan(B0)
        # affine map rho -> boundary values: B = B0 + sum_j rho_j J_j
        self.J=np.array([np.nan_to_num(self.V.boundary(np.eye(M)[j])-B0) for j in range(M)])
        self.B0=np.nan_to_num(B0)
    def S(self,rho):
        V=self.V
        H0=self.H
        S,H,G=V.solve(rho,H0,1e-13)
        self.H=H; Gb=np.where(self.mask,0,G)
        dS=np.tensordot(self.J,Gb,axes=([1,2],[0,1]))
        return S,dS
def run(x,M=64,eps=0.02,rho0=None,Hinit=None):
    O=Outer(M,eps); O.H=Hinit
    du=2.0/M; half=M//2
    def unpack(r): return np.concatenate([r,r[::-1]])
    def f(r):
        rho=unpack(r); S,dS=O.S(rho)
        phi=S-2*x*L2*rho.sum()*du/2
        g=dS-2*x*L2*du/2; g=g[:half]+g[half:][::-1]
        return -phi,-g
    r0=rho0 if rho0 is not None else np.full(half,0.2)
    res=minimize(f,r0,jac=True,method='L-BFGS-B',bounds=[(0,1)]*half,options=dict(maxiter=500,ftol=1e-13,gtol=1e-9))
    rho=unpack(res.x); sigma=1-2*rho.sum()*du/2
    return sigma,rho,-res.fun,O.H,res
if __name__=='__main__':
    M=int(sys.argv[1]); eps=float(sys.argv[2]); out={}
    for x in [float(v) for v in (sys.argv[3].split(",") if len(sys.argv)>3 else ["0.3","0.4","0.5","0.7"])]:
        s,rho,phi,H,res=run(x,M,eps)
        out[x]=(s,rho,phi); print(M,eps,x,'sigma=%.4f'%s,'nit',res.nit,flush=True)
    pickle.dump(out,open(f"az/vp_M{M}_e{eps}_{len(sys.argv)}.pkl","wb"))
