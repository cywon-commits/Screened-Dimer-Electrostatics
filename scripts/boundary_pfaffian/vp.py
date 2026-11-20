"""Variational principle (CKP) for the Aztec diamond with boundary-monomer (dent) density rho on the edges.
Rotated coords u=x+y, v=y-x in [-1,1]^2 (diamond |x|+|y|<=1); slope constraint |h_u|,|h_v|<=1.
lnZ/N^2 = max_H  int (ent_CKP/2)(grad H) dx dy  - 2 x ln2 int rho(delta) d delta."""
import numpy as np
from scipy.optimize import minimize
from ckp import ent
class VP:
    def __init__(self,M):
        self.M=M; self.du=2.0/M; g=np.linspace(-1,1,M+1); self.g=g
        self.H=None
    def boundary(self,rho):
        """rho: array of length M (per boundary segment, same on each edge, oriented W->S, S->E, E->N, N->W)."""
        M=self.M; du=self.du; c=np.concatenate([[0],np.cumsum((1-2*rho)*du)])   # rise along an edge
        B=np.full((M+1,M+1),np.nan)  # index [iu, iv]
        # SW: u=-1, v from 1 to -1 : H = c(k), v index M-k
        for k in range(M+1): B[0,M-k]=c[k]
        top=c[-1]
        for k in range(M+1): B[k,0]=top-c[k]          # SE: v=-1, u from -1 to 1
        for k in range(M+1): B[M,k]=c[k]              # NE: u=1, v from -1 to 1  (E->N rises)
        for k in range(M+1): B[M-k,M]=top-c[k]        # NW: v=1, u from 1 to -1 (N->W falls)
        return B
    def slopes(self,H):
        du=self.du; h00=H[:-1,:-1]; h10=H[1:,:-1]; h01=H[:-1,1:]; h11=H[1:,1:]
        return [((h10-h00)/du,(h11-h10)/du),((h11-h01)/du,(h01-h00)/du),((h10-h00)/du,(h01-h00)/du),((h11-h01)/du,(h11-h10)/du)]
    def S_and_grad(self,H,kap=200.0):
        ep=getattr(self,'eps',0.0); sc=1.0+ep
        du=self.du; w=du*du/4/2/2   # triangle (u,v)-area du^2/2, weight 1/2, dxdy=dudv/2, ent per square = ent_CKP/2
        S=0.0; G=np.zeros_like(H)
        sl=self.slopes(H)
        # node maps for each triangle's (hu, hv) = (a-b)/du , (c-d)/du
        maps=[(('10','00'),('11','10')),(('11','01'),('01','00')),(('10','00'),('01','00')),(('11','01'),('11','10'))]
        idx={'00':(slice(0,-1),slice(0,-1)),'10':(slice(1,None),slice(0,-1)),'01':(slice(0,-1),slice(1,None)),'11':(slice(1,None),slice(1,None))}
        for (hu,hv),((a,b),(c,d)) in zip(sl,maps):
            huc=np.clip(hu/sc,-1,1); hvc=np.clip(hv/sc,-1,1)
            s=huc-hvc; t=huc+hvc
            e,es,et=ent(s,t,True)
            eu=(es+et)/sc; ev=(-es+et)/sc
            pu=np.maximum(np.abs(hu)-sc,0); pv=np.maximum(np.abs(hv)-sc,0)
            S+=w*(e.sum()-kap*(pu**2+pv**2).sum())
            gu=w*(eu-2*kap*pu*np.sign(hu))/du; gv=w*(ev-2*kap*pv*np.sign(hv))/du
            G[idx[a]]+=gu; G[idx[b]]-=gu; G[idx[c]]+=gv; G[idx[d]]-=gv
        return S,G
    def solve(self,rho,H0=None,tol=1e-12):
        B=self.boundary(rho); mask=np.isnan(B)
        if H0 is None:
            # harmonic-ish initial guess: bilinear blend of boundary
            H0=np.where(mask,0,B); 
            for _ in range(2000):
                H0[1:-1,1:-1]=0.25*(H0[2:,1:-1]+H0[:-2,1:-1]+H0[1:-1,2:]+H0[1:-1,:-2])
        H=np.where(mask,H0,B)
        def f(z):
            H[mask]=z; S,G=self.S_and_grad(H); return -S,-G[mask]
        r=minimize(f,H[mask],jac=True,method='L-BFGS-B',options=dict(maxiter=20000,maxcor=30,ftol=tol,gtol=1e-11))
        H[mask]=r.x; S,G=self.S_and_grad(H)
        self.H=H; self.Gb=G   # boundary entries of G = dS/dH_boundary
        return S,H,G
if __name__=='__main__':
    for M in (32,64,128):
        V=VP(M); S,H,G=V.solve(np.zeros(M)); print(M,'lnM/N^2 =',S,' (exact ln2/2 =',np.log(2)/2,')')

def refine(H):
    M=H.shape[0]-1; from scipy.interpolate import RegularGridInterpolator as RGI
    g=np.linspace(-1,1,M+1); f=RGI((g,g),H); g2=np.linspace(-1,1,2*M+1)
    U,Vv=np.meshgrid(g2,g2,indexing='ij'); return f(np.stack([U,Vv],-1))
def solve_ml(rho_fun,Ms=(16,32,64,128),tol=1e-13,verbose=False,eps=0.0):
    H=None
    for M in Ms:
        V=VP(M); V.eps=eps; w=-1+V.du*(np.arange(M)+.5); rho=rho_fun(w)
        H0=refine(H) if H is not None else None
        S,H,G=V.solve(rho,H0,tol)
        if verbose: print(M,S,flush=True)
    return S,V,rho
