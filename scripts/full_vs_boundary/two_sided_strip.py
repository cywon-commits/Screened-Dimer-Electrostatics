import sys, numpy as np, pickle, os, time
import tm_diag_b as T
from scipy.sparse.linalg import LinearOperator, eigs
def lam(W,z,rin):
    T.RIN=rin; n=1<<W
    R=LinearOperator((n,n),matvec=lambda v: T.apply_O(T.apply_E(np.real(v).reshape(-1),z,W),z,W),dtype=float)
    # power iteration first (positive vector), then Arnoldi started from it
    v=np.ones(n)/np.sqrt(n)
    for _ in range(60):
        w=R.matvec(v); v=w/np.linalg.norm(w)
    vals,vecs=eigs(R,k=2,which='LR',v0=v,tol=1e-12,ncv=24)
    k=np.argmax(vals.real); psi=np.real(vecs[:,k]); psi*=np.sign(psi.sum())
    res=np.linalg.norm(R.matvec(psi)-vals[k].real*psi)/np.linalg.norm(psi)
    return vals[k].real, res, psi.min()/psi.max()
f='data/full_vs_boundary/strip_dg2.pkl'; out=pickle.load(open(f,'rb')) if os.path.exists(f) else {}
rhos=[float(a) for a in sys.argv[2].split(',')]
for W in [int(a) for a in sys.argv[1].split(',')]:
    for r in rhos:
        if (W,r) in out: continue
        t0=time.time(); t=np.tan(np.pi*r/2); z=t**(W/2)
        a=lam(W,z,1.0); b=lam(W,z,0.0)
        out[(W,r)]=(np.log(a[0])-np.log(b[0]),np.log(a[0]),np.log(b[0]),a[1],b[1],a[2],b[2])
        print(W,r,'%.6f'%out[(W,r)][0],'res %.1e %.1e  minpsi %.1e %.1e  %.0fs'%(a[1],b[1],a[2],b[2],time.time()-t0),flush=True)
        pickle.dump(out,open(f,'wb'))
