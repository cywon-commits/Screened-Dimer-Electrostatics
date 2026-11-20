import numpy as np, pickle
from scipy.optimize import brentq
from flint import arb, ctx
import bpf_exact as BX
cache={}
def Qb(N,z):
    if N not in cache: cache[N]=BX.build(N,BX.Rmat(N))
    ctx.prec=6*N+300; za=arb(str(z)) if False else arb(z); h=arb(10)**(-25)
    d=(BX.lndet(N,cache[N],za*(1+h))-BX.lndet(N,cache[N],za*(1-h)))/(2*h)
    return N/2-float(d.mid())/4
rows=[]
s=pickle.load(open('data/summary.pkl','rb'))
for r in s:
    if r[0] in (12,16,20) and 0<r[1]<=0.1: rows.append((int(r[0]),float(r[1]),r[0]/2-float(r[2]),'exact'))
for r in pickle.load(open('data/full_vs_boundary.pkl','rb')): rows.append((r[0],r[1],r[2],'MC'))
out=[]
print('  N     z        Q_un(full)   x=N*/N    g = z_eff/z    N/xi')
for N,z,Q,src in rows:
    if not (0.02*N < Q < 0.49*N): continue
    f=lambda lg: Qb(N,z*np.exp(lg))-Q
    try: lg=brentq(f,-1.5,4.5,xtol=1e-4)
    except Exception as e: print(N,z,'no root'); continue
    x=2*np.log2(1/z)/N; n_bulk=0.9*z**(8/7); xi=1/np.sqrt(3.04*max(n_bulk,1e-300))
    out.append((N,z,Q,x,np.exp(lg),N/xi,src)); print('%3d  %.0e   %7.3f     %.3f     %.3f        %.2f   %s'%(N,z,Q,x,np.exp(lg),N/xi,src),flush=True)
pickle.dump(out,open('data/geff.pkl','wb'))
