import numpy as np, sys, time, pickle, bpf_big as B
def tri_top(N,c):
    x,y=c; k=x+y-(N-1); b=x-y; u=k-N
    if abs(u)>abs(b): return 1.0 if u<0 else 0.0
    if abs(b)>abs(u): return 0.0
    return 0.5 if u<0 else 0.0
out={}
for N in [int(a) for a in sys.argv[1:]]:
    t=time.time(); P=B.BoundaryPf(N)
    wt=np.array([tri_top(N,c) for c in P.bd]); sg=np.array([1.0 if (c[0]+c[1]-(N-1))%2==0 else -1.0 for c in P.bd])
    res=[]
    for z in np.geomspace(1e-14,0.1,79):
        lnZ,n=P.logZ_and_n(np.full(len(P.bd),z))
        Q=np.sum(wt*sg*n); res.append((z,Q,lnZ,n.sum()))
    out[N]=res; print(N,'done %.0fs'%(time.time()-t),flush=True)
    pickle.dump(out,open('az/bscan_%s.pkl'%'_'.join(sys.argv[1:]),'wb'))
