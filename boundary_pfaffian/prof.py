import numpy as np, pickle, bpf_big as B
out={}
for N in (40,56):
    P=B.BoundaryPf(N)
    for x in (0.3,0.4,0.5,0.7):
        z=2.0**(-x*N/2); lnZ,n=P.logZ_and_n(np.full(4*N,z))
        out[(N,x)]=(np.array(P.bd),n)
        print(N,x,'tot n/N',round(n.sum()/N,3),flush=True)
pickle.dump(out,open('az/profiles.pkl','wb'))
