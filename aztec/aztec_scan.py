import sys, time, numpy as np, aztec as A, pickle
N=int(sys.argv[1]); zs=[float(x) for x in sys.argv[2:]]
for z in zs:
    t=time.time(); lnZ,m,l,r=A.Aztec(N,z).run()
    inc=A.bond_maps(N,l,r)
    phi=[];charge=0;mon=0
    data={}
    for k in range(2*N+1):
        for i in range(len(m[k])):
            c=A.cell_coords(N,k,i); f=max(inc.get(c,(0,0))[1],m[k][i]); phi.append(f)
            data[c]=(m[k][i]*(1 if k%2==0 else -1), m[k][i], f)
    phi=np.array(phi)
    pickle.dump(data,open(f'az/N{N}_z{z}.pkl','wb'))
    print(f"N {N} z {z} lnZ {lnZ:.6f} frozen(>0.99) {np.mean(phi>0.99):.4f} frozen(>0.9) {np.mean(phi>0.9):.4f} <n> {np.mean([v[1] for v in data.values()]):.5f} t {time.time()-t:.1f}s",flush=True)
