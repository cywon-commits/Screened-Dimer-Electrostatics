"""Dilute limit (one vacancy pair): total pair weight of the full model over that of the boundary model,
sum_{a,b} M(A\\{a,b})/M(A) = sum |K^-1(a,b)|, for the Aztec diamond of order N."""
import numpy as np, sys
def aztec(N):
    cells=[(i,j) for i in range(-N,N) for j in range(-N,N) if abs(i+.5)+abs(j+.5)<=N]
    return cells
def run(N):
    C=aztec(N); B=[c for c in C if (c[0]+c[1])%2==0]; W=[c for c in C if (c[0]+c[1])%2]; wi={c:k for k,c in enumerate(W)}
    K=np.zeros((len(B),len(W)))
    for k,(i,j) in enumerate(B):
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            n=(i+d[0],j+d[1])
            if n in wi: K[k,wi[n]]=1.0 if d[1]==0 else (-1.0)**i
    Ki=np.abs(np.linalg.inv(K))           # [white, black]
    depth=lambda c: N-(abs(c[0]+.5)+abs(c[1]+.5))        # 0 on the outermost cells
    db=np.array([depth(c) for c in B]); dw=np.array([depth(c) for c in W])
    Sfull=Ki.sum(); Sb=Ki[np.ix_(dw==0,db==0)].sum()
    # depth profile: one vacancy at depth d, the other anywhere / on the edge
    prof=[Ki[np.ix_(dw==0,db==d)].sum()/Sb for d in range(0,8)]
    one_free=Ki[np.ix_(dw==0,np.ones(len(B),bool))].sum()/Sb       # black anywhere, white on the edge
    return Sfull/Sb, one_free, prof
for N in (4,6,8,10,12,14,16,18,20,22,24,26):
    r,one,prof=run(N)
    print(f'N={N:2d}  S_full/S_bdry={r:8.4f}  sqrt={np.sqrt(r):.4f}  one-vacancy factor={one:.4f}  (one)^2={one*one:.4f}  depth profile {np.round(prof,4)}',flush=True)
