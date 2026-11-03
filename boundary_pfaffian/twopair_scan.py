import numpy as np, twopair as T, sys, pickle
from math import comb, log
res={}
for N in [int(a) for a in sys.argv[1:]]:
    c=N//2; ds=sorted(set([2,4,8,16,32,64,128,256,512]+[int(N*f) for f in (0.05,0.1,0.15,0.2,0.3)]))
    ds=[d for d in ds if d<=0.6*N]
    need=set()
    for d in ds: need|={c-d//2,c-d//2+d}
    B=T.binom_rows(N,need)
    row=[]
    for d in ds:
        I=[c-d//2,c-d//2+d]; J=I[:]
        lnd=T.lnratio_det(N,I,J,B)
        l11=T.lnbig(T.S(N,I[0],J[0],B))-N*log(2); l22=T.lnbig(T.S(N,I[1],J[1],B))-N*log(2)
        l12=T.lnbig(T.S(N,I[0],J[1],B))-N*log(2); l21=T.lnbig(T.S(N,I[1],J[0],B))-N*log(2)
        row.append((d,lnd,l11,l22,l12,l21))
        print(N,d,'F=%.4f'%(lnd-l11-l22),flush=True)
    res[N]=row
pickle.dump(res,open('az/twopair_%s.pkl'%'_'.join(sys.argv[1:]),'wb'))
