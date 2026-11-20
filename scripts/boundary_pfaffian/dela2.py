import numpy as np, twopair as T
from xkernel import lnR_exact
N=400
for a0 in (0.5,0.35,0.2,0.12):
    i=int(a0*N); r1=lnR_exact(N,i,i)
    out=[]
    for d in (2,5,10,20,40):
        i2=i+d
        if i2>=N: continue
        w=T.lnratio_det(N,[i,i2],[i,i2])-r1-lnR_exact(N,i2,i2)
        out.append((d, round(w,3), round(w-2*np.log(d/N),3)))
    print(a0,out)
