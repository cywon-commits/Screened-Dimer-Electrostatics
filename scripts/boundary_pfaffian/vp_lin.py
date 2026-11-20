import numpy as np, vp, vp_outer as VO
from xkernel import E
M=64; du=2/M; w=-1+du*(np.arange(M)+.5); dl=w/2
res={}
for eps in (0.02,0.01,0.005):
    O=VO.Outer(M,eps); S0,g=O.S(np.zeros(M)); res[eps]=g/(2*du)
u=2*res[0.01]-res[0.02]; u2=2*res[0.005]-res[0.01]
for k in range(0,M//2,3):
    a=0.5-abs(dl[k]); print('delta %.3f  a %.3f  VP u: %.4f (eps.01/.005 extrap %.4f)   E/2 %.4f  crossing H-ln2/2 %.4f'%(dl[k],a,u[k],u2[k],E(a,a)[0]/2, -(a*np.log(a)+(1-a)*np.log(1-a))-np.log(2)/2))
