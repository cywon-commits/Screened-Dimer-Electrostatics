import numpy as np, pickle, vp_outer as VO, tf_solve as T, sys
M=64; du=2/M; w=-1+du*(np.arange(M)+.5); d=w/2; h=du/2
K=T.Kinf(); D=d[:,None]-d[None,:]
L=np.log(np.abs(D)+np.eye(M)); np.fill_diagonal(L,np.log(h/2)-1)
W=L+(K(D)+np.log(2))/2
a1=pickle.load(open('az/vp_M64_e0.01.pkl','rb')); a2=pickle.load(open('az/vp_M64_e0.02.pkl','rb'))
for x in (0.9,0.7,0.5,0.3):
  if x in a1: rho=np.clip(2*a1[x][1]-a2[x][1],0,1)
  else: continue
  tf=4*(rho@W@rho)*h*h
  cs={}
  for eps in (0.02,0.01):
    O=VO.Outer(M,eps); vals=[]
    for lam in (0.9,1.0,1.1):
      S,_=O.S(lam*rho); vals.append(S)
    cs[eps]=(vals[0]-2*vals[1]+vals[2])/0.01
  print('x',x,'d2S/dlam2  VP eps.02 %.4f eps.01 %.4f extrap %.4f | linear log-gas %.4f'%(cs[0.02],cs[0.01],2*cs[0.01]-cs[0.02],tf),flush=True)
