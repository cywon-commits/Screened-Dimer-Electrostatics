import numpy as np, pickle, vp_outer as VO, tf_solve as T
M=64; du=2/M; w=-1+du*(np.arange(M)+.5); d=w/2; h=du/2
K=T.Kinf(); D=d[:,None]-d[None,:]
L=np.log(np.abs(D)+np.eye(M)); np.fill_diagonal(L,np.log(h/2)-1)
W=L+(K(D)+np.log(2))/2
u=T.H(0.5+d)-np.log(2)/2
a1=pickle.load(open('az/vp_M64_e0.01.pkl','rb')); a2=pickle.load(open('az/vp_M64_e0.02.pkl','rb'))
for x in (0.7,0.5):
  rho=2*a1[x][1]-a2[x][1]; rho=np.clip(rho,0,1)
  lin=4*(u@rho)*h; quad=4*0.5*(rho@W@rho)*h*h
  print('x',x,'TF: linear coef %.4f quad coef %.4f'%(lin,quad))
  out={}
  for eps in (0.02,0.01):
    O=VO.Outer(M,eps); S0,_=O.S(np.zeros(M)); r=[]
    for lam in (0.25,0.5,0.75,1.0):
      S,_=O.S(lam*rho); r.append(S-S0)
    out[eps]=np.array(r)
  dS=2*out[0.01]-out[0.02]
  for lam,v in zip((0.25,0.5,0.75,1.0),dS): print('   lam %.2f  VP dS %.4f   TF %.4f'%(lam,v,lam*lin+lam**2*quad))
