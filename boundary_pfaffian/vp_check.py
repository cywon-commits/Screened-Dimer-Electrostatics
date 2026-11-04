import numpy as np, pickle, vp_t as vp
L2=np.log(2); M=64; du=2/M; wg=-1+du*(np.arange(M)+.5)
a1=pickle.load(open('az/vp_M64_e0.01.pkl','rb')); a2=pickle.load(open('az/vp_M64_e0.02.pkl','rb'))
b=pickle.load(open('az/vp2_M64_e0.0_0.5,0.2.pkl','rb'))
x=0.5
cands={'old(eps-extrap)':np.clip(2*a1[x][1]-a2[x][1],0,1),'new(eps=0)':b[x][1]}
for lam in (0.95,1.05): cands['old x%.2f'%lam]=np.clip(lam*cands['old(eps-extrap)'],0,1)
for name,rho in cands.items():
    S,V,_=vp.solve_ml(lambda w: np.interp(w,wg,rho),Ms=(16,32,64,128),eps=0.0)
    phi=S-2*x*L2*rho.sum()*du/2
    print(name,'sigma %.4f'%(1-rho.sum()*du),'Phi %.6f'%phi,flush=True)
