import numpy as np, pickle, vp_t as vp, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from ckp import ent
M=64; du=2/M; wg=-1+du*(np.arange(M)+.5)
a1=pickle.load(open('az/vp_M64_e0.01.pkl','rb')); a2=pickle.load(open('az/vp_M64_e0.02.pkl','rb'))
fig,ax=plt.subplots(1,4,figsize=(16,4.2)); out={}
for k,x in enumerate((None,0.7,0.5,0.3)):
    rho=np.zeros(M) if x is None else np.clip(2*a1[x][1]-a2[x][1],0,1)
    S,V,_=vp.solve_ml(lambda w: np.interp(w,wg,rho),Ms=(16,32,64,128),eps=0.0)
    sl=V.slopes(V.H); E=np.mean([ent(np.clip(a,-1,1)-np.clip(b,-1,1),np.clip(a,-1,1)+np.clip(b,-1,1)) for a,b in sl],0)
    Mm=V.M; g=-1+(2/Mm)*(np.arange(Mm)+.5); U,Vv=np.meshgrid(g,g,indexing='ij'); X=(U-Vv)/2; Y=(U+Vv)/2
    out[x]=(X,Y,E)
    c=ax[k].tricontourf(X.ravel(),Y.ravel(),E.ravel()/0.5831,levels=np.linspace(0,1,21),cmap='magma')
    ax[k].tricontour(X.ravel(),Y.ravel(),E.ravel(),levels=[0.01],colors='c',linewidths=1.2)
    t=np.linspace(0,2*np.pi,200); ax[k].plot(np.cos(t)/np.sqrt(2),np.sin(t)/np.sqrt(2),'w:',lw=1)
    ax[k].set_aspect('equal'); ax[k].set_title(('plain AD' if x is None else f'x={x}, σ∞={1-rho.sum()*du:.3f}')+'\nlocal entropy / max',fontsize=9)
    ax[k].set_xticks([]); ax[k].set_yticks([])
    print(x,'frozen area fraction %.3f'%(E<0.01).mean(),flush=True)
plt.tight_layout(); plt.savefig('az_vp_maps.png',dpi=120); pickle.dump(out,open('az/vp_maps.pkl','wb'))
