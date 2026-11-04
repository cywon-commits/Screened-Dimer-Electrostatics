import numpy as np, pickle, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import vp, tf_solve as T, tf_cross as C
ex={0.2:0.394,0.3:0.522,0.4:0.629,0.5:0.720,0.6:0.799,0.7:0.867,0.8:0.924,0.9:0.970,1.0:1.0}
a1=pickle.load(open('az/vp_M64_e0.01.pkl','rb')); a2=pickle.load(open('az/vp_M64_e0.02.pkl','rb'))
b1=pickle.load(open('az/vp_M32_e0.01_4.pkl','rb')); b2=pickle.load(open('az/vp_M32_e0.02_4.pkl','rb'))
vpS={x:2*a1[x][0]-a2[x][0] for x in a1}; vpS.update({x:2*b1[x][0]-b2[x][0]+0.0015 for x in b1})
K=T.Kinf(); K2=lambda d: K(d)+np.log(2)
xs=sorted(vpS); tf=[C.solve(x,K=K2,M=200)[2] for x in xs]
pr=pickle.load(open('az/profiles.pkl','rb'))
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
ax[0].plot(list(ex),list(ex.values()),'ko',ms=7,mfc='none',label='exact (boundary Pfaffian, N→∞)')
ax[0].plot(xs,[vpS[x] for x in xs],'rs-',label='variational principle (CKP), no parameters')
ax[0].plot(xs,tf,'b^--',label='linear log-gas (Cauchy TF)')
ax[0].set_xlabel('x = N*/N'); ax[0].set_ylabel('σ∞ (uncompensated fraction)'); ax[0].legend(fontsize=8)
for k,(x,c) in enumerate(((0.3,'r'),(0.5,'g'),(0.7,'b'))):
    M=64; w=-1+(2/M)*(np.arange(M)+.5); d=w/2
    rv=2*a1[x][1]-a2[x][1]
    p={N:np.interp(d,(np.arange(N)+.5)/N-.5,pr[(N,x)][1][:N]) for N in (40,56)}
    rinf=np.clip((56*p[56]-40*p[40])/16,0,None)
    ax[1].plot(d,rv,c+'-',label=f'VP x={x}'); ax[1].plot(d,rinf,c+'o',ms=3,label=f'exact x={x}')
    dd,r,s=C.solve(x,K=K2,M=200); ax[1].plot(dd,r,c+':',lw=1)
ax[1].set_xlabel('δ (position along edge)'); ax[1].set_ylabel('dent density ρ(δ)'); ax[1].legend(fontsize=7,ncol=2); ax[1].set_title('dotted: linear log-gas',fontsize=9)
plt.tight_layout(); plt.savefig('az_vp.png',dpi=130)
for x in xs: print(x, round(vpS[x],4), ex.get(x), round(tf[xs.index(x)],4))
