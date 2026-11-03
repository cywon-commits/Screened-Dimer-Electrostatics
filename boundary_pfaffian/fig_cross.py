import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from xkernel import E, lnR_exact, H
import twopair as T
fig,ax=plt.subplots(1,2,figsize=(10,4))
a=np.linspace(0.02,0.5,60); ac=1-1/np.sqrt(2)
ax[0].plot(a,[E(t,t)[0] for t in a],'k-',label='saddle E(a,a)')
ax[0].plot(a,2*H(a)-np.log(2),'b--',label='crossing 2H(a)-ln2')
ax[0].plot(a,2*a*np.log(1+np.sqrt(2)),'r:',label='Delannoy 2a ln(1+√2)')
for N,m in ((400,'o'),(1600,'s')):
    pts=a[::6]; ax[0].plot(pts,[lnR_exact(N,int(t*N),int(t*N))/N for t in pts],m,mfc='none',label=f'exact lnR/N, N={N}')
ax[0].axvline(ac,color='gray',lw=.5); ax[0].set_xlabel('a (dent position from shared corner)'); ax[0].set_ylabel('(1/N) ln R'); ax[0].legend(fontsize=7); ax[0].set_ylim(-0.1,0.75)
N=400
for a0,c in ((0.5,'b'),(0.35,'g'),(0.2,'r'),(0.12,'m')):
    i=int(a0*N); r1=lnR_exact(N,i,i); ds=[1,2,3,5,8,12,20,30,40,60]
    w=[T.lnratio_det(N,[i,i+d],[i,i+d])-r1-lnR_exact(N,i+d,i+d) for d in ds]
    ax[1].semilogx(ds,w,'o-',color=c,label=f'a={a0}')
ax[1].semilogx(ds,2*np.log(np.array(ds)/N)+3*np.log(2),'k:',label='Cauchy 2ln(d/N)+K(0)')
ax[1].set_xlabel('separation d (sites), N=400'); ax[1].set_ylabel('two-pair interaction'); ax[1].legend(fontsize=7); ax[1].set_ylim(-12,1)
plt.tight_layout(); plt.savefig('az_crossover.png',dpi=130)
