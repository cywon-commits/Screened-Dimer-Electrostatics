import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import tm_diag_b as T
from loz_cap import lnZ
from ckp import L
fig,ax=plt.subplots(1,2,figsize=(10.5,4.2))
ts=np.linspace(0.02,1,100); ax[0].plot(ts,(2/np.pi)*np.arctan(ts),'k-',label=r'$\rho=(2/\pi)\arctan(z^{2/W})$')
T.RIN=0.0
for W,m in ((6,'v'),(10,'s'),(14,'o')):
    tt=[t for t in (0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0) if t**(W/2)>2e-4]
    ax[0].plot(tt,[T.densities(W,t**(W/2))[2][0] for t in tt],m,mfc='none',label=f'exact TM, W={W}')
ax[0].set_xlabel(r'$z^{2/W}=2^{-N^*/W}$'); ax[0].set_ylabel('edge vacancy density ρ'); ax[0].set_title('domino staircase strip (edge-type frozen state): no threshold',fontsize=9); ax[0].legend(fontsize=8)
s=lambda r:(2*L(np.pi*r/2)+L(np.pi*(1-r)))/np.pi
xs=np.linspace(0,1,200); r=(2/np.pi)*np.arccos(2**(xs-1)); ax[1].plot(xs,s(r)-r*xs*np.log(2),'k-',label=r'$s(\rho^*)-\rho^* x\ln 2$,  $\rho^*=(2/\pi)\arccos 2^{x-1}$')
M=12
for x in (0.1,0.2,0.35,0.5,0.65,0.8,0.95):
    f={N:lnZ(M,N,2.0**(-x*N/2))/(M*N) for N in (16,32,64)}
    ax[1].plot(x,f[32],'s',color='C0',mfc='none',label='exact TM, N=32' if x==0.1 else None)
    ax[1].plot(x,f[64],'o',color='C1',mfc='none',label='N=64' if x==0.1 else None)
    ax[1].plot(x,2*f[64]-f[32],'x',color='C3',label='1/N extrapolation' if x==0.1 else None)
ax[1].set_xlabel(r'$x=N^*/N$'); ax[1].set_ylabel('free energy per site'); ax[1].set_title('lozenge capacitor (vertex-type frozen state): threshold x=1, √ onset',fontsize=9); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('figures/working/az_isotherms.png',dpi=130)
