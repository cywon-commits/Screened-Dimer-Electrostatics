import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from img_test import prof
W=20; z=0.002
def fitalpha(r,W,skip=1):
    def res(p):
        a,de,do,Ae,Ao=p; out=[]
        for x in range(skip,W//2):
            dd=de if x%2==0 else do
            out.append(np.log(r[x])-np.log(Ae if x%2==0 else Ao)-a*np.log(np.sin(np.pi*(x+dd)/(W+2*dd))))
        return out
    return least_squares(res,[0,0.5,0.5,1,1],bounds=([-2,-0.9,-0.9,0.01,0.01],[2,5,5,10,10])).x[0]
fig,ax=plt.subplots(1,3,figsize=(14,4))
for zb,c in ((0.002,'C0'),(0.02,'C1'),(1.0,'C3'),(10.0,'C2')):
    zs=np.full(W,z); zs[0]=zs[-1]=zb; n=prof(W,zs); r=n/n[W//2]
    d=np.arange(1,W//2); env=[(r[k]+r[k+1])/2 for k in range(1,W//2-1)]
    ax[0].plot(np.arange(1,W//2-1)+0.5,env,'o-',color=c,label=f'z_b={zb}')
x=np.linspace(1,10,50); ax[0].plot(x,(np.sin(np.pi*x/W))**0.25,'k--',lw=1,label='sin(πd/W)^(+1/4)')
ax[0].plot(x,(np.sin(np.pi*x/W))**-0.25,'k:',lw=1,label='sin(πd/W)^(−1/4)')
ax[0].set_xlabel('distance from edge d'); ax[0].set_ylabel('n(d)/n(center), even-odd averaged'); ax[0].set_title(f'bulk z={z}, W={W}',fontsize=9); ax[0].legend(fontsize=7)
zbs=np.geomspace(0.002,30,22); al=[]
for zb in zbs:
    zs=np.full(W,z); zs[0]=zs[-1]=zb; n=prof(W,zs); al.append(fitalpha(n/n[W//2],W))
ax[1].semilogx(zbs,al,'o-'); ax[1].axhline(0.25,color='k',ls='--',lw=.8); ax[1].axhline(-0.25,color='k',ls=':',lw=.8)
ax[1].set_ylim(-0.3,0.32); ax[1].text(0.05,0.27,'Dirichlet (same-sign image)',fontsize=8); ax[1].text(0.05,-0.285,'Neumann (opposite-sign image)',fontsize=8)
ax[1].set_xlabel('edge-site fugacity z_b (bulk z=0.002)'); ax[1].set_ylabel('near-edge exponent α'); ax[1].set_title('W=20',fontsize=9)
Ws=np.array([10,12,14,16,18,20]); zst=np.array([0.1715,0.1740,0.1359,0.1383,0.1231,0.1307])
ax[2].loglog(Ws,zst,'o'); ax[2].loglog(Ws,np.sqrt(0.3/Ws),'k--',label='z_b* = (0.3/W)^(1/2)')
ax[2].set_xlabel('strip width W'); ax[2].set_ylabel('crossover z_b* (α=0)'); ax[2].legend(fontsize=8)
plt.tight_layout(); plt.savefig('az_image_sign.png',dpi=130)
