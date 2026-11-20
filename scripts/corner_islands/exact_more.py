import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import brentq
from exact_sol import solve, cheb
L2=np.log(2)
# (a) x -> 1 asymptotics: compare with log-gas closed form  1-sigma = 2 ln2 (1-x) / [1 - 4 ln2 + ln(16/(1-sigma))]
print('x->1:   J_M-1      1-x        1-sigma     log-gas formula for 1-sigma   ratio   | c1/16')
for e in (1e-1,1e-2,1e-3,1e-4,1e-5):
    s=solve(1+e,m=20000); x=s['x']; t=1-s['sigma']
    g=lambda y: y-2*L2*(1-x)/(1-4*L2+np.log(16/y))
    try: y=brentq(g,1e-14,0.5)
    except Exception: y=np.nan
    print('     %8.1e  %.3e  %.4e  %.4e   %.4f   | %.4e'%(e,1-x,t,y,t/y,e/32))
# (b) frozen-island boundary (arctic curve): envelope of v + T u = Psi(T), T = tan|phi| on the arc
def arctic(JM, m=600):
    c0,c1,a=cheb(JM); n=np.arange(len(a))
    T=np.linspace(1e-4,1-1e-4,m); J=1-8*T**2/(1+T**2)**2
    Bp=(c0-J)/c1; ep=Bp-np.sqrt(Bp*Bp-1)
    Th=-2**-0.25*(1-J)**-0.25*(a[None,:]*ep[:,None]**n[None,:]).sum(1)
    Psi=1+T+np.sqrt(T)*Th
    dP=np.gradient(Psi,T); return dP, Psi-T*dP      # (u,v)
fig,ax=plt.subplots(1,1,figsize=(5.2,5.2))
for x,c in ((0.999,'k'),(0.9,'C4'),(0.7,'C3'),(0.5,'C2'),(0.3,'C1'),(0.15,'C0')):
    lj=brentq(lambda l: solve(1+np.exp(l))['x']-x,-12,11,xtol=1e-10); JM=1+np.exp(lj); u,v=arctic(JM)
    if x==0.999: print('x~1 check: (u-1)^2+(v-1)^2 on the arc: min %.4f max %.4f (arctic circle = 1)'%(((u-1)**2+(v-1)**2)[5:-5].min(),((u-1)**2+(v-1)**2)[5:-5].max()))
    # corner S at (0,0); octant 0<=v<=u; reflect to the quarter and to all four corners (edges along u and v axes, square [0,2]^2)
    U=np.concatenate([u,v[::-1]]); V=np.concatenate([v,u[::-1]])
    for (a_,b_) in ((U,V),(2-U,V),(U,2-V),(2-U,2-V)):
        ax.plot(a_,b_,'-',color=c,label=f'x={x}' if (a_ is U and b_ is V) else None)
ax.plot([0,2,2,0,0],[0,0,2,2,0],'k-',lw=.6); ax.set_aspect('equal'); ax.legend(fontsize=8,loc='center'); ax.set_xticks([]); ax.set_yticks([])
ax.set_title('frozen-corner boundaries from the exact solution\n(Aztec diamond drawn as a square, edges = sides)',fontsize=9)
plt.tight_layout(); plt.savefig('figures/working/az_arctic.png',dpi=130)
