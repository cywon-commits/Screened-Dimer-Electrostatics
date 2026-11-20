import numpy as np, pickle, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import brentq
from exact_sol import solve
sm=pickle.load(open('data/smallx.pkl','rb')); pr=pickle.load(open('data/profiles.pkl','rb'))
fig,ax=plt.subplots(1,3,figsize=(15.5,4.3))
ls=np.linspace(-9,9.5,80); S=[solve(1+np.exp(l)) for l in ls]
X=np.array([s['x'] for s in S]); Sg=np.array([s['sigma'] for s in S]); Ds=np.array([s['dstar'] for s in S])
ax[0].plot(X,Sg,'k-',label='exact solution (hodograph, no parameters)')
Ns=(40,48,56)
for x in (0.1,0.15,0.2,0.3,0.4,0.5):
    v=np.array([sm[(N,x)][0] for N in Ns]); ax[0].plot(x,np.polyfit(1/np.array(Ns,float),v,1)[-1],'ro',mfc='none',label='boundary Pfaffian, N→∞' if x==0.1 else None)
for x,v in ((0.6,0.799),(0.7,0.867),(0.8,0.924),(0.9,0.970),(1.0,1.0)): ax[0].plot(x,v,'ro',mfc='none')
from scipy.special import gamma
A=8*np.log(2)/np.pi**2; Bc=4*np.sqrt(np.pi)*np.pi*np.e/(gamma(0.25)**2*np.log(2)); xx=np.linspace(0.03,1,100)
ax[0].plot(xx,A*xx*np.log(Bc/xx),'b:',label='small-x law')
ax[0].set_xlabel('x = N*/N'); ax[0].set_ylabel('σ∞'); ax[0].legend(fontsize=8)
for x,c in ((0.2,'C1'),(0.3,'C2'),(0.5,'C3'),(0.7,'C4')):
    lj=brentq(lambda l: solve(1+np.exp(l))['x']-x,-9,11,xtol=1e-10); s=solve(1+np.exp(lj))
    u=s['u']; rho=(1-s['h'])/2
    d=np.concatenate([[0],u,(2-u)[::-1],[2]]); r=np.concatenate([[0],rho,rho[::-1],[0]])
    ax[1].plot(d/2-0.5,r,'-',color=c,label=f'exact solution x={x}')
    src=sm if (56,x) in sm else None
    n56=(sm[(56,x)][1] if src else pr[(56,x)][1][:56]); n40=(sm[(40,x)][1] if src else pr[(40,x)][1][:40])
    dl=(np.arange(56)+.5)/56-.5; p40=np.interp(dl,(np.arange(40)+.5)/40-.5,n40)
    ax[1].plot(dl,(56*n56-40*p40)/16,'o',color=c,ms=3)
ax[1].set_xlabel('δ (position along edge)'); ax[1].set_ylabel('vacancy density ρ(δ)'); ax[1].set_title('dots: boundary Pfaffian (N=40,56 extrapolated)',fontsize=9); ax[1].legend(fontsize=8)
ax[2].plot(X,Ds,'k-',label='frozen length per corner d* (half-edge = 1)'); ax[2].plot(xx,4*xx*np.log(2)/np.pi,'b:',label='isolated corner: 4x ln2/π')
ax[2].set_xlabel('x = N*/N'); ax[2].set_ylabel('d*'); ax[2].legend(fontsize=8)
plt.tight_layout(); plt.savefig('figures/working/az_exact.png',dpi=130)
