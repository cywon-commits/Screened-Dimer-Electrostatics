import numpy as np, pickle, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.special import gamma
out=pickle.load(open('az/smallx.pkl','rb'))
A=8*np.log(2)/np.pi**2; Bc=4*np.sqrt(np.pi)*np.pi*np.e/(gamma(0.25)**2*np.log(2))
# far field along the edge: charges eps=(-1)^(a+b) at (2a,2b) (corner S at origin, W at (2,0)); field d/du of sum eps ln|r-r_i|
n=np.arange(-400,401); Ag,Bg=np.meshgrid(n,n,indexing='ij'); sg=(-1.0)**(Ag+Bg); w=np.exp(-2e-4*(Ag**2+Bg**2))
def F(d): return np.array([(sg*w*(di-2*Ag)/((di-2*Ag)**2+(2*Bg)**2)).sum() for di in np.atleast_1d(d)])
def rho_pred(d,x):
    C=8*x*np.log(2)/np.pi; ds=C/2
    near=lambda t: np.where(t>ds,(2/np.pi)*np.arcsin(np.clip(ds/np.maximum(t,1e-12),0,1)),1.0)
    h=near(d)+near(2-d)+(C/np.pi)*(F(d)-1/d-1/(2-d))
    return np.clip((1-h)/2,0,None)
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
xs=np.linspace(0.02,1,200); ax[0].plot(xs,A*xs*np.log(Bc/xs),'k-',label=r'$\sigma_\infty=\frac{8\ln 2}{\pi^2}\,x\ln(B/x)$, B=6.645 (no parameters)')
Ns=(40,48,56)
for x in (0.05,0.1,0.15,0.2,0.3,0.4,0.5):
    s=np.array([out[(N,x)][0] for N in Ns]); ax[0].plot(x,np.polyfit(1/np.array(Ns,float),s,1)[-1],'ro',mfc='none',label='exact (boundary Pfaffian, N→∞)' if x==0.05 else None)
for x,v in ((0.6,0.799),(0.7,0.867),(0.8,0.924),(0.9,0.970),(1.0,1.0)): ax[0].plot(x,v,'ro',mfc='none')
ax[0].set_xlabel('x = N*/N'); ax[0].set_ylabel('σ∞'); ax[0].legend(fontsize=8)
print(' x   support edge d* (theory 4x ln2/pi, half-edge units) ; centre rho: theory, exact N=56, extrap')
for x,c in ((0.1,'C0'),(0.2,'C1'),(0.3,'C2'),(0.5,'C3')):
    d=np.linspace(0.01,1.99,300); rp=rho_pred(d,x); ax[1].plot(d/2-0.5,rp,'-',color=c,label=f'theory x={x}')
    N=56; ne=out[(56,x)][1]; n40=out[(40,x)][1]; dl=(np.arange(N)+.5)/N-.5
    p40=np.interp(dl,(np.arange(40)+.5)/40-.5,n40); ext=(56*ne-40*p40)/16
    ax[1].plot(dl,ext,'o',color=c,ms=3)
    print('%.1f  d*=%.3f   centre %.4f  %.4f  %.4f'%(x,4*x*np.log(2)/np.pi,rho_pred(np.array([1.0]),x)[0],ne[27:29].mean(),ext[27:29].mean()))
ax[1].set_xlabel('δ (position along edge)'); ax[1].set_ylabel('vacancy density ρ(δ)'); ax[1].set_title('lines: island + lattice far field; dots: exact (N=40,56 extrapolated)',fontsize=9); ax[1].legend(fontsize=8)
plt.tight_layout(); plt.savefig('az_smallx.png',dpi=130)
