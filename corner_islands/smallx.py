import numpy as np, pickle, bpf_big as Bp
from scipy.special import gamma
A=8*np.log(2)/np.pi**2; Bc=4*np.sqrt(np.pi)*np.pi*np.e/(gamma(0.25)**2*np.log(2))
print('A = 8 ln2/pi^2 = %.5f   B = %.4f'%(A,Bc))
xs=(0.03,0.05,0.1,0.15,0.2,0.3,0.4,0.5); Ns=(24,32,40,48,56); out={}
for N in Ns:
    P=Bp.BoundaryPf(N)
    for x in xs:
        lnZ,n=P.logZ_and_n(np.full(4*N,2.0**(-x*N/2))); out[(N,x)]=(1-n.sum()/(2*N), n[:N].copy())
    print('N',N,'done',flush=True)
pickle.dump(out,open('az/smallx.pkl','wb'))
print(' x    sigma_N (N=32,40,48,56)           extrap(1/N,1/N^2)   asymptotic A x ln(B/x)')
for x in xs:
    s=np.array([out[(N,x)][0] for N in Ns]); iN=1/np.array(Ns,float)
    c=np.polyfit(iN[1:],s[1:],2)[-1]; c1=np.polyfit(iN[2:],s[2:],1)[-1]
    print('%.2f  '%x+' '.join('%.4f'%v for v in s[1:])+'    %.4f (lin %.4f)    %.4f'%(c,c1,A*x*np.log(Bc/x)))
