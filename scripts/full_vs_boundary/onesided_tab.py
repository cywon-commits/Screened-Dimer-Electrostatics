import pickle, numpy as np, sys
from scipy.interpolate import CubicSpline
d=pickle.load(open('data/full_vs_boundary/onesided.pkl','rb'))
rs=np.array([0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.425,0.45])
def phi(c):
    lz,F,rho,ph=c
    bad=[k for k in range(1,len(rho)) if not (rho[k]>rho[k-1] and F[k]>F[k-1])]      # wrong eigenvalue at the smallest z
    k0=(max(bad)+3) if bad else 2
    keep=list(range(k0,len(rho)-2))
    return CubicSpline(rho[keep],ph[keep]),rho[keep[0]],rho[keep[-1]]
def dg(W,r):
    fb,lo1,hi1=phi(d[(W,'bdry')]); fo,lo2,hi2=phi(d[(W,'one')])
    return np.where((r>=max(lo1,lo2))&(r<=min(hi1,hi2)),fo(r)-fb(r),np.nan)
if __name__=='__main__':
    print('W   '+'  '.join(f'rho={r:.3f}' for r in rs))
    for W in sorted({k[0] for k in d if (k[0],'one') in d and (k[0],'bdry') in d}):
        print(f'{W:2d}  '+'  '.join(f'{v:9.6f}' for v in dg(W,rs)))
