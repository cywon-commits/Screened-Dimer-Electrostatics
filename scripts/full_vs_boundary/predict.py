"""D_inf(x) = (1/ln 2) d/dx [ 2 int_0^1 dg(rho_x(u)) du ] with dg from the one-sided strip (onesided.py)."""
import sys, pickle, numpy as np
import exact_sol as E
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
import onesided_tab as O   # run from the repository root: PYTHONPATH=src:scripts/full_vs_boundary
def sol(x):
    lj=brentq(lambda lj: E.solve(1+np.exp(lj))['x']-x,-12,18,xtol=1e-13); return E.solve(1+np.exp(lj))
def dg_table(W,rmax):
    r=np.concatenate([[0.0],np.arange(0.05,rmax+1e-9,0.01)]); v=np.full(len(r),np.nan); v[0]=0.0
    for w in sorted({k[0] for k in O.d if k[0]<=W and (k[0],'one') in O.d}):       # the largest available width wins
        new=O.dg(w,r[1:]); ok=~np.isnan(new); v[1:][ok]=new[ok]
    assert not np.isnan(v).any(),(W,r[np.isnan(v)])
    return CubicSpline(r,v,bc_type=((1,0.0),'not-a-knot'))
def I(s,dg):
    rho=np.clip((1-s['h'])/2,0,None); return 2*np.trapezoid(dg(rho),s['u'])
if __name__=='__main__':
    Ws=[int(a) for a in sys.argv[1].split(',')]; rmax=float(sys.argv[2]); h=0.01
    print('x       rho_M   '+'   '.join(f'W={W}' for W in Ws))
    for x in [float(a) for a in sys.argv[3].split(',')]:
        sp,sm,s0=sol(x+h),sol(x-h),sol(x)
        print(f'{x:.4f}  {(1-s0["h"][-1])/2:.3f}   '+'  '.join(f'{(I(sp,dg_table(W,rmax))-I(sm,dg_table(W,rmax)))/(2*h)/np.log(2):+.4f}' for W in Ws),flush=True)
