"""Exact solution of the corner-charge problem for the Aztec diamond with boundary vacancies (all x).
Hodograph plane: zeta = (z-i)/(z+i) in the sector -pi/4 < arg < 0, J = (rho^4 + rho^-4)/2 on the real side.
One parameter J_M > 1 (image of the edge midpoint M).  F(phi) = rho+1/rho = 2 cosh(arccosh(s)/4), s = c0 - c1 cos(phi).
  support part SM (v=0):  u(phi0) = 1 + (rho+1/rho) sum a_n sin(n phi0) / (2 sqrt(c1) sin(phi0/2))
  antimirror MO (u=1):    v(J)    = 1 - sum a_n eps^n / (rho+1/rho),  eps = sqrt(B^2-1)-B, B = (J-c0)/c1
  x ln2 = int_MO ln t dv, t=(1+rho)/(1-rho);   sigma = int_0^1 h du, h = (4/pi) arctan(rho) on SM, 1 on the frozen part.
"""
import numpy as np
def cheb(JM, n=3000):
    c0=(1+JM)/2; c1=(JM-1)/2; k=np.arange(n); phi=np.pi*(k+.5)/n
    F=2*np.cosh(np.arccosh(c0-c1*np.cos(phi))/4)
    from scipy.fft import dct
    a=dct(F,type=2)/n; a[0]/=2
    keep=np.abs(a)>1e-15*np.abs(a).max(); last=np.max(np.where(keep)[0])+1
    return c0,c1,a[:last]
def rho_of_J(J): return (J-np.sqrt(J*J-1))**0.25
def solve(JM, m=4000):
    c0,c1,a=cheb(JM); n=np.arange(len(a))
    # SM: phi0 in (0,pi)
    p=np.pi*(np.arange(m)+.5)/m; J=c0-c1*np.cos(p); r=rho_of_J(J)
    S=(a[None,:]*np.sin(np.outer(p,n))).sum(1)
    u=1+(r+1/r)*S/(2*np.sqrt(c1)*np.sin(p/2))
    dstar=1+2*(n*a).sum()/np.sqrt(c1)
    h=(4/np.pi)*np.arctan(r)
    uu=np.concatenate([[dstar],u,[1.0]]); hh=np.concatenate([[1.0],h,[(4/np.pi)*np.arctan(rho_of_J(JM))]])
    sigma=dstar+np.trapezoid(hh,uu)
    # MO: rho in (0,rho_M)
    rM=rho_of_J(JM); q=(np.arange(m)+.5)/m; rr=rM*np.sin(np.pi*q/2)**2   # cluster near 0 and rho_M
    rr=np.concatenate([[1e-9],rr,[rM]]); JJ=(rr**4+rr**-4)/2; B=np.maximum((JJ-c0)/c1,1.0); e=-1/(B+np.sqrt(B*B-1))
    v=1-(a[None,:]*e[:,None]**n[None,:]).sum(1)/(rr+1/rr)
    lt=np.log((1+rr)/(1-rr)); x=-np.trapezoid(lt,v)/np.log(2)
    return dict(x=x,sigma=sigma,dstar=dstar,u=uu,h=hh,v=v,rho_MO=rr)
if __name__=='__main__':
    from scipy.special import gamma
    A=8*np.log(2)/np.pi**2; Bc=4*np.sqrt(np.pi)*np.pi*np.e/(gamma(0.25)**2*np.log(2))
    print('   J_M        x       sigma     d*      | small-x law A x ln(B/x)   d* law 4x ln2/pi')
    for JM in (1.0001,1.01,1.1,1.5,3,10,50,300,3000,1e5,1e7):
        s=solve(JM); x=s['x']
        print('%9.4g  %.5f  %.5f  %.5f  |  %.5f   %.5f'%(JM,x,s['sigma'],s['dstar'],A*x*np.log(Bc/x),4*x*np.log(2)/np.pi))
