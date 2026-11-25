"""Exact solution in complete elliptic integrals (parameter p_M = rho_M^2 + rho_M^-2 > 2, k^2 = (p_M-2)/(2 p_M)).
With P = rho^2 + rho^-2 and  I(P) = (2/sqrt(2 p_M)) [ Pi(n1,k)/(p_M-P) + Pi(n2,k)/(p_M+P) ],  n1=(p_M-2)/(p_M-P), n2=(p_M-2)/(p_M+P)
(principal value for n1>1:  Pi(n,k) = K(k) - Pi(k^2/n, k)):
   edge SM (2<P<p_M):   u = 1 - (1/pi) sqrt((P+2)(p_M^2-P^2)) I(P),   vacancy density = 1/2 - (2/pi) arctan(rho)
   line MO (P>p_M):     v = 1 + (1/pi) sqrt((P-2)(P^2-p_M^2)) I(P)
   d* = u(P=2),   x ln2 = int_0^{rho_M} v 2 drho/(1-rho^2),   sigma = h_M + int_{h_M}^1 u dh,  h = (4/pi) arctan(rho)."""
import numpy as np
from scipy.special import elliprf, elliprj, ellipk, ellipe
from scipy.integrate import quad
def Pi(n,m):                       # complete third kind, parameter m = k^2, principal value for n>1
    if n>1: return ellipk(m)-Pi(m/n,m)
    return elliprf(0,1-m,1)+(n/3)*elliprj(0,1-m,1,1-n)
class Sol:
    def __init__(s,rM):
        s.rM=rM; s.pM=rM**2+rM**-2; s.m=(s.pM-2)/(2*s.pM); s.c=2/np.sqrt(2*s.pM)
    def I(s,P):
        pM=s.pM; return s.c*(Pi((pM-2)/(pM-P),s.m)/(pM-P)+Pi((pM-2)/(pM+P),s.m)/(pM+P))
    def u(s,rho):
        P=rho**2+rho**-2; return 1-np.sqrt((P+2)*(s.pM**2-P**2))*s.I(P)/np.pi
    def v(s,rho):
        P=rho**2+rho**-2; return 1+np.sqrt((P-2)*(P**2-s.pM**2))*s.I(P)/np.pi
    def dstar(s):
        pM,m=s.pM,s.m
        return 1-(2/np.pi)*np.sqrt(pM**2-4)*s.c*((ellipk(m)-ellipe(m)/(1-m))/(pM-2)+Pi((pM-2)/(pM+2),m)/(pM+2))
    def x(s):
        f=lambda r: s.v(r)*2/(1-r*r)
        return quad(f,0,s.rM,epsabs=1e-12,epsrel=1e-12,limit=400)[0]/np.log(2)
    def sigma(s):
        hM=4/np.pi*np.arctan(s.rM)
        f=lambda h: s.u(np.tan(np.pi*h/4))
        return hM+quad(f,hM,1,epsabs=1e-12,epsrel=1e-12,limit=400)[0]
if __name__=='__main__':
    from exact_sol import solve
    print('rho_M     x (elliptic)   x (series)    sigma (ell)   sigma (ser)   d* (ell)    d* (ser)')
    for rM in (0.1,0.3,0.5,0.7,0.9):
        S=Sol(rM); c=solve((rM**4+rM**-4)/2,m=20000)
        print('%.2f     %.8f    %.8f    %.8f    %.8f   %.8f  %.8f'%(rM,S.x(),c['x'],S.sigma(),c['sigma'],S.dstar(),c['dstar']))
    from scipy.special import gamma
    varpi=gamma(0.25)**2/(2*np.sqrt(2*np.pi))
    print('small rho_M:  x ln2 / rho_M  ->  pi/varpi = %.6f'%(np.pi/varpi))
    for rM in (0.05,0.02,0.01,0.003,0.001):
        S=Sol(rM); x=S.x(); print('   rho_M %.3f  x %.6f  x ln2/rho_M %.6f   d*/x %.6f (4 ln2/pi = %.6f)'%(rM,x,x*np.log(2)/rM,S.dstar()/x,4*np.log(2)/np.pi))

# --- closed relations between the differentials (scripts/corner_islands/DERIVATION.md, section 5) ---
def dsigma_dx(rM):
    """d sigma / d x = (2 ln2/pi) K(k1')/K(k1), k1 = 2/p_M."""
    from scipy.special import ellipk
    k1 = 2 / (rM**2 + rM**-2); return 2 * np.log(2) / np.pi * ellipk(1 - k1 * k1) / ellipk(k1 * k1)
def a0(rM):
    """Mean of F on the cut, frak F at the centre: a0 = (4 sqrt(2 p_M)/pi) [E(k) - K(k)/2], k^2 = (p_M-2)/(2 p_M)."""
    from scipy.special import ellipk, ellipe
    pM = rM**2 + rM**-2; m = (pM - 2) / (2 * pM); return 4 * np.sqrt(2 * pM) / np.pi * (ellipe(m) - ellipk(m) / 2)
