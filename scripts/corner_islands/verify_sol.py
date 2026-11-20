"""A-posteriori verification of the hodograph solution in the interior of the fundamental triangle S-M-O.
For zeta in the sector {|zeta|<1, -pi/4<arg<0}: build (u,v)(zeta), slopes h(zeta) from the entropy function itself
(solve grad ent(h) = -(ln|z|, ln|w|)/2), and test
  (1) the map zeta -> (u,v) is one-to-one onto the liquid part of the triangle (Jacobian of one sign, image inside),
  (2) h is a gradient: d_v h_u = d_u h_v,
  (3) Euler-Lagrange: d_u X + d_v Y = 0 with (X,Y) = (ln|z|, ln|w|) = -2 grad ent(h)."""
import numpy as np, sys
from scipy.optimize import brentq, fsolve
from exact_sol import solve, cheb
from ckp import ent
def uv(zeta, JM, cc):
    c0,c1,a=cc; n=np.arange(len(a))
    t=-zeta**4; J=-(t+1/t)/2
    B=(c0-J)/c1; sq=np.sqrt(B*B-1+0j); e=B-sq; e=np.where(np.abs(e)<1,e,B+sq); sqB=B-e
    Phi2=-(2**0.75/(2*c1))*(a*e[...,None]**n).sum(-1)/sqB
    X2=(1-J)**0.25*(JM-J)**0.5
    mu=(1-zeta**2)/(1+zeta**2)
    Psi=np.exp(-1j*np.pi/4)*np.sqrt(mu)*X2*Phi2+1-1j*mu
    lam=1j*mu; u=-Psi.imag/lam.imag; v=Psi.real+lam.real*u
    z=1j*(1+zeta)/(1-zeta); w=(1+z)/(z-1)
    return u,v,z,w
def slopes(z,w,sx,sy):
    X=np.log(abs(z)); Y=np.log(abs(w))
    def f(h):
        hu,hv=np.clip(h,-0.999999,0.999999); e,es,et=ent(hu-hv,hu+hv,True)
        return [es+et-sx*X/2, -es+et-sy*Y/2]
    return fsolve(f,[0.1,0.1],xtol=1e-13)
x=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
JM=1+np.exp(brentq(lambda l: solve(1+np.exp(l))['x']-x,-12,11,xtol=1e-12)); cc=cheb(JM)
# sign convention fixed on the edge SM: h_v = 0 and |h_u| = (4/pi) arctan(rho)
rho=0.9; u,v,z,w=uv(np.array([rho-1e-9j]),JM,cc)
for sx in (1,-1):
    for sy in (1,-1):
        h=slopes(z[0],w[0],sx,sy)
        print('signs',sx,sy,'h on SM =',np.round(h,5),' expected |h_u| =',round(4/np.pi*np.arctan(rho),5),' u,v =',round(float(u[0]),5),round(float(v[0]),6))

def fields(zeta,sy):
    u,v,z,w=uv(np.array([zeta]),JM,cc); h=slopes(z[0],w[0],1,sy)
    return np.array([u[0],v[0],h[0],h[1],np.log(abs(z[0])),np.log(abs(w[0]))])
rs=np.linspace(0.08,0.97,12); ps=-np.pi/4*np.linspace(0.06,0.94,9); d=1e-5
for sy in (1,-1):
    worst_c=worst_e=0; jac=[]; inside=True; pts=[]
    for r in rs:
        for p in ps:
            zt=r*np.exp(1j*p)
            f0=fields(zt,sy); fr=(fields(zt+d,sy)-fields(zt-d,sy))/(2*d); fi=(fields(zt+1j*d,sy)-fields(zt-1j*d,sy))/(2*d)
            Jm=np.array([[fr[0],fi[0]],[fr[1],fi[1]]])          # d(u,v)/d(Re zeta, Im zeta)
            G=np.linalg.solve(Jm.T,np.array([fr[2:],fi[2:]]))    # rows: d/du, d/dv of (hu,hv,X,Y)
            curl=G[1,0]-G[0,1]; el=G[0,2]+G[1,3]; scale=np.abs(G).max()
            worst_c=max(worst_c,abs(curl)/scale); worst_e=max(worst_e,abs(el)/scale); jac.append(np.linalg.det(Jm))
            u,v=f0[:2]; inside&=(v>-1e-9 and v<=u+1e-9 and u<1+1e-9); pts.append((u,v))
    jac=np.array(jac)
    print('x=%.2f sign_y=%+d: max |curl h|/scale = %.1e   max |div(X,Y)|/scale = %.1e   Jacobian sign constant: %s (min %.2e, max %.2e)   image inside triangle S-M-O: %s'%(x,sy,worst_c,worst_e,np.all(jac>0) or np.all(jac<0),jac.min(),jac.max(),inside))
