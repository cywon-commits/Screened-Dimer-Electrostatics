import numpy as np
from scipy.optimize import minimize
from ckp import ent
n=400; th=2*np.pi*(np.arange(n)+.5)/n; T,P=np.meshgrid(th,th)
def R(X,Y):
    z=np.exp(X+1j*T); w=np.exp(Y+1j*P); return np.mean(np.log(np.abs(1+z+w-z*w)))
def entR(a,b):
    r=minimize(lambda q: -(R(q[0],q[1])-a*q[0]-b*q[1]),[0.1,0.1],method='Nelder-Mead',options=dict(xatol=1e-9,fatol=1e-12))
    # ent = min over (X,Y) of R - aX - bY ; R convex => use minimize of that
    r=minimize(lambda q: (R(q[0],q[1])-a*q[0]-b*q[1]),[0.0,0.0],method='Nelder-Mead',options=dict(xatol=1e-10,fatol=1e-13,maxiter=4000))
    return r.fun
print('R(0,0)',R(0,0),2*0.915965594/np.pi)
for hu,hv in ((0.3,0.0),(0.3,-0.4),(0.6,0.2),(-0.5,0.5)):
    a,b=(1+hu)/2,(1+hv)/2
    print(hu,hv,'ronkin',round(entR(a,b),6),'ckp(s=hu-hv,t=hu+hv)',round(float(ent(hu-hv,hu+hv)),6),'ckp(s=hu+hv,t=hu-hv)',round(float(ent(hu+hv,hu-hv)),6),'ckp(2hu,2hv)?',round(float(ent(2*hu,2*hv)) if abs(2*hu)+abs(2*hv)<=2 else -1,6))
