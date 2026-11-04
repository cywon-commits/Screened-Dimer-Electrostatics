"""Cohn-Kenyon-Propp local entropy for domino tilings, slope (s,t) with |s|+|t|<=2 (Thurston heights).
ent(s,t) = (1/pi) sum_i L(pi p_i),  cos(pi A) = [cos(pi s/2) - cos(pi t/2)]/2,
p_a=(A+s/2)/2, p_b=(A-s/2)/2, p_c=(1-A+t/2)/2, p_d=(1-A-t/2)/2.   L = Lobachevsky.
Returns entropy per domino (ent(0,0)=2G/pi); divide by 2 for per square."""
import numpy as np
# Lobachevsky L(x) = -int_0^x ln(2 sin t) dt on [0,pi], tabulated
_x=np.linspace(0,np.pi,200001); _h=_x[1]
_f=-np.log(2*np.sin(np.clip(_x,1e-300,None))); _f[0]=0; _f[-1]=0
# integrate with singular-endpoint correction: near 0, -ln(2t) integrable; use midpoint rule on cells
_mid=(_x[:-1]+_x[1:])/2; _g=-np.log(2*np.sin(_mid))
_L=np.concatenate([[0],np.cumsum(_g*_h)])
# correct the first cell exactly: int_0^h -ln(2t) dt = -h ln(2h) + h
_L[1:]+= (-_h*np.log(2*_h)+_h) - _g[0]*_h
_L[-1:]+=0
def L(x): return np.interp(x,_x,_L)
def ent(s,t,grad=False):
    s=np.asarray(s,float); t=np.asarray(t,float)
    cA=np.clip((np.cos(np.pi*s/2)-np.cos(np.pi*t/2))/2,-1,1); A=np.arccos(cA)/np.pi
    p=np.stack([(A+s/2)/2,(A-s/2)/2,(1-A+t/2)/2,(1-A-t/2)/2]); p=np.clip(p,0,1)
    e=L(np.pi*p).sum(0)/np.pi
    if not grad: return e
    sA=np.maximum(np.sin(np.pi*A),1e-12)
    As=np.sin(np.pi*s/2)/(4*sA); At=-np.sin(np.pi*t/2)/(4*sA)
    lg=-np.log(2*np.maximum(np.sin(np.pi*p),1e-15))      # dL/dx at pi p, times pi /pi
    ps=np.stack([(As+.5)/2,(As-.5)/2,-As/2,-As/2]); pt=np.stack([At/2,At/2,(-At+.5)/2,(-At-.5)/2])
    return e,(lg*ps).sum(0),(lg*pt).sum(0)
if __name__=='__main__':
    G=0.915965594177219
    print(ent(0,0),2*G/np.pi, ent(2,0),ent(1,1),ent(0.999,0.999))
    e,gs,gt=ent(0.3,-0.4,True); h=1e-6
    print(gs,(ent(0.3+h,-0.4)-ent(0.3-h,-0.4))/2/h, gt,(ent(0.3,-0.4+h)-ent(0.3,-0.4-h))/2/h)
