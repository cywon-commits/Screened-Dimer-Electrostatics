"""Crossover (crossing <-> Delannoy) single-pair exponent for corner pairs.
R(i,j) = sum_{k<N} C(k,i)C(k,j)/2^(k+1),  i=pN, j=qN.
(1/N) ln R -> E(p,q) = max_{kappa in [max(p,q),1]} kappa[H(p/kappa)+H(q/kappa)] - kappa ln2
kappa=1 (top level) <=> crossing regime D=1/2-(1-p)(1-q) >= 0 ; interior kappa <=> Delannoy regime.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from math import lgamma, log
def H(a): a=np.clip(a,1e-15,1-1e-15); return -(a*np.log(a)+(1-a)*np.log(1-a))
def E(p,q):
    f=lambda k: -(k*(H(p/k)+H(q/k))-k*np.log(2))
    lo=max(p,q)+1e-12
    r=minimize_scalar(f,bounds=(lo,1.0),method='bounded',options=dict(xatol=1e-12))
    return max(-r.fun,-f(1.0)), r.x
def lnR_exact(N,i,j):
    ks=np.arange(max(i,j),N)
    lt=np.array([lgamma(k+1)-lgamma(i+1)-lgamma(k-i+1)+lgamma(k+1)-lgamma(j+1)-lgamma(k-j+1)-(k+1)*log(2) for k in ks])
    m=lt.max(); return m+np.log(np.exp(lt-m).sum())
if __name__=='__main__':
    for a in (0.1,0.2,0.25,0.29,0.3,0.4,0.5):
        e,k=E(a,a); ec=2*H(a)-np.log(2)
        ex=[lnR_exact(N,int(a*N),int(a*N))/N for N in (400,1600)]
        print(f'a={a}: E={e:.4f} kappa*={k:.3f} crossing-formula={ec:.4f}  exact/N(400,1600)={ex[0]:.4f},{ex[1]:.4f}')
