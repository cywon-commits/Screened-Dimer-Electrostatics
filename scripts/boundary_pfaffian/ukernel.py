"""Uniform (crossing <-> Delannoy) two-point function of boundary dents in the Aztec diamond.
R(i,j) = sum_{k<N} T(k),  T(k) = C(k,i)C(k,j)/2^(k+1).
Near the top level k = N-1 - l:  ln T ~ ln T(N-1) - alpha*l - l^2/(2 N s),
   alpha = -ln[2(1-p)(1-q)] = -ln(1-2D),   1/s = p/(1-p) + q/(1-q),   p=i/N, q=j/N.
Uniform sum:  R ~ T(N-1) * G,  G = sqrt(pi N s/2) exp(alpha^2 N s/2) erfc(alpha sqrt(N s/2)) + 1/2.
Crossing (alpha*sqrt(N)>>1): G -> 1/alpha ~ 1/(2D)  (Cauchy).  D -> 0: G ~ sqrt(pi N s/2) (finite).
Delannoy (alpha<0): G ~ sqrt(2 pi N s) exp(alpha^2 N s/2)  (interior saddle, Gaussian approx).
Effective D:  D_eff = 1/(2G)  replaces D in the Cauchy kernel -ln|D|.
"""
import numpy as np
from math import lgamma, log
from scipy.special import erfcx
def lnT_top(N,i,j):
    k=N-1
    if i>k or j>k: return -np.inf
    return (lgamma(k+1)-lgamma(i+1)-lgamma(k-i+1))+(lgamma(k+1)-lgamma(j+1)-lgamma(k-j+1))-(k+1)*log(2)
def lnG(N,p,q):
    arg=2*(1-p)*(1-q); alpha=-np.log(arg)
    s=1.0/(p/(1-p)+q/(1-q)) if p>0 or q>0 else 1e9
    y=alpha*np.sqrt(N*s/2)
    # sqrt(pi N s/2) * exp(y^2) erfc(y) = sqrt(pi N s/2)*erfcx(y)  (stable for y>0); for y<0 erfcx grows like 2exp(y^2)
    return np.log(np.sqrt(np.pi*N*s/2)*erfcx(y)+0.5)
def lnR_uniform(N,i,j):
    p,q=(i+0.5)/N,(j+0.5)/N
    return lnT_top(N,i,j)+lnG(N,p,q)
def lnD_eff(N,i,j):
    p,q=(i+0.5)/N,(j+0.5)/N
    return -np.log(2)-lnG(N,p,q)
