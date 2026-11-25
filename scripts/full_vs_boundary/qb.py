"""Exact boundary-model Q_un(N,z) = N/2 - <vacancies on one edge>, from the boundary Pfaffian (src/bpf_exact.py)."""
import sys, time, pickle, os, numpy as np
from flint import arb, ctx
import bpf_exact as B
_G={}
def G(N):
    if N not in _G: _G[N]=B.build(N,B.Rmat(N))
    return _G[N]
def Qb(N,z,prec=None):
    ctx.prec=prec or int(6*N)+400
    za=arb(z); h=arb(10)**(-30)
    d=(B.lndet(N,G(N),za*(1+h))-B.lndet(N,G(N),za*(1-h)))/(2*h)      # total number of vacancies (4 edges)
    return N/2-float(d.mid())/4
if __name__=='__main__':
    for N in (16,32,64):
        t=time.time(); G(N); t1=time.time(); q=Qb(N,10**(-N/16)); print(N,q,'build %.1fs eval %.1fs'%(t1-t,time.time()-t1),flush=True)
