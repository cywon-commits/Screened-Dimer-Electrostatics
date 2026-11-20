import numpy as np, twopair as T, dentdet as D
from math import log
SGN={(0,1):1,(2,1):-1,(2,3):1,(0,3):1}
def corner_frac(N,ep,pp,em,pm):
    i,j=D.corner_idx(N,ep,pp,em,pm); return i/N, j/N
def Dfun(p,q): return 0.5-(1-p)*(1-q)
def lndet_exact(N,plus,minus,cache):
    M=[]
    for (ep,pp) in plus:
        row=[]
        for (em,pm) in minus:
            key=(ep,pp,em,pm)
            if key not in cache: cache[key]=T.S(N,*D.corner_idx(N,ep,pp,em,pm))
            row.append(SGN[(ep,em)]*cache[key])
        M.append(row)
    m=len(M); A=[r[:] for r in M]; sign=1; prev=1
    for k in range(m-1):
        if A[k][k]==0:
            for r in range(k+1,m):
                if A[r][k]!=0: A[k],A[r]=A[r],A[k]; sign=-sign; break
        for r in range(k+1,m):
            for c in range(k+1,m): A[r][c]=(A[r][c]*A[k][k]-A[r][k]*A[k][c])//prev
        prev=A[k][k]
    det=sign*A[m-1][m-1]
    return T.lnbig(abs(det))-m*N*log(2)
