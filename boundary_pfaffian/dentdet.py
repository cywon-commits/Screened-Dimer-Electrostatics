import twopair as T
from math import log
def corner_idx(N, ep, pp, em, pm):
    """indices (i,j) from the shared corner for + dent on edge ep at pos pp and - dent on edge em at pos pm."""
    pair = {(0,1):('bottom'),(2,1):('right'),(2,3):('top'),(0,3):('left')}[(ep,em)]
    if pair=='bottom': return N-1-pp, pm
    if pair=='right':  return pp, N-1-pm
    if pair=='top':    return N-1-pp, pm
    if pair=='left':   return pp, N-1-pm
def lndet(N, plus, minus):
    """plus/minus: lists of (edge,pos). Returns ln det[R] (sign-insensitive)."""
    from fractions import Fraction
    M=[[T.S(N,*corner_idx(N,ep,pp,em,pm)) for (em,pm) in minus] for (ep,pp) in plus]
    m=len(M); A=[r[:] for r in M]; sign=1; prev=1
    for k in range(m-1):
        if A[k][k]==0:
            for r in range(k+1,m):
                if A[r][k]!=0: A[k],A[r]=A[r],A[k]; sign=-sign; break
        for r in range(k+1,m):
            for c in range(k+1,m): A[r][c]=(A[r][c]*A[k][k]-A[r][k]*A[k][c])//prev
        prev=A[k][k]
    det=sign*A[m-1][m-1]
    return T.lnbig(abs(det))-m*N*log(2), det
