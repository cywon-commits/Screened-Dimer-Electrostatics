import asym as A, sys
from math import log
def S_val(N,a,b,c,d,cache):
    f=lambda x: int(round(x*N))
    # positions as fractions along edges in P.bd order (pos index / N)
    plus=[(0,f(a)),(2,f(c))]; minus=[(1,f(b)),(3,f(d))]
    ld=A.lndet_exact(N,plus,minus,cache)
    lab=A.lndet_exact(N,[plus[0]],[minus[0]],cache); lcd=A.lndet_exact(N,[plus[1]],[minus[1]],cache)
    FB=ld-lab-lcd
    pa,qd=A.corner_frac(N,0,f(a),3,f(d)); pc,qb=A.corner_frac(N,2,f(c),1,f(b))
    return FB, FB+log(abs(A.Dfun(pa,qd)))+log(abs(A.Dfun(pc,qb)))+log(2)
