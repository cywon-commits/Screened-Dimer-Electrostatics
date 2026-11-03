import gmpy2, sys, time
from gmpy2 import mpz
from math import log
def Smat(N,I,J):
    # S(i,j)=sum_k C(k,i)C(k,j)2^(N-1-k)
    B={}
    for i in set(I)|set(J): B[i]=[gmpy2.comb(k,i) for k in range(N)]
    return [[sum(B[i][k]*B[j][k]<<(N-1-k) for k in range(max(i,j),N)) for j in J] for i in I]
def lndet(A):
    A=[r[:] for r in A]; m=len(A); prev=mpz(1); sign=1
    for k in range(m-1):
        if A[k][k]==0:
            for r in range(k+1,m):
                if A[r][k]!=0: A[k],A[r]=A[r],A[k]; sign=-sign; break
        akk=A[k][k]
        for r in range(k+1,m):
            ark=A[r][k]; Ar=A[r]; Ak=A[k]
            for c in range(k+1,m): Ar[c]=(Ar[c]*akk-ark*Ak[c])//prev
        prev=akk
    d=sign*A[m-1][m-1]
    return float(gmpy2.log(abs(gmpy2.mpfr(d,precision=64)))) if d!=0 else float('nan')
def lnR(N,i,j): return lndet([[Smat(N,[i],[j])[0][0]]])-N*log(2)
if __name__=='__main__':
    N=int(sys.argv[1]); a0=float(sys.argv[2]); i0=int(a0*N)
    for m in [int(v) for v in sys.argv[3:]]:
        t=time.time(); I=list(range(i0,i0+m)); S=Smat(N,I,I)
        single=sum(lndet([[S[a][a]]]) for a in range(m))
        print(N,a0,m,round((lndet(S)-single)/m,3),'%.0fs'%(time.time()-t),flush=True)
