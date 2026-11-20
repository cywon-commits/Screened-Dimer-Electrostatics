"""Boundary-monomer Aztec diamond in arbitrary precision.
|K^-1| between boundary cells of adjacent edges = R(i,j) = sum_{k<N} C(k,i)C(k,j)/2^(k+1), i,j counted from the shared corner
(Ciucu-Krattenthaler Thm 2.1).  Z^2/Z0^2 = det(I + W S W Gt),  Gt_bc = R_bc (-1)^(e_b+e_c) sgn(c-b),  S_ab = sgn(b-a)."""
import sys, numpy as np
from flint import arb, arb_mat, fmpz, fmpq, ctx
def Rmat(N):
    binom=[[fmpz(1) if i==0 else fmpz(0) for i in range(N)] for k in range(N)]
    for k in range(1,N):
        for i in range(1,k+1): binom[k][i]=binom[k-1][i-1]+(binom[k-1][i] if i<=k-1 else 0)
    R=[[fmpq(0)]*N for _ in range(N)]
    for i in range(N):
        for j in range(i,N):
            s=fmpz(0)
            for k in range(j,N): s+=binom[k][i]*binom[k][j]*(fmpz(1)<<(N-1-k))
            R[i][j]=R[j][i]=fmpq(s,fmpz(1)<<N)
    return R
def build(N,R):
    m=4*N; Gt=[[fmpq(0)]*m for _ in range(m)]
    for a in range(m):
        ea,pa=divmod(a,N)
        for eb in ((ea+1)%4,(ea-1)%4):
            for pb in range(N):
                b=eb*N+pb
                r=R[N-1-pa][pb] if eb==(ea+1)%4 else R[pa][N-1-pb]
                Gt[a][b]=r*((-1)**(ea+eb))*(1 if b>a else -1)
    return Gt
def lndet(N,Gt,z):
    m=4*N; M=arb_mat(m,m)
    # (W S W Gt)_ac = z^2 sum_b sgn(b-a) Gt_bc  -> prefix sums over b
    z2=z*z
    for c in range(m):
        col=[Gt[b][c] for b in range(m)]; tot=sum(col,fmpq(0)); pre=fmpq(0)
        for a in range(m):
            # sum_{b>a} - sum_{b<a}
            val=(tot-pre-col[a])-pre; pre+=col[a]
            M[a,c]=z2*arb(val.p)/arb(val.q)+(1 if a==c else 0)
    return M.det().log()/2
def sigma(N,x,Gt,prec):
    ctx.prec=prec
    z=arb(2)**(arb(-x)*N/2); h=arb(10)**(-30)
    d=(lndet(N,Gt,z*(1+h))-lndet(N,Gt,z*(1-h)))/(2*h)     # = sum of n
    return 1-float(d.mid())/(2*N), float(d.rad())
if __name__=='__main__':
    N=int(sys.argv[1]); xs=[float(v) for v in sys.argv[2].split(',')]
    R=Rmat(N); Gt=build(N,R)
    for x in xs:
        s,rad=sigma(N,x,Gt,int(6*N)+400); print('N %d x %.2f sigma %.8f (rad %.1e)'%(N,x,s,rad),flush=True)
