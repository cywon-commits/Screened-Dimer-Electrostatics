"""Row-frame transfer matrix for monomer-dimer on a strip of width W.
periodic=True: cylinder; False: open strip (straight domino edges at y=0 and y=W-1 in the transverse direction).
Per-site fugacity zs[x]. T is symmetric.
Functions:
  gap_Q(W, periodic, Q): z=0, top eigenvalue in sector |Q| (Q = sum_{x in s} (-1)^x)  -> m_Q = ln(L0/LQ)
  profile(W, z, periodic): n(x) = z_x d ln L0 / d z_x  via  psi.(T psi - T_{z_x=0} psi)/L0
"""
import sys, numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

def make_T(W, zs, periodic):
    N = 1 << W
    zs = np.asarray(zs, float)
    idx = np.arange(N)
    maskB = (((idx & 1) == 0) & (((idx >> (W - 1)) & 1) == 0)).astype(float)
    def site(a0, a1, x, z, forced=False):
        A0 = a0.reshape(-1, 2, 1 << x); A1 = a1.reshape(-1, 2, 1 << x)
        n0 = np.empty_like(A0)
        if forced:
            n0[:, 0, :] = A0[:, 0, :]; n0[:, 1, :] = 0
            return n0.reshape(-1), None
        n1 = np.zeros_like(A1)
        np.multiply(A0[:, 0, :], z, out=n0[:, 0, :]); n0[:, 0, :] += A0[:, 1, :]; n0[:, 0, :] += A1[:, 0, :]
        n0[:, 1, :] = A0[:, 0, :]
        n1[:, 0, :] = A0[:, 0, :]
        return n0.reshape(-1), n1.reshape(-1)
    zero = np.zeros(N)
    def T(v, zz=None):
        zz = zs if zz is None else zz
        a0, a1 = v, zero
        for x in range(W): a0, a1 = site(a0, a1, x, zz[x])
        # open strip: a horizontal dimer may not be left pending after the last site
        out = a0
        if periodic:
            b0, b1 = zero, v * maskB
            for x in range(W - 1): b0, b1 = site(b0, b1, x, zz[x])
            b0, _ = site(b0, b1, W - 1, zz[W - 1], forced=True)
            out = out + b0
        return out
    return T

def Qof(W):
    idx = np.arange(1 << W); Q = np.zeros(1 << W, int)
    for x in range(W): Q += ((idx >> x) & 1) * (1 if x % 2 == 0 else -1)
    return Q

def top(W, T, mask=None, k=1):
    N = 1 << W
    if mask is None:
        op = LinearOperator((N, N), matvec=lambda v: T(np.asarray(v).reshape(-1)), dtype=float)
        v0 = np.ones(N)
    else:
        m = mask.astype(float)
        op = LinearOperator((N, N), matvec=lambda v: m * T(m * np.asarray(v).reshape(-1)), dtype=float)
        v0 = m * np.random.default_rng(3).random(N)
    vals, V = eigsh(op, k=k, which='LM', v0=v0, tol=1e-13)
    o = np.argsort(-np.abs(vals)); return vals[o], V[:, o]

def gap_Q(W, periodic, Qs=(1, 2)):
    T = make_T(W, np.zeros(W), periodic); Q = Qof(W)
    L0 = abs(top(W, T, np.abs(Q) == 0)[0][0])
    return {q: np.log(L0 / abs(top(W, T, np.abs(Q) == q)[0][0])) for q in Qs}

def profile(W, z, periodic=False):
    zs = np.full(W, z); T = make_T(W, zs, periodic)
    vals, V = top(W, T); L0 = vals[0]; psi = V[:, 0]
    Tpsi = T(psi)
    n = np.empty(W)
    for x in range(W):
        z0 = zs.copy(); z0[x] = 0
        n[x] = psi @ (Tpsi - T(psi, z0)) / L0
    return L0, n

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'gap':
        for W in [int(w) for w in sys.argv[2:]]:
            gc = gap_Q(W, True); go = gap_Q(W, False)
            print(f"W {W:2d}  cyl: m1*W/2pi {gc[1]*W/(2*np.pi):.6f} m2*W/2pi {gc[2]*W/(2*np.pi):.6f}   "
                  f"open: m1*W'/pi {go[1]*(W+1)/np.pi:.6f} m2*W'/pi {go[2]*(W+1)/np.pi:.6f}", flush=True)
    else:
        z = float(sys.argv[2])
        for W in [int(w) for w in sys.argv[3:]]:
            L0, n = profile(W, z)
            print(f"z {z} W {W} L0 {L0:.8e} n " + ' '.join(f'{v:.10f}' for v in n), flush=True)
