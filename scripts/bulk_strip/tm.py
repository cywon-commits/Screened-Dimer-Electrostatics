"""Row-to-row transfer matrix of the monomer-dimer model on a periodic strip of width W.
State s (W bits): bit x = site x of the current row already covered by a vertical dimer from below.
T is symmetric; sectors by transverse momentum p (0 or pi) via projection.
m1 = ln(L0/|L_pi|)  (charge, C-odd),  m2 = ln(L0/|L_0^(1)|)  (C-even),  n = z dlnL0/dz / W.
"""
import sys, numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

def make_T(W, z):
    N = 1 << W
    def site(a0, a1, x, last_forced=False):
        A0 = a0.reshape(-1, 2, 1 << x); A1 = a1.reshape(-1, 2, 1 << x)
        n0 = np.empty_like(A0); n1 = np.zeros_like(A1)
        if last_forced:   # site covered by wrap dimer: only free & h=0 allowed, new bit 0
            n0[:, 0, :] = A0[:, 0, :]; n0[:, 1, :] = 0
            return n0.reshape(-1), n1.reshape(-1)
        n0[:, 0, :] = A0[:, 1, :] + A1[:, 0, :] + z * A0[:, 0, :]
        n0[:, 1, :] = A0[:, 0, :]
        n1[:, 0, :] = A0[:, 0, :]
        return n0.reshape(-1), n1.reshape(-1)
    idx = np.arange(N)
    maskB = ((idx & 1) == 0) & ((idx >> (W - 1)) & 1 == 0)
    def T(v):
        v = np.asarray(v).reshape(-1)
        a0, a1 = v.copy(), np.zeros(N)              # pass A: no wrap dimer
        for x in range(W): a0, a1 = site(a0, a1, x)
        out = a0
        b0, b1 = np.zeros(N), np.where(maskB, v, 0.0)  # pass B: wrap dimer (W-1,0)
        for x in range(W - 1): b0, b1 = site(b0, b1, x)
        b0, b1 = site(b0, b1, W - 1, last_forced=True)
        return out + b0
    return T

def rot_perm(W):
    s = np.arange(1 << W)
    return ((s << 1) | (s >> (W - 1))) & ((1 << W) - 1)

def sector_op(W, z, sign):
    T = make_T(W, z); perm = rot_perm(W); N = 1 << W
    def P(v):
        acc = v.copy(); w = v
        for j in range(1, W):
            w = w[perm]; acc += (sign ** j) * w
        return acc / W
    return LinearOperator((N, N), matvec=lambda v: P(T(P(np.asarray(v).reshape(-1)))), dtype=float)

def levels(W, z, sign, k):
    op = sector_op(W, z, sign)
    rng = np.random.default_rng(1); v0 = rng.standard_normal(1 << W)
    vals = eigsh(op, k=k, which='LM', v0=v0, tol=1e-12, return_eigenvectors=False)
    return sorted(vals, key=lambda x: -abs(x))

if __name__ == '__main__':
    z = float(sys.argv[1]); Ws = [int(w) for w in sys.argv[2:]]
    for W in Ws:
        e0 = levels(W, z, +1, 3)
        L0 = e0[0]
        dz = 1e-4 * z
        Lp = levels(W, z + dz, +1, 1)[0]; Lm = levels(W, z - dz, +1, 1)[0]
        n = z * (np.log(abs(Lp)) - np.log(abs(Lm))) / (2 * dz) / W
        epi = levels(W, z, -1, 2)
        m1 = np.log(abs(L0) / abs(epi[0])); m1b = np.log(abs(L0) / abs(epi[1]))
        m2 = np.log(abs(L0) / abs(e0[1])); m2b = np.log(abs(L0) / abs(e0[2]))
        print(f"z {z} W {W:2d} n {n:.10f} m1 {m1:.10f} m2 {m2:.10f} | L0 {L0:.6e} Lpi {epi[0]:+.6e} {epi[1]:+.6e} L0sec {e0[1]:+.6e} {e0[2]:+.6e} | m1^2/n {m1*m1/n:.6f} m2/m1 {m2/m1:.6f} next_pi {m1b:.5f} next_0 {m2b:.5f}", flush=True)
