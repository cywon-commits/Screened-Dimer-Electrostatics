"""Faster strip transfer matrix (see tm.py). Changes:
 - sector operator P*T (P applied once; P commutes with T, PT symmetric)
 - bit rotation via reshape/transpose instead of gather
 - density by Hellmann-Feynman: n = z <psi| dT/dz |psi> / (L0 W), dT/dz by central difference of T*psi
usage: python3 tm2.py z W1 W2 ...
"""
import sys, time, numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

def make_T(W, z):
    N = 1 << W
    idx = np.arange(N)
    maskB = (((idx & 1) == 0) & (((idx >> (W - 1)) & 1) == 0)).astype(float)
    def site(a0, a1, x, forced=False):
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
    def T(v):
        a0, a1 = v, zero
        for x in range(W): a0, a1 = site(a0, a1, x)
        b0, b1 = zero, v * maskB
        for x in range(W - 1): b0, b1 = site(b0, b1, x)
        b0, _ = site(b0, b1, W - 1, forced=True)
        return a0 + b0
    return T

def rot(v, W):
    # new[s] = v[rotl(s)] equivalently a cyclic relabelling; any consistent cyclic shift works for projection
    return v.reshape(2, -1).T.reshape(-1)

def P(v, W, sign):
    acc = v.copy(); w = v
    for j in range(1, W):
        w = rot(w, W); acc += (sign ** j) * w
    return acc / W

def sector(W, z, sign, k, vecs=False):
    T = make_T(W, z); N = 1 << W
    op = LinearOperator((N, N), matvec=lambda v: P(T(np.asarray(v).reshape(-1)), W, sign), dtype=float)
    v0 = P(np.random.default_rng(7).standard_normal(N), W, sign)
    if vecs:
        vals, V = eigsh(op, k=k, which='LM', v0=v0, tol=1e-13)
        o = np.argsort(-np.abs(vals)); return vals[o], V[:, o]
    vals = eigsh(op, k=k, which='LM', v0=v0, tol=1e-13, return_eigenvectors=False)
    return vals[np.argsort(-np.abs(vals))], None

if __name__ == '__main__':
    z = float(sys.argv[1])
    for W in [int(w) for w in sys.argv[2:]]:
        t0 = time.time()
        e0, V = sector(W, z, +1, 3, vecs=True)
        psi = V[:, 0]; L0 = e0[0]
        dz = 1e-4 * z
        dT = (make_T(W, z + dz)(psi) - make_T(W, z - dz)(psi)) / (2 * dz)
        n = z * psi @ dT / (L0 * W)
        epi, _ = sector(W, z, -1, 2)
        m1 = np.log(abs(L0) / abs(epi[0])); m2 = np.log(abs(L0) / abs(e0[1]))
        print(f"z {z} W {W:2d} n {n:.10f} m1 {m1:.10f} m2 {m2:.10f} | L0 {L0:.6e} Lpi {epi[0]:+.6e} {epi[1]:+.6e} L0sec {e0[1]:+.6e} {e0[2]:+.6e} | m1^2/n {m1*m1/n:.6f} m2/m1 {m2/m1:.6f} t {time.time()-t0:.0f}s", flush=True)
