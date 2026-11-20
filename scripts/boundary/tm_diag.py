"""Diagonal-frame transfer matrix: strip bounded by two monochrome staircase edges.
Diagonals t = x+y; each diagonal has W sites, all of colour (t mod 2).
Even diagonal site i at transverse position p=2i, odd diagonal site j at p=2j+1, p in [0, 2W-1].
Left edge p=0 (even colour), right edge p=2W-1 (odd colour).
E: even diag -> odd diag;  O: odd diag -> even diag;  T2 = O E (not symmetric).
State bits = sites of the diagonal already matched from the previous diagonal.
Densities by forward/backward vectors with left & right Perron vectors.
"""
import sys, numpy as np
from scipy.sparse.linalg import LinearOperator, eigs

def V(a, i):  # view with bit i as axis 1
    return a.reshape(-1, 2, 1 << i)

def E_step(A0, A1, i, z, W):
    X0, X1 = V(A0, i), V(A1, i)
    N0 = np.zeros_like(X0); N1 = np.zeros_like(X1)
    N0[:, 0] = X0[:, 1] + z * X0[:, 0]
    N1[:, 0] = X0[:, 0]
    if i >= 1:
        N0[:, 1] = X1[:, 1] + z * X1[:, 0] + X0[:, 0]
        N1[:, 1] = X1[:, 0]
    return N0.reshape(-1), N1.reshape(-1)

def E_stepT(G0, G1, i, z, W):
    H0, H1 = V(G0, i), V(G1, i)
    Y0 = np.zeros_like(H0); Y1 = np.zeros_like(H1)
    Y0[:, 1] = H0[:, 0]
    Y0[:, 0] = z * H0[:, 0] + H1[:, 0]
    if i >= 1:
        Y0[:, 0] += H0[:, 1]
        Y1[:, 1] = H0[:, 1]
        Y1[:, 0] = z * H0[:, 1] + H1[:, 1]
    return Y0.reshape(-1), Y1.reshape(-1)

def O_step(A0, A1, j, z, W):
    X0, X1 = V(A0, j), V(A1, j)
    N0 = np.zeros_like(X0); N1 = np.zeros_like(X1)
    N0[:, 0] = X0[:, 1] + z * X0[:, 0]
    N0[:, 1] = X0[:, 0] + X1[:, 1] + z * X1[:, 0]
    if j < W - 1:
        N1[:, 0] = X0[:, 0]
        N1[:, 1] = X1[:, 0]
    return N0.reshape(-1), N1.reshape(-1)

def O_stepT(G0, G1, j, z, W):
    H0, H1 = V(G0, j), V(G1, j)
    Y0 = np.zeros_like(H0); Y1 = np.zeros_like(H1)
    Y0[:, 1] = H0[:, 0]
    Y0[:, 0] = z * H0[:, 0] + H0[:, 1]
    Y1[:, 1] = H0[:, 1]
    Y1[:, 0] = z * H0[:, 1]
    if j < W - 1:
        Y0[:, 0] += H1[:, 0]
        Y1[:, 0] += H1[:, 1]
    return Y0.reshape(-1), Y1.reshape(-1)

def mono(A0, A1, i, z):  # only the monomer term at site i (same form for E and O)
    X0, X1 = V(A0, i), V(A1, i)
    N0 = np.zeros_like(X0)
    N0[:, 0] = z * X0[:, 0]; N0[:, 1] = z * X1[:, 0]
    return N0.reshape(-1), np.zeros_like(A1)

def E_fin(N0, N1):
    return np.concatenate([N0.reshape(-1, 2)[:, 0], N1.reshape(-1, 2)[:, 0]])

def E_finT(g):
    h = g.size // 2
    G0 = np.zeros(g.size); G1 = np.zeros(g.size)
    G0.reshape(-1, 2)[:, 0] = g[:h]; G1.reshape(-1, 2)[:, 0] = g[h:]
    return G0, G1

def apply_E(v, z, W, store=False):
    A0, A1 = v, np.zeros_like(v); F = []
    for i in range(W):
        if store: F.append((A0, A1))
        A0, A1 = E_step(A0, A1, i, z, W)
    return (E_fin(A0, A1), F) if store else E_fin(A0, A1)

def apply_ET(g, z, W, store=False):
    G0, G1 = E_finT(g); B = [None] * W
    for i in reversed(range(W)):
        if store: B[i] = (G0, G1)
        G0, G1 = E_stepT(G0, G1, i, z, W)
    return (G0, B) if store else G0

def apply_O(v, z, W, store=False):
    A0, A1 = v, np.zeros_like(v); F = []
    for j in range(W):
        if store: F.append((A0, A1))
        A0, A1 = O_step(A0, A1, j, z, W)
    return (A0, F) if store else A0

def apply_OT(g, z, W, store=False):
    G0, G1 = g, np.zeros_like(g); B = [None] * W
    for j in reversed(range(W)):
        if store: B[j] = (G0, G1)
        G0, G1 = O_stepT(G0, G1, j, z, W)
    return (G0, B) if store else G0

def perron(W, z):
    N = 1 << W
    R = LinearOperator((N, N), matvec=lambda v: apply_O(apply_E(np.real(v).reshape(-1), z, W), z, W), dtype=float)
    L = LinearOperator((N, N), matvec=lambda v: apply_ET(apply_OT(np.real(v).reshape(-1), z, W), z, W), dtype=float)
    lr, vr = eigs(R, k=1, which='LR', v0=np.ones(N), tol=1e-13)
    ll, vl = eigs(L, k=1, which='LR', v0=np.ones(N), tol=1e-13)
    psi = np.real(vr[:, 0]); phi = np.real(vl[:, 0])
    psi *= np.sign(psi.sum()); phi *= np.sign(phi.sum())
    return np.real(lr[0]), np.real(ll[0]), psi, phi

def densities(W, z):
    Lam, Lam2, psi, phi = perron(W, z)
    norm = Lam * (phi @ psi)
    # even sites: phi^T O [E with mono_i] psi
    g = apply_OT(phi, z, W)
    Epsi, F = apply_E(psi, z, W, store=True)
    _, B = apply_ET(g, z, W, store=True)
    ne = np.empty(W)
    for i in range(W):
        M0, M1 = mono(*F[i], i, z)
        G0, G1 = B[i]
        ne[i] = (G0 @ M0 + G1 @ M1) / norm
    _, F2 = apply_O(Epsi, z, W, store=True)
    _, B2 = apply_OT(phi, z, W, store=True)
    no = np.empty(W)
    for j in range(W):
        M0, M1 = mono(*F2[j], j, z)
        G0, G1 = B2[j]
        no[j] = (G0 @ M0 + G1 @ M1) / norm
    return Lam, Lam2, ne, no

if __name__ == '__main__':
    z = float(sys.argv[1])
    for W in [int(w) for w in sys.argv[2:]]:
        Lam, Lam2, ne, no = densities(W, z)
        print(f"z {z} W {W} Lam {Lam:.10e} LamL {Lam2:.10e}", flush=True)
        print("ne " + ' '.join(f'{v:.10f}' for v in ne), flush=True)
        print("no " + ' '.join(f'{v:.10f}' for v in no), flush=True)
