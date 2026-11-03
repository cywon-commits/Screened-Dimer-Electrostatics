"""Exact partition function and boundary-monomer densities for the Aztec diamond of order N with
monomer fugacities w_a allowed only on boundary cells (all other cells must be covered by dimers).

Z(w) = Pf(K + B),  B_ab = U_a U_b w_a w_b sgn(b - a) on boundary cells in counter-clockwise order.
U_a: signs fixed from the inverse Kasteleyn matrix (rule verified against exact transfer for N <= 7).
Large N: det(K + B) = det K * det(I + B G),  G = K^{-1} restricted to boundary cells (sparse solves).
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla

def region(N):
    cells = [(x, y) for y in range(2 * N) for x in range(2 * N) if abs(x - N + .5) + abs(y - N + .5) <= N]
    return cells, {c: i for i, c in enumerate(cells)}

def kasteleyn_sparse(N):
    cells, idx = region(N); r, c, v = [], [], []
    for (x, y), i in idx.items():
        j = idx.get((x + 1, y))
        if j is not None: r += [i, j]; c += [j, i]; v += [1, -1]
        j = idx.get((x, y + 1))
        if j is not None: s = (-1) ** x; r += [i, j]; c += [j, i]; v += [s, -s]
    n = len(cells)
    return cells, idx, sp.csc_matrix((v, (r, c)), shape=(n, n))

def boundary(N, cells, idx):
    bd = [c for c in cells if sum((c[0] + dx, c[1] + dy) in idx for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1))) < 4]
    ang = np.array([np.arctan2(c[1] - N + .5, c[0] - N + .5) for c in bd])
    return [bd[k] for k in np.argsort(ang)]

class BoundaryPf:
    def __init__(self, N):
        self.N = N
        cells, idx, K = kasteleyn_sparse(N)
        self.cells, self.idx = cells, idx
        bd = boundary(N, cells, idx); self.bd = bd; m = len(bd); assert m == 4 * N
        ii = np.array([idx[c] for c in bd]); self.ii = ii
        lu = spla.splu(K)
        E = np.zeros((len(cells), m)); E[ii, np.arange(m)] = 1.0
        X = lu.solve(E)                        # columns = K^{-1} e_b
        self.G = X[ii, :]                      # G[a,b] = (K^{-1})_{ii[a], ii[b]}
        d = lu.U.diagonal(); self.logdetK = np.sum(np.log(np.abs(d)))  # |det K|
        G = self.G; e = np.arange(m) // N
        # U_a U_b = sign(G_ab) (-1)^(e_a+e_b) sign(b-a) for opposite-colour pairs (verified N<=12);
        # propagate from U_0 = 1 always through the largest |G_ab| to stay clear of round-off.
        U = np.zeros(m); U[0] = 1.0; known = np.zeros(m, bool); known[0] = True
        absG = np.abs(G)
        while not known.all():
            sub = absG[np.ix_(known, ~known)]
            k = np.unravel_index(np.argmax(sub), sub.shape)
            a = np.where(known)[0][k[0]]; b = np.where(~known)[0][k[1]]
            U[b] = U[a] * np.sign(G[a, b]) * (-1) ** (e[a] + e[b]) * np.sign(b - a)
            known[b] = True
        self.U = U
        a, b = np.meshgrid(np.arange(m), np.arange(m), indexing='ij')
        self.S = np.sign(b - a).astype(float)

    def logZ_and_n(self, w):
        """w: fugacities of the boundary cells (in self.bd order). Returns lnZ and monomer density n_a."""
        v = self.U * w
        B = self.S * np.outer(v, v)
        M = np.eye(len(v)) + B @ self.G
        sgn, ld = np.linalg.slogdet(M)
        lnZ = 0.5 * (self.logdetK + ld)
        Y = np.linalg.solve(M.T, self.G.T).T          # d logdet M = tr(Y dB)
        # d/dv_a of tr(Y B) with B_ij = S_ij v_i v_j
        grad = (Y.T * self.S).sum(axis=1) * 0  # placeholder replaced below
        SY = self.S * Y.T                      # (S_ij Y_ji)
        grad = SY @ v + (self.S.T * Y).T @ v   # sum_j S_aj Y_ja v_j + sum_i S_ia Y_ai v_i
        n = 0.5 * v * grad                     # w d lnZ/dw = v d/dv (since v = U w)
        return lnZ, n
