"""Exact grand-canonical monomer-dimer statistics on the Aztec diamond of order N.

Diagonal frame: cells (x,y) of the Aztec diamond, a=x+y (diagonal index), b=x-y.
Diagonals k=0..2N alternate short (N cells, k even) and long (N+1 cells, k odd); colour = k mod 2.
Short cell i (b=-N+1+2i) touches long cells i, i+1; long cell j touches short cells j-1, j.
State on a diagonal = bits of cells already covered by a dimer from the previous diagonal.
Every cell is processed once as a 'source': covered-from-below / monomer (z) / dimer to next-diagonal
left or right neighbour. Last diagonal: free cells must be monomers.

Outputs per cell: P(mono), P(left bond), P(right bond); from these all bond probabilities.
"""
import numpy as np

def V(a, i): return a.reshape(-1, 2, 1 << i)

def step(A0, A1, i, z, lok, rok, terms=('b', 'm', 'l', 'r')):
    X0, X1 = V(A0, i), V(A1, i)
    N0 = np.zeros_like(X0); N1 = np.zeros_like(X1)
    if 'b' in terms: N0[:, 0] += X0[:, 1]; N0[:, 1] += X1[:, 1]
    if 'm' in terms: N0[:, 0] += z * X0[:, 0]; N0[:, 1] += z * X1[:, 0]
    if 'l' in terms and lok: N0[:, 1] += X0[:, 0]
    if 'r' in terms and rok: N1[:, 0] += X0[:, 0]; N1[:, 1] += X1[:, 0]
    return N0.reshape(-1), N1.reshape(-1)

def stepT(G0, G1, i, z, lok, rok):
    H0, H1 = V(G0, i), V(G1, i)
    Y0 = np.zeros_like(H0); Y1 = np.zeros_like(H1)
    Y0[:, 1] = H0[:, 0]
    Y0[:, 0] = z * H0[:, 0]
    if lok: Y0[:, 0] += H0[:, 1]
    if rok: Y0[:, 0] += H1[:, 0]
    Y1[:, 1] = H0[:, 1]
    Y1[:, 0] = z * H0[:, 1]
    if rok: Y1[:, 0] += H1[:, 1]
    return Y0.reshape(-1), Y1.reshape(-1)

class Aztec:
    def __init__(self, N, z, zmap=None):
        self.N, self.z = N, z
        self.zmap = zmap  # optional list of per-diagonal arrays of site fugacities
        self.K = 2 * N + 1
        self.L = [N if k % 2 == 0 else N + 1 for k in range(self.K)]

    def flags(self, k, i):
        N = self.N
        if k % 2 == 0:   # short -> long : neighbours i, i+1 always exist
            return True, True
        return i >= 1, i <= N - 1   # long -> short : neighbours i-1 (left), i (right)

    def transfer(self, k, v, store=False):
        n = self.L[k]; A0, A1 = v, np.zeros_like(v); F = []
        for i in range(n):
            if store: F.append((A0, A1))
            A0, A1 = step(A0, A1, i, self.zs(k, i), *self.flags(k, i))
        out = np.concatenate([A0, A1]) if k % 2 == 0 else A0.reshape(-1, 2)[:, 0].copy()
        return (out, F) if store else out

    def transferT(self, k, g, store=False):
        n = self.L[k]
        if k % 2 == 0:
            h = g.size // 2; G0, G1 = g[:h].copy(), g[h:].copy()
        else:
            G0 = np.zeros(2 * g.size); G0.reshape(-1, 2)[:, 0] = g; G1 = np.zeros_like(G0)
        B = [None] * n
        for i in reversed(range(n)):
            if store: B[i] = (G0, G1)
            G0, G1 = stepT(G0, G1, i, self.zs(k, i), *self.flags(k, i))
        return (G0, B) if store else G0

    def zs(self, k, i):
        return self.z if self.zmap is None else self.zmap[k][i]

    def final_vec(self):
        n = self.L[-1]; idx = np.arange(1 << n); k = self.K - 1
        w = np.ones(1 << n)
        for i in range(n):
            w *= np.where((idx >> i) & 1, 1.0, self.zs(k, i))
        return w

    def run(self):
        """Returns lnZ and per-diagonal arrays mono[k][i], left[k][i], right[k][i]."""
        K, L = self.K, self.L
        fwd = [np.zeros(1 << L[0])]; fwd[0][0] = 1.0
        scale = [0.0]
        for k in range(K - 1):
            v = self.transfer(k, fwd[-1]); s = np.abs(v).max(); fwd.append(v / s); scale.append(scale[-1] + np.log(s))
        f = self.final_vec()
        Zr = fwd[-1] @ f; lnZ = np.log(Zr) + scale[-1]
        bwd = [None] * K; bwd[-1] = f.copy(); bsc = [0.0] * K
        for k in range(K - 2, -1, -1):
            g = self.transferT(k, bwd[k + 1]); s = np.abs(g).max(); bwd[k] = g / s; bsc[k] = bsc[k + 1] + np.log(s)
        mono, left, right = [], [], []
        for k in range(K):
            n = L[k]; mk = np.zeros(n); lk = np.zeros(n); rk = np.zeros(n)
            if k == K - 1:
                norm = fwd[k] @ bwd[k]
                idx = np.arange(1 << n)
                for i in range(n):
                    freei = ((idx >> i) & 1) == 0
                    mk[i] = (fwd[k] * bwd[k] * freei).sum() / norm
            else:
                _, F = self.transfer(k, fwd[k], store=True)
                _, B = self.transferT(k, bwd[k + 1], store=True)
                # rescale: fwd[k] -> transfer -> . bwd[k+1] has scale factor relative to norm
                ref = self.transfer(k, fwd[k]) @ bwd[k + 1]
                for i in range(n):
                    lo, ro = self.flags(k, i); G0, G1 = B[i]
                    for t, arr in (('m', mk), ('l', lk), ('r', rk)):
                        M0, M1 = step(*F[i], i, self.zs(k, i), lo, ro, terms=(t,))
                        arr[i] = (G0 @ M0 + G1 @ M1) / ref
            mono.append(mk); left.append(lk); right.append(rk)
        return lnZ, mono, left, right

def cell_coords(N, k, i):
    """(x, y) of cell i on diagonal k (Aztec cells 0..2N-1)."""
    a = N - 1 + k
    b = (-N + 1 + 2 * i) if k % 2 == 0 else (-N + 2 * i)
    return (a + b) // 2, (a - b) // 2

def bond_maps(N, left, right):
    """Return dict cell->(sum of incident bond probs, max incident bond prob)."""
    inc = {}
    def add(c, p):
        s, m = inc.get(c, (0.0, 0.0)); inc[c] = (s + p, max(m, p))
    K = 2 * N + 1
    for k in range(K - 1):
        for i in range(len(left[k])):
            c = cell_coords(N, k, i)
            if k % 2 == 0: jl, jr = i, i + 1
            else: jl, jr = i - 1, i
            if left[k][i] > 0:
                add(c, left[k][i]); add(cell_coords(N, k + 1, jl), left[k][i])
            if right[k][i] > 0:
                add(c, right[k][i]); add(cell_coords(N, k + 1, jr), right[k][i])
    return inc
