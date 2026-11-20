"""Lozenge 'capacitor': frozen parallelogram strip (N layers, ring of M sites) melted by vacancy pairs.
k non-intersecting paths (steps 0 or +1 per layer) with free ends, fugacity z per end.
Prediction (variational): rho = (2/pi) arccos(2^(x-1)), x = 2 log2(1/z)/N ; rho = k/M."""
import numpy as np, sys
def layer(v, M):
    # many-body transfer for non-intersecting 0/+1 steps on a ring; sum over carry into site 0
    out = np.zeros_like(v)
    idx = np.arange(1 << M)
    for c0 in (0, 1):
        # A[c] : amplitude with carry c entering current site; bits < i already new, bits >= i old
        A = [v.copy() if c0 == 0 else np.zeros_like(v), v.copy() if c0 == 1 else np.zeros_like(v)]
        for i in range(M):
            X0 = A[0].reshape(-1, 2, 1 << i); X1 = A[1].reshape(-1, 2, 1 << i)
            N0 = np.zeros_like(X0); N1 = np.zeros_like(X1)
            # old s_i = 0: new t_i = carry_in, carry_out = 0
            N0[:, 0] += X0[:, 0]; N0[:, 1] += X1[:, 0]
            # old s_i = 1, stays (needs carry_in = 0): t_i = 1, carry_out 0
            N0[:, 1] += X0[:, 1]
            # old s_i = 1, moves: t_i = carry_in, carry_out = 1
            N1[:, 0] += X0[:, 1]; N1[:, 1] += X1[:, 1]
            A = [N0.reshape(-1), N1.reshape(-1)]
        out += A[c0]          # carry leaving last site must equal carry entering site 0
    return out
def lnZ(M, N, z):
    pc = np.array([bin(i).count('1') for i in range(1 << M)])
    b = z ** pc.astype(float); v = b.copy(); ls = 0.0
    for _ in range(N):
        v = layer(v, M); s = v.max(); v /= s; ls += np.log(s)
    return ls + np.log(v @ b)
if __name__ == '__main__':
    M = int(sys.argv[1])
    for x in (0.2, 0.5, 0.8, 0.95):
        pred = (2 / np.pi) * np.arccos(2 ** (x - 1)); row = []
        for N in (8, 16, 32, 64):
            z = 2.0 ** (-x * N / 2); h = 1e-4
            rho = (lnZ(M, N, z * np.exp(h)) - lnZ(M, N, z * np.exp(-h))) / (2 * h) / (2 * M)
            row.append((N, rho))
        print(f'M {M} x {x}: predicted rho {pred:.4f} | ' + '  '.join(f'N={n}: {r:.4f}' for n, r in row), flush=True)
