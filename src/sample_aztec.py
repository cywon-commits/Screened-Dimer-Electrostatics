"""Equilibrium sample of the Aztec diamond with boundary vacancies (fugacity z per vacancy).
Vacancy set: Metropolis on subsets S of the 4N boundary cells with weight z^|S| M(G\\S), M = number of perfect matchings (Kasteleyn determinant).
Tiling given S: exact sequential sampling from the inverse Kasteleyn matrix."""
import sys, pickle, numpy as np
def geometry(N):
    cells = [(x, y) for y in range(2 * N) for x in range(2 * N) if abs(x - N + .5) + abs(y - N + .5) <= N]
    edge = [[(N - 1 - p, p), (N + p, p), (N + p, 2 * N - 1 - p), (N - 1 - p, 2 * N - 1 - p)][e] for e in range(4) for p in range(N)]
    return cells, edge
def kast(cells, S=()):
    S = set(S); B = [c for c in cells if (c[0] + c[1]) % 2 == 0 and c not in S]; W = [c for c in cells if (c[0] + c[1]) % 2 == 1 and c not in S]
    iw = {c: i for i, c in enumerate(W)}; K = np.zeros((len(B), len(W)), complex)
    for i, (x, y) in enumerate(B):
        for dx, dy, w in ((1, 0, 1), (-1, 0, 1), (0, 1, 1j), (0, -1, 1j)):
            j = iw.get((x + dx, y + dy));
            if j is not None: K[i, j] = w
    return B, W, K
def lnM(cells, S):
    B, W, K = kast(cells, S)
    return -np.inf if len(B) != len(W) else np.linalg.slogdet(K)[1]
def run(N, x, nstep, rng, burn=0, every=0):
    z = 2.0 ** (-x * N / 2); cells, edge = geometry(N)
    ev = [c for c in edge if sum(c) % 2 == 0]; od = [c for c in edge if sum(c) % 2 == 1]
    S = set(); l = lnM(cells, S); occ = np.zeros(4 * N); nacc = ns = 0; idx = {c: i for i, c in enumerate(edge)}
    for t in range(nstep):
        a = ev[rng.integers(len(ev))]; b = od[rng.integers(len(od))]
        if (a in S) == (b in S):                      # create or annihilate a +/- pair
            S2 = S ^ {a, b}; dz = 2 * np.log(z) * (1 if a not in S else -1)
        else:                                         # move a vacancy within its sublattice
            L = ev if rng.random() < .5 else od; a = L[rng.integers(len(L))]; b = L[rng.integers(len(L))]
            if (a in S) == (b in S): S2 = None
            else: S2 = S ^ {a, b}; dz = 0.0
        if S2 is not None:
            l2 = lnM(cells, S2)
            if np.log(rng.random()) < l2 - l + dz: S, l = S2, l2; nacc += 1
        if t >= burn:
            ns += 1
            for c in S: occ[idx[c]] += 1
    return cells, S, occ / max(ns, 1), nacc / nstep
def tiling(cells, S, rng):
    cells = [c for c in cells if c not in S]; dom = []
    while cells:
        B, W, K = kast(cells); Ki = np.linalg.inv(K); p = np.abs(K[0] * Ki[:, 0]); p /= p.sum(); j = rng.choice(len(W), p=p)
        dom.append((B[0], W[j])); cells = [c for c in cells if c != B[0] and c != W[j]]
    return dom
if __name__ == "__main__":
    mode = sys.argv[1]; rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    if mode == "check":                               # compare the edge profile with the exact boundary-Pfaffian values
        N, x = 24, 0.5; cells, S, occ, acc = run(N, x, 60000, rng, burn=5000)
        ref = pickle.load(open("../data/smallx.pkl", "rb"))[(N, x)][1]; mc = occ.reshape(4, N).mean(0)
        print("acceptance %.2f" % acc); print("exact ", np.round(ref[:N], 3)); print("MC    ", np.round(mc, 3)); print("max |diff| %.3f, mean total %.3f vs %.3f" % (abs(mc - ref[:N]).max(), mc.sum(), ref[:N].sum()))
    else:
        N, x = 16, 0.5; cells, S, occ, acc = run(N, x, 40000, rng)
        pickle.dump(dict(N=N, x=x, vac=sorted(S), dominoes=tiling(cells, S, rng)), open("../data/fig1_sample.pkl", "wb")); print("N", N, "x", x, "vacancies", sorted(S), "acc %.2f" % acc)
