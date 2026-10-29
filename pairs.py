"""Fixed opposite-colour monomer pair in the z=0 Aztec diamond.
R(p,q) = Z_{pq}/Z_0.  One run per p: z_p = 1, z_q = eps on all opposite-colour cells, 0 elsewhere.
Then m_q = eps Z_pq / Z and Z_0/Z = 1 - m_p  =>  R(p,q) = m_q / (eps (1 - m_p)).
usage: python3 pairs.py N [boundary_only]
"""
import sys, time, numpy as np, aztec as A

def row(N, kp, ip, eps=1e-9):
    L = A.Aztec(N, 0).L
    cp = kp % 2
    zm = [np.full(n, eps if (k % 2) != cp else 0.0) for k, n in enumerate(L)]
    zm[kp][ip] = 1.0
    lnZ, m, l, r = A.Aztec(N, 0, zm).run()
    mp = m[kp][ip]
    R = {}
    for k, n in enumerate(L):
        if k % 2 == cp: continue
        for i in range(n):
            R[(k, i)] = m[k][i] / (eps * (1 - mp))
    return R

def boundary_cells(N, colour):
    L = A.Aztec(N, 0).L; out = []
    for k, n in enumerate(L):
        if k % 2 != colour: continue
        for i in range(n):
            if k in (0, 2 * N) or i in (0, n - 1): out.append((k, i))
    return out

if __name__ == '__main__':
    N = int(sys.argv[1]); bonly = len(sys.argv) > 2
    L = A.Aztec(N, 0).L
    ps = boundary_cells(N, 0) if bonly else [(k, i) for k, n in enumerate(L) if k % 2 == 0 for i in range(n)]
    best = (0, None, None); t = time.time(); S = 0.0
    allR = {}
    for p in ps:
        R = row(N, *p)
        q, v = max(R.items(), key=lambda kv: kv[1])
        allR[p] = R
        S += sum(R.values())
        if v > best[0]: best = (v, p, q)
    v, p, q = best
    print(f"N {N} {'boundary-p' if bonly else 'all-p'}  maxR {v:.6e} ln {np.log(v):.4f}  p {p} xy {A.cell_coords(N,*p)}  q {q} xy {A.cell_coords(N,*q)}  sumR(over scanned p) {S:.6e}  t {time.time()-t:.0f}s", flush=True)
    import pickle; pickle.dump(allR, open(f'az/pairs_N{N}{"_b" if bonly else ""}.pkl', 'wb'))
