#!/usr/bin/env python3
"""z ladder for aztec_pt with a target neighbour swap acceptance.

Swap acceptance between ln z values a distance D apart is acc ~ erfc(D * sigma_n / 2), where
sigma_n^2 = Var(n) = d<n>/d ln z (n = number of monomers); checked on the N=32 run to ~0.01.
Var(n) is taken from the measured curve of the largest reference size M <= N (ladder_ref_N<M>.txt, written
by `ladder.py ref OUTDIR M`) scaled extensively by the number of cells; below z=1e-2 (where N < xi and the fluctuations grow faster than the area) the scale factor is
squared, to stay on the safe side.  Every interval between required z points is cut into the fewest equal-acceptance pieces
whose predicted acceptance is >= target.

usage: ladder.py N ZMIN ZMAX TARGET   -> prints the ladder (ascending)
       ladder.py ref OUTDIR M          -> writes ladder_ref_N<M>.txt from finished N=M runs in OUTDIR
"""
import sys, os, glob, re, numpy as np
from scipy.special import erfcinv
REQ = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1]

def ncells(N): return 2 * N * (N + 1)

HERE = os.path.dirname(os.path.abspath(__file__))

def make_ref(outdir, M):
    d = {}
    for f in glob.glob(os.path.join(outdir, f'N{M}_s*_r*.txt')):
        if f.endswith('_hist.txt'): continue
        z = float(open(f).readline().split()[-1])
        d.setdefault(z, []).extend(float(l.split()[4]) for l in open(f) if l[0] == 'B')
    z = np.array(sorted(d)); n = np.array([np.mean(d[k]) for k in z]) * ncells(M)
    np.savetxt(os.path.join(HERE, f'ladder_ref_N{M}.txt'), np.c_[z, n, np.gradient(n, np.log(z))],
               header=f'z <n_total> Var(n)=d<n>/dlnz  (N={M}, NC={ncells(M)})')

def ladder(N, zmin, zmax, target):
    refs = sorted(int(re.search(r'N(\d+)', os.path.basename(f))[1]) for f in glob.glob(os.path.join(HERE, 'ladder_ref_N*.txt')))
    M = max([m for m in refs if m <= N] or refs[:1])
    z, n, var = np.loadtxt(os.path.join(HERE, f'ladder_ref_N{M}.txt'), unpack=True)
    lz = np.log(z); scale = ncells(N) / ncells(M)
    def sig(l):  # flat extrapolation outside the reference range; extra factor `scale` below z=1e-2
        w = np.clip((np.log(1e-2) - l) / np.log(3), 0, 1)
        return np.sqrt(scale ** (1 + w) * np.interp(l, lz, var))
    step = 2 * erfcinv(target)  # required int sigma dlnz per interval
    anchors = sorted({zmin, zmax} | {x for x in REQ if zmin < x < zmax})
    out = [anchors[0]]
    for a, b in zip(anchors[:-1], anchors[1:]):
        g = np.linspace(np.log(a), np.log(b), 2001); s = sig(g)
        cum = np.concatenate([[0], np.cumsum(0.5 * (s[1:] + s[:-1]) * np.diff(g))])
        k = max(1, int(np.ceil(cum[-1] / step - 1e-9)))
        out += [float(np.exp(np.interp(cum[-1] * j / k, cum, g))) for j in range(1, k)] + [b]
    return out

if __name__ == '__main__':
    if sys.argv[1] == 'ref': make_ref(sys.argv[2], int(sys.argv[3])); sys.exit()
    N, a, b, t = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    print(' '.join(f'{x:.6g}' for x in ladder(N, a, b, t)))
