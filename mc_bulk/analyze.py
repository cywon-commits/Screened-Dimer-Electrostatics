import sys, glob, re
import numpy as np
from scipy.optimize import curve_fit

def load(fn):
    L = int(re.search(r'_L(\d+)', fn).group(1)); z = float(re.search(r'z([\d.]+)_', fn).group(1))
    n, Cq, Cd, Sq, D = [], [], [], [], []
    for line in open(fn):
        t = line.split()
        if not t or t[0] == '#': continue
        if t[0] == 'B': n.append(float(t[3])); D.append(float(t[5]))
        elif t[0] == 'Cq': Cq.append([float(v) for v in t[1:]])
        elif t[0] == 'Cd': Cd.append([float(v) for v in t[1:]])
        elif t[0] == 'Sq': Sq.append([float(v) for v in t[1:]])
    k = min(len(n), len(Cq), len(Cd), len(Sq))
    return L, z, np.array(n[:k]), np.array(D[:k]), np.array(Cq[:k]), np.array(Cd[:k]), np.array(Sq[:k])

RMAX = {}
def fitmass(C, L, r0, key=None):
    H = L // 2
    rm = RMAX.get(key, H)
    if rm - r0 < 3: return np.nan
    r = np.arange(r0, rm + 1); y = C[r0:rm + 1]
    if np.any(y <= 0): return np.nan
    f = lambda r, lnA, m: lnA + np.log(np.cosh(m * (r - H)))
    m0 = max(0.05, np.log(y[0] / y[1])) if y[1] > 0 else 0.2
    try:
        p, _ = curve_fit(f, r, np.log(y), p0=[np.log(y[-1]), m0], maxfev=20000)
        return p[1]
    except Exception:
        return np.nan

def jack(blocks, func):
    B = len(blocks[0]); full = func(*[b.mean(0) for b in blocks])
    vals = np.array([func(*[np.delete(b, i, 0).mean(0) for b in blocks]) for i in range(B)])
    err = np.sqrt((B - 1) * np.nanvar(vals, axis=0))
    return full, err

print(f"{'z':>7} {'L':>4} {'n':>10} {'r0':>3} {'m1':>14} {'m1^2/n':>14} {'m2':>14} {'m2/m1':>12} {'pi*S/k^2(k1..3)':>24}")
for fn in sorted(glob.glob(sys.argv[1] if len(sys.argv)>1 else 'runs/z*_L*.txt'), key=lambda f: -float(re.search(r'z([\d.]+)_', f).group(1))):
    L, z, n, D, Cq, Cd, Sq = load(fn)
    if len(n) < 4: continue
    H = L // 2
    B=len(n)
    for key,arr in (('q',-Cq),('d',Cd-D[:,None]**2)):
        mu=arr.mean(0); er=arr.std(0)/np.sqrt(B)
        bad=np.where(mu[1:]<4*er[1:])[0]
        RMAX[key]=(bad[0] if len(bad) else H)
    xi_est = 1 / np.sqrt(3.04 * n.mean())
    print(fn, 'rmax_q',RMAX['q'],'rmax_d',RMAX['d'])
    for fac in (1.0, 1.5, 2.0, 2.5):
        r0 = max(2, int(round(fac * xi_est)))
        
        def m1f(Cq_): return fitmass(-Cq_, L, r0, 'q')
        def ratio(n_, Cq_): return fitmass(-Cq_, L, r0, 'q') ** 2 / n_
        def m2f(Cd_, D_): return fitmass(Cd_ - D_ ** 2, L, r0, 'd')
        def m21(Cd_, D_, Cq_): return fitmass(Cd_ - D_ ** 2, L, r0, 'd') / fitmass(-Cq_, L, r0, 'q')
        m1, dm1 = jack([Cq], m1f)
        R, dR = jack([n, Cq], ratio)
        m2, dm2 = jack([Cd, D], m2f)
        q, dq = jack([Cd, D, Cq], m21)
        ks = 2 * np.pi * np.arange(1, 4) / L; khat2 = 4 * np.sin(ks / 2) ** 2
        sl = np.pi * Sq.mean(0)[:3] / khat2
        print(f"{z:7.3f} {L:4d} {n.mean():10.6f} {r0:3d} {m1:8.4f}({dm1:.4f}) {R:8.4f}({dR:.4f}) {m2:8.4f}({dm2:.4f}) {q:6.3f}({dq:.3f}) " + ' '.join(f'{v:.3f}' for v in sl))
print("predictions: m1^2/n = 7 sin(pi/7) = %.4f (DH pi=%.4f), m2/m1 = 2cos(pi/14) = %.4f, pi*S/k^2 -> 1" % (7*np.sin(np.pi/7), np.pi, 2*np.cos(np.pi/14)))
