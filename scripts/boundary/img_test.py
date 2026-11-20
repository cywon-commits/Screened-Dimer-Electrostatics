import numpy as np, sys, time
from tm_open import make_T, top
def prof(W, zs):
    T = make_T(W, zs, False); vals, V = top(W, T); L0 = vals[0]; psi = V[:, 0]; Tpsi = T(psi)
    n = np.empty(W)
    for x in range(W):
        z0 = np.array(zs, float); z0[x] = 0; n[x] = psi @ (Tpsi - T(psi, z0)) / L0
    return n
if __name__=="__main__":
  W = int(sys.argv[1])
  for z in [float(v) for v in sys.argv[2:]]:
    t=time.time(); n = prof(W, np.full(W, z))
    print(f"W {W} z {z}  n_mid {n[W//2]:.6f}  ratio n(d)/n_mid d=0..{W//2-1}: " + ' '.join(f'{v:.4f}' for v in n[:W//2]/n[W//2]), f'({time.time()-t:.0f}s)', flush=True)
