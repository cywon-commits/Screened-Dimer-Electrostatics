import numpy as np, sys, time
from img_test import prof
W = int(sys.argv[1]); zb_fac=float(sys.argv[2])
for z in [float(v) for v in sys.argv[3:]]:
    zs=np.full(W,z); zs[0]=zs[-1]=z*zb_fac
    t=time.time(); n = prof(W, zs); r=n/n[W//2]
    env=[(r[2*k]+r[2*k+1])/2 for k in range(W//4)]
    xi=1/np.sqrt(3.04*n[W//2])
    print(f"W {W} z {z} zb/z {zb_fac} n_mid {n[W//2]:.5f} xi~{xi:.1f} | r(d) " + ' '.join(f'{v:.3f}' for v in r[:W//2]) + " | env(pairs) "+' '.join(f'{v:.3f}' for v in env), f'({time.time()-t:.0f}s)', flush=True)
