import numpy as np, sys, io, contextlib, warnings
warnings.filterwarnings('ignore')
sys.argv=['v',sys.argv[1]]
with contextlib.redirect_stdout(io.StringIO()): import verify_sol as V
x=V.x
print('x =',x)
# diagonal ray: u=v and h_u=h_v
m=0
for r in (0.2,0.5,0.8,0.95):
    zt=r*np.exp(-1j*np.pi/4*(1-1e-9)); f=V.fields(zt,1); m=max(m,abs(f[0]-f[1]),abs(f[2]-f[3]))
print('  diagonal SO: max(|u-v|,|h_u-h_v|) = %.1e'%m)
# arc (frozen boundary): slopes -> vertex (-1,-1)
for eps in (1e-2,1e-3,1e-4):
    f=V.fields((1-eps)*np.exp(-1j*np.pi/8),1); print('  near the arc (1-|zeta|=%.0e): h = (%.4f, %.4f)'%(eps,f[2],f[3]))
# antimirror MO: u=1, h_v=0 ; centre O: h=0
f=V.fields(0.5*V.solve(V.JM)['rho_MO'][-1]-1e-9j,1); print('  on MO: u-1 = %.1e, h_v = %.1e'%(f[0]-1,f[3]))
f=V.fields(1e-3*np.exp(-1j*np.pi/8),1); print('  near O: (u,v) = (%.4f, %.4f), h = (%.4f, %.4f)'%tuple(f[:4]))
