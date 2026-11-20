import numpy as np, pickle
from scipy.optimize import brentq
from exact_sol import solve
out=pickle.load(open('az/smallx.pkl','rb')); Ns=(40,48,56)
old={0.6:0.799,0.7:0.867,0.8:0.924,0.9:0.970}
print('  x     exact solution sigma    Pfaffian extrap (lin / quad)     d*')
res={}
for x in (0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9):
    lj=brentq(lambda l: solve(1+np.exp(l))['x']-x,-9,11,xtol=1e-10); s=solve(1+np.exp(lj)); res[x]=s
    if (40,x) in out:
        v=np.array([out[(N,x)][0] for N in Ns]); iN=1/np.array(Ns,float)
        ex='%.4f / %.4f'%(np.polyfit(iN,v,1)[-1],np.polyfit(1/np.array((32,40,48,56.)),np.array([out[(N,x)][0] for N in (32,40,48,56)]),2)[-1])
    else: ex='%.3f (earlier fit)'%old[x]
    print('%.2f    %.5f               %s        %.4f'%(x,s['sigma'],ex,s['dstar']),flush=True)
pickle.dump({x:(r['sigma'],r['dstar'],r['u'],r['h']) for x,r in res.items()},open('az/exact_sol.pkl','wb'))
