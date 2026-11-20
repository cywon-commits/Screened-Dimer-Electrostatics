import numpy as np, re, glob, pickle
from scipy.optimize import brentq
from exact_sol import solve
d={}
for f in glob.glob('data/bx_*.log'):
    for l in open(f):
        m=re.match(r'N (\d+) x ([\d.]+) sigma ([\d.]+)',l)
        if m: d[(int(m[1]),float(m[2]))]=float(m[3])
pickle.dump(d,open('data/bx_exact.pkl','wb'))
print('  x    N used          Pfaffian N->inf (quadratic in 1/N)    exact solution    difference')
for x in (0.1,0.2,0.3,0.5,0.7,0.9):
    Ns=sorted(N for (N,xx) in d if xx==x and N>=64); v=np.array([d[(N,x)] for N in Ns]); iN=1/np.array(Ns,float)
    c=np.polyfit(iN[-3:],v[-3:],2)[-1]
    s=solve(1+np.exp(brentq(lambda l: solve(1+np.exp(l))['x']-x,-12,11,xtol=1e-12)),m=8000)['sigma']
    print('%.1f   %-14s  %.6f                           %.6f       %+.1e'%(x,Ns,c,s,c-s))
