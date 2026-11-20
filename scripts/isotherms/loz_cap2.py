import numpy as np
from loz_cap import lnZ
from ckp import L
s=lambda r:(2*L(np.pi*r/2)+L(np.pi*(1-r)))/np.pi
print('free energy per site  lnZ/(MN)  vs  s(rho*) - rho* x ln2')
for M in (12,16):
  for x in (0.2,0.5,0.8,0.95):
    r=(2/np.pi)*np.arccos(2**(x-1)); pred=s(r)-r*x*np.log(2)
    vals=[lnZ(M,N,2.0**(-x*N/2))/(M*N) for N in ((16,32,64) if M==12 else (16,32))]
    print(f'M {M} x {x}: pred {pred:.5f} | '+' '.join(f'{v:.5f}' for v in vals),flush=True)
