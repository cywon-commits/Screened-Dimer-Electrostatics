import numpy as np, tm_diag_b as T
T.RIN=1.0
for t in (0.3,0.5,0.7):
    pred=(2/np.pi)*np.arctan(t); row=[]
    for W in (6,8,10,12):
        z=t**(W/2); Lam,_,ne,no=T.densities(W,z); h=W//2
        Q=ne[:h].sum()-no[:h].sum(); nb=(ne[h-1]+no[h-1]+ne[h]+no[h])/4
        row.append((W,z,ne[0],Q,nb))
    print('full model  t=%.1f  predicted adsorbed charge per edge cell %.4f'%(t,pred))
    for W,z,a,Q,nb in row: print('   W %2d z %.2e  edge-site density %.4f   total adsorbed charge Q %.4f   mid-strip monomer density %.2e'%(W,z,a,Q,nb))
    Ws=np.array([r[0] for r in row[-3:]]); Qs=np.array([r[3] for r in row[-3:]]); print('   1/W extrapolation of Q: %.4f'%np.polyfit(1/Ws,Qs,1)[1],flush=True)
