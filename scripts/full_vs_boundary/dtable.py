"""D(N,z) = Q_un(full) - Q_un(boundary, exact) at the same (N,z), for every converged full-model point."""
import sys, glob, csv, pickle, os, numpy as np
from qb import Qb
R=''
rows=[]
for f in sorted(glob.glob(R+'taskA_results/summaries/*/summary.csv')+glob.glob(R+'taskB_results/summaries/*/summary.csv')):
    run=f.split('/')[-2]
    if run in ('B64b','diag32c','diag32h'): continue
    for r in csv.DictReader(open(f)):
        N=int(r['N']); z=float(r['z']); q=float(r['Q_un']); e=float(r['err'])
        if e<0.03 and z<=0.2: rows.append((N,z,q,e,run))
for N in (12,16,20):
    f=f'data/full_vs_boundary/full_exact_{N}.pkl'
    if os.path.exists(f):
        for z,(q,n) in pickle.load(open(f,'rb')).items(): rows.append((N,z,q,0.0,'exact'))
cache='data/full_vs_boundary/dtable.pkl'; D=pickle.load(open(cache,'rb')) if os.path.exists(cache) else {}
rows.sort()
for N,z,q,e,run in rows:
    key=(N,round(np.log(z),8),run)
    if key not in D:
        D[key]=(N,z,q,e,Qb(N,z),run); pickle.dump(D,open(cache,'wb'))
    N,z,q,e,qb,run=D[key]; x=2*np.log2(1/z)/N
    print(f'{N:3d} z={z:.4e} x={x:.4f} N*={x*N:5.1f} Qfull={q:8.4f}+-{e:.4f} Qbdry={qb:8.4f} D={q-qb:+.4f}  {run}',flush=True)
