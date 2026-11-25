"""Working figure: (a) edge free energy dg(rho) from the one-sided strip, (b) D = Q_un(full) - Q_un(boundary) against 1/N
with the predicted limits D_inf(x) at 1/N = 0.  Also writes dg.csv and D_table.csv.
Run from the repository root: PYTHONPATH=src:scripts/full_vs_boundary python3 scripts/full_vs_boundary/fig_full_vs_boundary.py"""
import pickle, csv, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import onesided_tab as O
P='data/full_vs_boundary/'
rs=np.arange(0.05,0.4001,0.01); Ws=sorted({k[0] for k in O.d if (k[0],'one') in O.d})
tab={W:O.dg(W,rs) for W in Ws}
with open(P+'dg.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['rho']+[f'dg_W{W}' for W in Ws])
    for i,r in enumerate(rs): w.writerow([f'{r:.2f}']+['' if np.isnan(tab[W][i]) else f'{tab[W][i]:.6f}' for W in Ws])
D=pickle.load(open(P+'dtable.pkl','rb')); rows=sorted(D.values())
with open(P+'D_table.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['N','z','x','Q_un_full','err','Q_un_boundary','D','source'])
    for N,z,q,e,qb,run in rows: w.writerow([N,f'{z:.10g}',f'{2*np.log2(1/z)/N:.5f}',f'{q:.5f}',f'{e:.5f}',f'{qb:.5f}',f'{q-qb:.5f}',run])
pred={0.7:-0.3882,0.5:-0.8504,0.4152:-1.1554,0.3:-1.7798}          # predict.py 16 0.4 ...
extra={0.4152:[(64,-1.467,0.038)]}   # N = 64, z = 1e-4: Q_un = 18.475 +- 0.038 (taskA_results/reports/interim4_N64pl.md), Q_un(boundary) = 19.942
fig,ax=plt.subplots(1,2,figsize=(9.6,3.9))
a=ax[0]
for W in Ws: a.plot(rs,tab[W],label=f'W = {W}',lw=1)
t=np.tan(np.pi*rs/2); a.plot(rs,0.545*t*t,'k--',lw=0.8,label=r'$0.545\,\tan^2(\pi\rho/2)$')
a.set_xlabel(r'edge vacancy density $\rho$'); a.set_ylabel(r'$\delta g(\rho)$'); a.legend(fontsize=7); a.set_title('(a) edge free energy, one-sided strip',fontsize=9)
a=ax[1]; col={0.7:'C0',0.5:'C1',0.4152:'C2',0.3:'C3'}
for x in pred:
    pts=[(N,q-qb,e) for N,z,q,e,qb,run in rows if abs(2*np.log2(1/z)/N-x)<2e-4]+extra.get(x,[])
    best={}
    for N,dv,e in pts:
        if N not in best or e<best[N][1]: best[N]=(dv,e)
    Ns=np.array(sorted(best)); a.errorbar(1/Ns,[best[N][0] for N in Ns],[best[N][1] for N in Ns],fmt='o-',ms=3.5,lw=0.8,color=col[x],label=f'x = {x}')
    a.plot([0],[pred[x]],marker='*',ms=10,color=col[x],clip_on=False)
a.set_xlim(0,0.09); a.set_xlabel('1/N'); a.set_ylabel(r'$Q_{\rm un}^{\rm full}-Q_{\rm un}^{\rm boundary}$'); a.legend(fontsize=7); a.set_title('(b) data and predicted limits (stars)',fontsize=9)
plt.tight_layout(); plt.savefig('figures/working/full_vs_boundary.png',dpi=150)
