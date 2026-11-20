import sys, numpy as np
nbz={0.0:0,0.0001:0,0.0003:0,0.001:0.0002,0.003:0.0007,0.01:0.00267,0.03:0.00918,0.1:0.03543}
def psi(fn, z):
    a=np.loadtxt(fn); h=a[:,2:]; P=h/h.sum(1,keepdims=True); phi=P.max(1)
    nb=nbz.get(z,0); pb=max(nb,(1-nb)/4)
    return np.mean((phi-pb)/(1-pb))
def blocks(fn):
    Q=[float(l.split()[2]) for l in open(fn) if l.startswith('B')]
    Q=np.array(Q); return Q.mean(), Q.std()/np.sqrt(len(Q)), len(Q)
if __name__=='__main__':
    pre=sys.argv[1]; z=float(sys.argv[2])
    m,e,k=blocks(pre+'.txt'); print(f"Q_top {m:.4f} +- {e:.4f} ({k} blocks)  Psi {psi(pre+'_hist.txt',z):.4f}")
