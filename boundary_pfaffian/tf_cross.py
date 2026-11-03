import numpy as np, pickle
from scipy.optimize import minimize
import tf_solve as T
from xkernel import E
def solve(x, M=400, K=None, ufun='cross'):
    h=1.0/M; d=-0.5+h*(np.arange(M)+0.5)
    if ufun=='cross': u=T.H(0.5+d)-np.log(2)/2
    else: u=np.array([E(min(0.5+di,0.5-di),min(0.5+di,0.5-di))[0]/2 for di in d])
    f=u-x*np.log(2)/2
    D=d[:,None]-d[None,:]
    L=np.log(np.abs(D)+np.eye(M)); np.fill_diagonal(L,np.log(h/2)-1)
    G=(L+K(D)/2)*h*h
    obj=lambda r: -(f@r*h+0.5*r@G@r); grad=lambda r: -(f*h+G@r)
    res=minimize(obj,np.clip(f,0,None),jac=grad,bounds=[(0,1)]*M,method='L-BFGS-B',options=dict(maxiter=5000,ftol=1e-14,gtol=1e-12))
    return d,res.x,1-2*res.x.sum()*h
if __name__=='__main__':
    K=T.Kinf(); ex=dict((round(float(a),2),b) for a,b,*_ in pickle.load(open('az/tf_compare.pkl','rb')))
    for x in (0.3,0.4,0.5,0.6,0.8):
        s1=solve(x,K=K)[2]; s2=solve(x,K=K,ufun='full')[2]
        print(f'x={x}: exact {ex[x]:.3f}  TF(crossing u) {s1:.3f}  TF(crossover u) {s2:.3f}')
