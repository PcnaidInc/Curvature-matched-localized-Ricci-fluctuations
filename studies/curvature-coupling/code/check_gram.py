import os
for x in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[x]='1'
import numpy as np,sympy as s,json,sys
from pathlib import Path
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
import coupling as c
ROOT=Path(__file__).resolve().parents[1]
def main(source):
 cfg=json.loads((source/'inputs/r060_fine.json').read_text());bg=c.core.Background(.2);win=c.core.window(bg,cfg)
 r,A,rp,*_=bg.at(win['t']).T;x=win['tau']/(5*cfg['sigma']);qw=win['qw'];xx,w=leggauss(cfg['nk']);k=(xx+1)*cfg['K']/2;kw=w*cfg['K']/2;b=cfg['bproper']/np.sqrt(2/cfg['center']-1)
 weights=[kw[:,None]*kw[None,:]*np.exp(-b*b*(k[:,None]+sg*k[None,:])**2/2)for sg in(1,-1)]
 errors=[]
 for l,p in [(0,0),(3,4),(35,39),(80,83)]:
  f,d,_=c.load_hist(source,cfg['name'],l,x,r,rp);ff,dd,_=c.load_hist(source,cfg['name'],p,x,r,rp)
  D=d.T@((qw/A)[:,None]*dd);K=f.T@((qw/A)[:,None]*ff);R=f.T@((qw/r**2)[:,None]*ff);H=f.T@((qw*c.second_kernel_ratio(win['tau']/cfg['sigma'],cfg['sigma']))[:,None]*ff)
  lp=l*(l+1)+p*(p+1);g=c.gram_pair(D,K,R,H,k,weights,lp)
  for J in range(abs(l-p),min(l+p,32)+1,2):
   pred=c.to_basis(g,J);direct=np.zeros((3,3))
   for sg,wt in zip((1,-1),weights):
    B=D+sg*k[:,None]*k[None,:]*K-(lp-J*(J+1))*R/2
    E=D-(H+B)/6;v=[B,D,E]
    for a in range(3):
     for bb in range(3):direct[a,bb]+=np.sum(wt*np.real(v[a]*v[bb].conj()))
   scale=np.sqrt(np.outer(np.diag(direct),np.diag(direct)))
   errors.append(float(np.max(abs(direct-pred)/np.maximum(scale,1e-300))))
 # Fully independent 4D Minkowski temporal-Gaussian reference, no spatial smearing.
 u,mu=s.symbols('u mu',real=True);b0=-u*(1-u)*(1-mu);d0=-u*(1-u);e0=d0-(-1+b0)/6;v=[b0,d0,e0]
 exact=s.Matrix(3,3,lambda i,j:s.simplify(3*s.integrate(u*(1-u)*s.integrate(v[i]*v[j],(mu,-1,1)),(u,0,1))))
 # radial sum integral48, phase factor1/16 ->3. Coefficient multiplies 1/(pi^4 sigma^8).
 n=64;t,tw=leggauss(n);q=(t+1)*14/2;qw=tw*14/2;mm,mw=leggauss(20);numeric=np.zeros((3,3))
 for ki,wi in zip(q,qw):
  for kj,wj in zip(q,qw):
   B=ki*kj*(-1+mm);D=np.full_like(mm,-ki*kj);H=-(ki+kj)**2;E=D-(H+B)/6; vv=np.array([B,D,E]);
   numeric+=(wi*wj*ki*kj*np.exp(-(ki+kj)**2/2)/16)*(vv*mw)@vv.T
 ex=np.array(exact).astype(float);den=np.sqrt(np.outer(np.diag(ex),np.diag(ex)));flaterr=float(np.max(abs(numeric-ex)/den))
 out=dict(gram_direct_parent_normalized_max=max(errors),flat_exact_matrix_coefficient_over_pi4=[[str(exact[i,j])for j in range(3)]for i in range(3)],flat_numerical_error=flaterr,
    flat_reference_scope='Massless Minkowski vacuum; time-only Gaussian unit exponent width, no spatial average. Independent reference, not the interior observable.',
    all_pass=max(errors)<1e-8 and flaterr<1e-10)
 (ROOT/'verification/GRAM_CHECK.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));
 if not out['all_pass']:raise SystemExit(1)
if __name__=='__main__':main(Path(sys.argv[1]))
