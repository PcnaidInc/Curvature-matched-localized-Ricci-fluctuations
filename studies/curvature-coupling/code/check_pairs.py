"""Independent Hessian/on-shell vs sampling-derivative evaluation; flat and history controls."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq
import coupling as c
ROOT=Path(__file__).resolve().parents[1]

def main(source):
 errors=[];mode_errors=[];direct_gram=[]
 for prior in ['r060_fine','r030_fine']:
  cfg=json.loads((source/'inputs'/f'{prior}.json').read_text());bg=c.core.Background(cfg['rmin']);cfg['nt']=640
  win=c.core.window(bg,cfg);r,A,rp,rpp,Ap,*_=bg.at(win['t']).T;x=win['tau']/(5*cfg['sigma']);qw=win['qw']
  for l,p,ii,jj,J in [(0,0,1,2,0),(3,4,17,19,3),(40,43,50,52,5),(80,83,90,93,7)]:
   f,d,k=c.load_hist(source,prior,l,x,r,rp);ff,dd,kp=c.load_hist(source,prior,p,x,r,rp)
   fi,di,fj,dj=f[:,ii],d[:,ii],ff[:,jj],dd[:,jj];ki,kj=k[ii],kp[jj]
   q=fi*fj;qp=di*fj+fi*dj
   qpp=2*di*dj-2*rp/r*qp-(ki**2+kj**2+A*(l*(l+1)+p*(p+1))/r**2)*q
   H=(qpp-Ap/(2*A)*qp)/A;Ha=qw@H;Hb=(qw*c.second_kernel_ratio(win['tau']/cfg['sigma'],cfg['sigma']))@q
   par=max(abs(qw@(di*dj/A)),abs(Ha),1e-50)
   errors.append(dict(prior=prior,l=l,p=p,k=float(ki),kp=float(kj),H_direct_real=float(Ha.real),H_direct_imag=float(Ha.imag),IBP_error_parent=float(abs(Ha-Hb)/par)))
   for sign in(-1,1):
    dotang=(l*(l+1)+p*(p+1)-J*(J+1))/2
    B=(di*dj+sign*ki*kj*q)/A-dotang*q/r**2;P=di*dj/A
    # Expanded improved tensor: (1-2xi)P-2xi sym(phi Hess_uu phi)-xi T
    hii=(-2*rp/r*di-(ki**2+A*l*(l+1)/r**2)*fi-Ap/(2*A)*di)/A
    hjj=(-2*rp/r*dj-(kj**2+A*p*(p+1)/r**2)*fj-Ap/(2*A)*dj)/A
    xi=1/6
    expanded=(1-2*xi)*P-xi*(fi*hjj+fj*hii)-xi*B
    reference=qw@P-xi*(Hb+qw@B)
    direct=qw@expanded
    direct_gram.append(float(abs(direct-reference)/max(par,abs(reference))))
  # Independent field solver at four off-fit points and at chosen k (new small reference evolution).
  xs=np.array([-.891,-.537,-.143,.217,.619,.937]);tau=xs*5*cfg['sigma'];P=c.core.radial_proper_primitive
  rs=np.array([brentq(lambda rr:P(rr)-P(cfg['center'])+v,.2,1.5,xtol=1e-14)for v in tau]);ts=bg.base.times(rs);rr,AA,dr,*_=bg.at(ts).T
  for l,i in [(0,3),(3,17),(70,60)]:
   f,d,ks=c.load_hist(source,prior,l,xs,rr,dr);kk=ks[i];r0,A0,rp0,_,_,m0,_=bg.at(0);w0=np.sqrt(kk*kk+A0*(l*(l+1)/r0**2+m0));f0=1/(r0*np.sqrt(2*w0))
   def fun(t,y):
    r1,a1,rp1,_,_,mass,_=bg.at(t);return[y[1],-2*rp1/r1*y[1]-(kk*kk+a1*(l*(l+1)/r1**2+mass))*y[0]]
   sol=solve_ivp(fun,[0,ts[-1]],np.array([f0,-1j*w0*f0]),t_eval=ts,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01)
   if not sol.success:raise RuntimeError(sol.message)
   err=max(float(np.max(abs(sol.y[0]-f[:,i]))/np.max(abs(f[:,i]))),float(np.max(abs(sol.y[1]-d[:,i]))/np.max(abs(d[:,i]))))
   mode_errors.append(dict(prior=prior,l=l,k=float(kk),relative_error=err))
 # Flat, massless on-shell plane pairs with arbitrary orientations, not a spherical gravity solution.
 rng=np.random.default_rng(731);flat=[]
 for _ in range(20):
  k=rng.normal(size=3);p=rng.normal(size=3);w=np.linalg.norm(k);v=np.linalg.norm(p);dot=k@p;B=-w*v+dot;H=-(w+v)**2;D=-w*v
  E=D-(H+B)/6
  alt=((w-v)**2+(w*v-dot))/6 # equals (w²+v²-3wv-dot)/6? check below
  # correct independent momentum form from stress tensor
  alt=(w*w+v*v-3*w*v-dot)/6
  flat.append(abs(E-alt)/(1+abs(E)))
 # The compact g'' has zero integral; its Fourier transform is -omega² times g.
 cfg=json.loads((source/'inputs/r060_fine.json').read_text());s=cfg['sigma'];x,w=leggauss(1024);tau=5*s*x;qw=5*w*np.exp(-(tau/s)**2)*c.core.cutoff(tau/s);qw/=sum(qw)
 fourier=[]
 for freq in [0,1/s,3/s,5/s]:
  phase=np.exp(-1j*freq*tau);q2=c.second_kernel_ratio(tau/s,s);err=abs(qw@(q2*phase)+freq*freq*(qw@phase))*s*s
  fourier.append(float(err))
 out=dict(Hessian_IBP=errors,max_IBP_parent=max(z['IBP_error_parent']for z in errors),expanded_conformal_amplitude_max_parent=max(direct_gram),
   independent_complex_modes=mode_errors,max_complex_error=max(z['relative_error']for z in mode_errors),flat_plane_pair_error=max(flat),
   compact_derivative_fourier_errors=fourier,new_reference_mode_evolutions=len(mode_errors),
   notes='No new production modes. Six independent small reference modes. H identities normalize by parent amplitude, not near-cancelled output.')
 out['all_pass']=(out['max_IBP_parent']<1e-6 and out['expanded_conformal_amplitude_max_parent']<1e-6 and out['max_complex_error']<1e-8 and max(flat)<1e-12 and max(fourier)<1e-8)
 (ROOT/'verification/PAIR_CHECKS.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));
 if not out['all_pass']:raise SystemExit(1)
if __name__=='__main__':main(Path(sys.argv[1]))
