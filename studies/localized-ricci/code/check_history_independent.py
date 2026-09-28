import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.chebyshev import chebvander
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import localized_trace as c
ROOT=Path(__file__).resolve().parents[1]
def main():
 bg=c.Background(.2);cfg=json.loads((ROOT/'inputs/r030_fine.json').read_text());x=np.array([-.943,-.711,-.397,-.063,.227,.583,.829,.953]);tau=x*5*cfg['sigma'];P=c.radial_proper_primitive;rc=cfg['center'];rs=np.array([brentq(lambda r:P(r)-P(rc)+t,.2,1.5,xtol=1e-14) for t in tau]);t=bg.base.times(rs);r,A,rp,*_=bg.at(t).T;errors=[]
 for l,i in [(0,3),(3,17),(80,85),(150,171)]:
  p=ROOT/'results/r030_fine_histories'/f'l{l:04d}.npz';z=np.load(p);nk=len(z['k']);kk=float(z['k'][i]);q=chebvander(x,int(z['degree']))@z['coefficients'];W=q[:,i];u=q[:,nk+i];ph=q[:,2*nk+i]+z['phase_offset'][i];f=np.exp(-1j*ph)/(r*np.sqrt(2*W));d=(u-1j*W-rp/r)*f
  r0,A0,rp0,rpp0,Ap0,m0,mp0=bg.at(0);w0=np.sqrt(kk*kk+A0*(l*(l+1)/r0**2+m0));f0=1/(r0*np.sqrt(2*w0))
  def rhs(s,y):
   rr,aa,dr,_,_,mm,_=bg.at(s);return [y[1],-2*dr/rr*y[1]-(kk*kk+aa*(l*(l+1)/rr**2+mm))*y[0]]
  sol=solve_ivp(rhs,[0,t[-1]],np.array([f0,-1j*w0*f0]),t_eval=t,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01);assert sol.success
  err=max(float(np.max(abs(f-sol.y[0]))/np.max(abs(f))),float(np.max(abs(d-sol.y[1]))/np.max(abs(d))))
  errors.append(dict(l=l,k=kk,relative_error=err))
 result=dict(independent_off_grid_modes=errors,all_pass=max(x['relative_error'] for x in errors)<1e-8,scope='Four bounded representative off-fit-point tests, not an error theorem over every possible history/window.')
 (ROOT/'verification/HISTORY_OFF_GRID_CHECKS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
