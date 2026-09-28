"""Bounded independent checks; none rewrites a production stress/covariance."""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[n]='1'
from pathlib import Path
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp,quad
from scipy.special import eval_legendre
from sympy import Matrix,diag,symbols,simplify
from sympy.physics.wigner import wigner_3j
import joint_inward as j
import localized_trace as c
from localized_trace_fast import compact_fourier
ROOT=Path(__file__).resolve().parents[1]

def main():
 out={};bg=c.Background(.2)
 samples=[]
 for name in ['r060_main','r045_main','r030_main']:
  cfg=json.loads((ROOT/'inputs'/f'{name}.json').read_text());win=c.window(bg,cfg)
  rc=cfg['center'];lk=rc**1.5/48**.25
  samples.append(dict(center=rc,support=win['support'],time_fraction=cfg['sigma']/lk,longitudinal_fraction=cfg['bproper']/lk,angular_fraction=rc*cfg['angular_deltas'][-1]/lk,norm=float(win['qw'].sum()),min_radius=float(win['r'].min())))
 out['sampling']=samples
 out['fractions_relative_error']=float(np.max(abs(np.array([[q['time_fraction'],q['longitudinal_fraction'],q['angular_fraction']] for q in samples])/np.array([samples[0]['time_fraction'],samples[0]['longitudinal_fraction'],samples[0]['angular_fraction']])-1)))
 out['full_support_valid']=all(q['support'][0]>.2 and q['support'][1]<1.5 for q in samples)
 rs=np.linspace(.21,1.45,31);ts=bg.base.times(rs);vals=bg.at(ts)
 out['analytic_background_overlap_abs']=float(max(np.max(abs(vals[:,0]-rs)),np.max(abs(vals[:,1]-(2/rs-1)))))
 x,w=leggauss(240);errors=[]
 for l,p,J in [(0,0,0),(2,3,3),(19,24,13),(160,171,21),(192,191,31)]:
  a=c.wigner000_squared(l,p,J);b=float(wigner_3j(l,p,J,0,0,0))**2;errors.append(abs(a/b-1))
 out['wigner_symbolic_relative']=max(errors)
 errors=[]
 for l,p,J in [(2,3,3),(5,8,7),(20,24,12)]:
  errors.append(abs((w@(eval_legendre(l,x)*eval_legendre(p,x)*eval_legendre(J,x))/2)/c.wigner000_squared(l,p,J)-1))
 out['Gaunt_quadrature_relative']=max(errors)
 # Two exact field/Einstein tensor identities.
 v=Matrix(symbols('v0:4'));g=diag(-1,1,1,1);s=(v.T*g*v)[0];T=v*v.T-g*s/2;tr=sum(g[i,i]*T[i,i] for i in range(4))
 out['trace_reversal_exact']=simplify(T-g*tr/2-v*v.T)==Matrix.zeros(4)
 out['rho_from_joint_exact']=simplify(T[0,0]-(v[0]**2-tr/2))==0
 A,r,rp,k,kp,L,P,J=symbols('A r rp k kp L P J',nonzero=True);f,ff,d,dd=symbols('f ff d dd')
 fpp=-2*rp/r*d-(k*k+A*L/r**2)*f;gpp=-2*rp/r*dd-(kp*kp+A*P/r**2)*ff
 raw=d*dd/A+k*kp*f*ff/A-(L+P-J*(J+1))*f*ff/(2*r*r)
 bybox=(fpp*ff+2*d*dd+f*gpp+2*rp/r*(d*ff+f*dd))/(2*A)+(k+kp)**2*f*ff/(2*A)+J*(J+1)*f*ff/(2*r*r)
 out['KG_trace_identity_exact']=simplify(raw-bybox)==0
 # Independent complex field f, f' against nonlinear Gaussian width and phase.
 cfg=json.loads((ROOT/'inputs/r030_main.json').read_text());cfg['nt']=32;win=c.window(bg,cfg);t=win['t'];rr,aa,dr,*_=bg.at(t).T
 ferr=[];wronsk=[];fields=[]
 for ell,kk in [(0,.8),(3,2.7),(31,45.),(80,220.)]:
  ll=ell*(ell+1);r0,A0,rp0,rpp0,Ap0,m0,mp0=bg.at(0);w0=np.sqrt(kk*kk+A0*(ll/r0**2+m0));f0=1/(r0*np.sqrt(2*w0))
  def rhs(z,y):
   r,A,rp,rpp,Ap,m,mp=bg.at(z);return [y[1],-2*rp/r*y[1]-(kk*kk+A*(ll/r**2+m))*y[0]]
  sol=solve_ivp(rhs,[0,t[-1]],np.array([f0,-1j*w0*f0]),t_eval=t,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01)
  gs=solve_ivp(lambda z,y:bg.rhs(z,y,np.array([kk]),np.array([ll]),1),[0,t[-1]],np.zeros(3),t_eval=t,method='DOP853',rtol=1e-12,atol=1e-14,max_step=.01)
  assert sol.success and gs.success
  u,de,ph=gs.y;W=np.sqrt(kk*kk+aa*ll/rr**2)+de;f=np.exp(-1j*ph)/(rr*np.sqrt(2*W));d=(u-1j*W-dr/rr)*f
  ferr.append(max(float(np.max(abs(f-sol.y[0]))/np.max(abs(f))),float(np.max(abs(d-sol.y[1]))/np.max(abs(d)))))
  wronsk.append(float(np.max(abs(rr**2*(sol.y[0]*sol.y[1].conj()-sol.y[1]*sol.y[0].conj())-1j))))
 out['independent_complex_field_relative']=max(ferr);out['independent_Wronskian_absolute']=max(wronsk)
 # The polynomial moment contraction compared with direct pair arrays, including rho.
 rng=np.random.default_rng(9026);nt=19;nk=11;f=rng.normal(size=(nt,nk))+1j*rng.normal(size=(nt,nk));d=rng.normal(size=(nt,nk))+1j*rng.normal(size=(nt,nk));ff=rng.normal(size=(nt,nk))+1j*rng.normal(size=(nt,nk));dd=rng.normal(size=(nt,nk))+1j*rng.normal(size=(nt,nk))
 rr=np.linspace(.4,.7,nt);aa=2/rr-1;qw=np.ones(nt)/nt;k=np.linspace(.1,8,nk);kw=np.ones(nk)/nk;b=.2;l,p,J=3,5,4
 weights=[kw[:,None]*kw[None,:]*np.exp(-b*b*(k[:,None]+sg*k[None,:])**2/2) for sg in (1,-1)]
 vv=j.physical_moments(f,d,ff,dd,rr,aa,qw,k,weights,l*(l+1)+p*(p+1));mat=j.evaluate(vv,J);direct=np.zeros((2,2));energy=0.
 for ii in range(nk):
  for jj in range(nk):
   for sg,weight in zip((1,-1),weights):
    D=np.sum(qw*d[:,ii]*dd[:,jj]/aa)
    B=np.sum(qw*(d[:,ii]*dd[:,jj]/aa+sg*k[ii]*k[jj]*f[:,ii]*ff[:,jj]/aa-(l*(l+1)+p*(p+1)-J*(J+1))/2*f[:,ii]*ff[:,jj]/rr**2))
    vec=np.array([B,D]);direct+=weight[ii,jj]*np.real(vec[:,None]*vec.conj()[None]);energy+=weight[ii,jj]*abs(D-B/2)**2
 out['joint_vs_direct_pairs_relative']=float(np.max(abs(mat-direct))/np.max(abs(direct)))
 out['energy_direct_combination_relative']=float(abs(energy/(np.array([-.5,1])@mat@np.array([-.5,1]))-1))
 ph1=np.exp(1j*rng.normal(size=nk));ph2=np.exp(1j*rng.normal(size=nk));vv2=j.physical_moments(f*ph1,d*ph1,ff*ph2,dd*ph2,rr,aa,qw,k,weights,l*(l+1)+p*(p+1))
 out['phase_invariance_relative']=float(np.max(abs(j.evaluate(vv2,J)-mat))/np.max(abs(mat)))
 # Time-only Minkowski Wick references; not an identical spatial window.
 x,w=leggauss(160);k=8*(x+1);kw=8*w;inte=np.sum(kw[:,None]*kw[None,:]*k[:,None]**3*k[None,:]**3*np.exp(-.5*(k[:,None]+k[None,:])**2))
 out['flat_trace_reference_relative']=float(abs((8/3)*inte/(16*np.pi**4)/(2/(35*np.pi**4))-1))
 out['flat_time_projection_reference_relative']=float(abs(inte/(8*np.pi**4)/(3/(70*np.pi**4))-1))
 # Finite product-space reference Fourier vs independent time quadrature at different radii.
 errs=[]
 for rc in [.6,.3]:
  cfg={'center':rc,'K':60,'L':6,'sigma':.02*(rc/.6)**1.5};fr=j.FrozenReference(cfg)
  x,w=leggauss(384);z=5*x;ww=5*w*np.exp(-z*z)*c.cutoff(z);ww/=ww.sum();ts=np.linspace(0,12,87);direct=np.cos(ts[:,None]*z[None])@ww
  errs.append(float(np.max(abs(compact_fourier(ts)-direct))))
 out['compact_Fourier_vs_quadrature_absolute']=max(errs)
 # Harmonic normalization / full-sphere monopole identity.
 out['monopole_identity_relative']=max(abs((2*l+1)*c.wigner000_squared(l,l,0)-1) for l in (0,3,31,127,191))
 numeric={k:v for k,v in out.items() if k.endswith(('relative','absolute','abs','error'))}
 out['all_pass']=bool(all(v<1e-8 for v in numeric.values()) and out['full_support_valid'] and out['trace_reversal_exact'] and out['KG_trace_identity_exact'] and out['rho_from_joint_exact'])
 (ROOT/'verification/INDEPENDENT_CHECKS.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()
