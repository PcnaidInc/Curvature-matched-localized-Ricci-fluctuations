"""New angularly localized trace covariance in the explicitly prepared scalar state.
Harmonic Gaunt factors sum all m,m' exactly. Mode integrals are never mean stresses.
Connected variance of Q_trace. At leading Einstein order deltaR=-8piG deltaT.
"""
from __future__ import annotations
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[v]='1'
import argparse,json,time,hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import gammaln
from numpy.polynomial.legendre import leggauss
from threadpoolctl import threadpool_limits
from prepared_geometry import PreparedBackground,tortoise
from reference_fluctuations import radial_proper_primitive,cutoff,atomic_json,sha,clock
threadpool_limits(limits=1)
ROOT=Path(__file__).resolve().parents[1]
class Background:
 def __init__(self,rmin=.2):
  self.base=PreparedBackground(rmin);self.rmin=rmin
  t=np.linspace(0,2,12001);v=np.array([self.base.at(z) for z in t])
  r=np.geomspace(1.5,rmin,16001);ti=self.base.times(r);h=2/r-1
  vi=np.stack([r,h,-h,-2*h/r**2,2*h/r**2,np.zeros_like(r),np.zeros_like(r)],axis=-1)
  self.cs=CubicSpline(np.r_[t[:-1],ti],np.vstack([v[:-1],vi]),axis=0)
 def at(self,t):return self.cs(t)
 def rhs(self,t,y,k,L,seed):
  r,A,rp,rpp,Ap,m,mp=self.at(t);q=k*k+A*(L/r**2+seed*m);qp=Ap*(L/r**2+seed*m)+A*(-2*L*rp/r**3+seed*mp)
  w=np.sqrt(q);u,d,phase=y.reshape(3,-1)
  return np.array([2*w*d+d*d+rpp/r-u*u,-2*u*(w+d)-qp/(2*w),w+d]).ravel()

def window(bg,cfg):
 x,w=leggauss(cfg['nt']);tau=5*cfg['sigma']*x;P=radial_proper_primitive;rc=cfg['center']
 rlo=brentq(lambda r:P(r)-P(rc)+5*cfg['sigma'],bg.rmin,rc,xtol=7e-15)
 rhi=brentq(lambda r:P(r)-P(rc)-5*cfg['sigma'],rc,1.5,xtol=7e-15)
 rr=np.array([brentq(lambda r:P(r)-P(rc)+s,rlo,rhi,xtol=7e-15) for s in tau]);qw=w*5*np.exp(-(tau/cfg['sigma'])**2)*cutoff(tau/cfg['sigma']);qw/=qw.sum()
 return dict(tau=tau,r=rr,t=bg.base.times(rr),qw=qw,support=[rlo,rhi])

def wigner000_squared(l,p,J):
 if J<abs(l-p) or J>l+p or (l+p+J)%2:return 0.
 s=(l+p+J)//2
 logabs=(gammaln(s+1)-gammaln(s-l+1)-gammaln(s-p+1)-gammaln(s-J+1)
    +.5*(gammaln(2*s-2*l+1)+gammaln(2*s-2*p+1)+gammaln(2*s-2*J+1)-gammaln(2*s+2)))
 return float(np.exp(2*logabs))

def modes(bg,win,cfg,k,l):
 """Store full unequal-time modes in execution cache; final state and covariance in release."""
 out=ROOT/'results'/cfg['name'];cache=out/'mode_cache';cache.mkdir(exist_ok=True)
 dep=''.join(sha(ROOT/'code'/n) for n in ('localized_trace.py','prepared_geometry.py','reference_fluctuations.py'))
 fp=hashlib.sha256((json.dumps(cfg,sort_keys=True)+dep).encode()).hexdigest()
 p=cache/f'l{l:04d}.npz';meta=out/'mode_records'/f'l{l:04d}.json';meta.parent.mkdir(exist_ok=True)
 if p.exists() and meta.exists():
  rec=json.loads(meta.read_text())
  if rec['fingerprint']!=fp or rec['cache_hash']!=sha(p):raise RuntimeError('mode hash/config mismatch')
  a=np.load(p);return a['f'],a['d'],False
 t0=time.perf_counter();ll=np.full_like(k,l*(l+1))
 sol=solve_ivp(lambda t,y:bg.rhs(t,y,k,ll,cfg['seed']),[0,float(win['t'][-1])],np.zeros(3*len(k)),t_eval=win['t'],method='DOP853',rtol=cfg['rtol'],atol=cfg['rtol']/100,max_step=.02)
 if not sol.success:raise RuntimeError(sol.message)
 u,de,phase=sol.y.reshape(3,len(k),-1);r,A,rp,rpp,Ap,m,mp=bg.at(win['t']).T
 W=np.sqrt(k[:,None]**2+A[None]*(l*(l+1)/r[None]**2+cfg['seed']*m[None]))+de
 if np.min(W)<=0:raise FloatingPointError('nonpositive mode width')
 f=(np.exp(-1j*phase)/(r[None]*np.sqrt(2*W))).T;d=(u-1j*W-rp[None]/r[None]).T*f
 tmp=p.with_suffix('.tmp.npz');np.savez_compressed(tmp,f=f,d=d);tmp.replace(p)
 final=out/'final_states'/f'l{l:04d}.npz';final.parent.mkdir(exist_ok=True);np.savez_compressed(final,k=k,u=u[:,-1],delta=de[:,-1],phase=phase[:,-1],last_time=win['t'][-1])
 rec=dict(l=l,labels=len(k),nfev=sol.nfev,seconds=time.perf_counter()-t0,normalization=float(np.max(abs(r[:,None]**2*(f*d.conj()-d*f.conj())-1j))),fingerprint=fp,cache_hash=sha(p),final_hash=sha(final),utc=clock());atomic_json(meta,rec)
 return f,d,True

def moments_for_pair(f,d,ff,dd,r,A,qw,k,kw,b,L,P):
 """Integrate time first; all allowed external J obtained from B_J=B_0+J(J+1) R/2."""
 D=d.T@((qw/A)[:,None]*dd);K=f.T@((qw/A)[:,None]*ff);R=f.T@((qw/r**2)[:,None]*ff)
 M=D-.5*(L+P)*R;kk=k[:,None]*k[None,:]
 v00=0.;v01=0.;v11=0.
 for s in (1,-1):
  q=k[:,None]+s*k[None,:];weight=kw[:,None]*kw[None,:]*np.exp(-b*b*q*q/2)
  B=M+s*kk*K
  v00+=np.sum(weight*abs(B)**2);v01+=np.sum(weight*np.real(B*R.conj()));v11+=np.sum(weight*abs(R)**2)
 return np.array([v00,v01,v11])

def static_moments(l,p,k,kw,b,cfg):
 """Frozen product R x R x S2 reference with same proper kernels; not an Einstein solution."""
 r=cfg['center'];A=2/r-1;L=l*(l+1);P=p*(p+1)
 wi=np.sqrt(k*k+A*L/r**2);wj=np.sqrt(k*k+A*P/r**2)
 # Fourier transform of normalized compact Gaussian, independently vs exact Gaussian checked elsewhere.
 x,w=leggauss(max(128,cfg['nt']));z=5*x;tw=5*w*np.exp(-z*z)*cutoff(z);tw/=tw.sum()
 freq=(wi[:,None]+wj[None,:])/np.sqrt(A)
 # For compact cutoff outside 4 sigma, difference from exact Gaussian tiny; use direct transform here.
 transform=np.zeros_like(freq)
 for zz,ww in zip(z,tw):transform+=ww*np.cos(freq*cfg['sigma']*zz)
 product=transform/(2*r*r*np.sqrt(wi[:,None]*wj[None,:]))
 R=product/r**2;M=(-wi[:,None]*wj[None,:]/A-.5*(L+P)/r**2)*product;kk=k[:,None]*k[None,:]/A
 out=np.zeros(3)
 for s in (1,-1):
  q=k[:,None]+s*k[None,:];weight=kw[:,None]*kw[None,:]*np.exp(-b*b*q*q/2)
  B=M+s*kk*product
  out+= [np.sum(weight*abs(B)**2),np.sum(weight*B*R),np.sum(weight*R*R)]
 return out

def run(cfg,stop_after=None):
 out=ROOT/'results'/cfg['name'];out.mkdir(parents=True,exist_ok=True);bg=Background(cfg.get('rmin',.2));win=window(bg,cfg);atomic_json(out/'config.json',cfg)
 x,w=leggauss(cfg['nk']);k=(x+1)*cfg['K']/2;kw=w*cfg['K']/2;b=cfg['bproper']/np.sqrt(2/cfg['center']-1)
 r,A,*_=bg.at(win['t']).T
 deps=['localized_trace.py','prepared_geometry.py','reference_fluctuations.py'];fp=hashlib.sha256((json.dumps(cfg,sort_keys=True)+''.join(sha(ROOT/'code'/p) for p in deps)).encode()).hexdigest()
 total=np.zeros(cfg['Jmax']+1);static=np.zeros_like(total);records=[];completed=0;new_modes=0;skip=0
 # Each outer-l is independent and resumable. Bounded rolling LRU cache.
 from collections import OrderedDict
 cache=OrderedDict()
 def fetch(l):
  nonlocal new_modes
  if l not in cache:
   f,d,new=modes(bg,win,cfg,k,l);cache[l]=(f,d);new_modes+=int(new)*len(k)
   if len(cache)>cfg['Jmax']+2:cache.popitem(last=False)
  return cache[l]
 for l in range(cfg['L']+1):
  p=out/f'pairblock_{l:04d}.npz';meta=p.with_suffix('.json')
  if p.exists() and meta.exists():
   rec=json.loads(meta.read_text())
   if rec['fingerprint']!=fp or rec['sha256']!=sha(p):raise RuntimeError('pairblock hash/config mismatch')
   dat=np.load(p);total+=dat['spectrum'];static+=dat['static_spectrum'];records.append(rec);skip+=1;continue
  t0=time.perf_counter();f,d=fetch(l);spec=np.zeros_like(total);refs=np.zeros_like(total);mom=[];lm=[]
  for pp in range(l,min(cfg['L'],l+cfg['Jmax'])+1):
   ff,dd=fetch(pp);v=moments_for_pair(f,d,ff,dd,r,A,win['qw'],k,kw,b,l*(l+1),pp*(pp+1));vs=static_moments(l,pp,k,kw,b,cfg)
   mom.append(np.r_[v,vs]);lm.append(pp)
   for J in range(abs(l-pp),min(l+pp,cfg['Jmax'])+1,2):
    a=J*(J+1)/2;factor=(2 if l!=pp else 1)*(2*l+1)*(2*pp+1)*(2*J+1)*wigner000_squared(l,pp,J)/(16*np.pi**4)
    spec[J]+=factor*(v[0]+2*a*v[1]+a*a*v[2]);refs[J]+=factor*(vs[0]+2*a*vs[1]+a*a*vs[2])
  tmp=p.with_suffix('.tmp.npz');np.savez_compressed(tmp,spectrum=spec,static_spectrum=refs,l=l,partners=np.array(lm),moments=np.array(mom));tmp.replace(p)
  rec=dict(l=l,pairs=len(lm),seconds=time.perf_counter()-t0,fingerprint=fp,sha256=sha(p),utc=clock());atomic_json(meta,rec);records.append(rec);total+=spec;static+=refs;completed+=1
  atomic_json(out/'STATUS.json',dict(complete=False,records=records,new_labels=new_modes,skipped_blocks=skip,support=win['support'],fingerprint=fp))
  np.savez(out/'PARTIAL_SPECTRA.npz',spectrum=total,static_spectrum=static,J=np.arange(len(total)))
  print('SAVED',cfg['name'],l,round(rec['seconds'],3),flush=True)
  if stop_after and completed>=stop_after:print('PLANNED STOP',flush=True);return
 atomic_json(out/'STATUS.json',dict(complete=True,records=records,new_labels=new_modes,skipped_blocks=skip,support=win['support'],fingerprint=fp))
 np.savez(out/'SPECTRA.npz',spectrum=total,static_spectrum=static,J=np.arange(len(total)))
 print('COMPLETE',cfg['name'],total[0],static[0],flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('config');ap.add_argument('--stop-after',type=int);args=ap.parse_args();run(json.loads(Path(args.config).read_text()),args.stop_after)
