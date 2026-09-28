"""New joint angular trace/pi^2 covariance on matched inward windows.
Physical oscillator code is reused with provenance; each campaign evolves new modes.
All m sums are exact through Gaunt weights, all retained field k,l are explicit.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import sys,json,hashlib,time,argparse
from pathlib import Path
from collections import OrderedDict
import numpy as np
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import leggauss
from threadpoolctl import threadpool_limits
import localized_trace as core
from localized_trace_fast import compact_fourier
threadpool_limits(limits=1)
ROOT=Path(__file__).resolve().parents[1]

def joint_moments(D,K,R,k,weights,lp):
    kk=k[:,None]*k[None,:]; M=D-lp*R/2
    v=np.zeros(6)
    for s,w in zip((1,-1),weights):
        B=M+s*kk*K
        v+=np.array([np.sum(w*abs(B)**2),np.sum(w*np.real(B*R.conj())),
            np.sum(w*abs(R)**2),np.sum(w*abs(D)**2),
            np.sum(w*np.real(B*D.conj())),np.sum(w*np.real(R*D.conj()))])
    return v

def physical_moments(f,d,ff,dd,r,A,qw,k,weights,lp):
    D=d.T@((qw/A)[:,None]*dd)
    K=f.T@((qw/A)[:,None]*ff)
    R=f.T@((qw/r**2)[:,None]*ff)
    return joint_moments(D,K,R,k,weights,lp)

class FrozenReference:
    def __init__(self,cfg):
        self.cfg=cfg
        r=cfg['center']; A=2/r-1
        vmax=2*np.sqrt(cfg['K']**2+A*cfg['L']*(cfg['L']+1)/r**2)*cfg['sigma']/np.sqrt(A)
        xs=np.linspace(0,max(1,vmax)*1.001,16385)
        self.cs=CubicSpline(xs,compact_fourier(xs),extrapolate=False)
        rng=np.random.default_rng(2060926);test=rng.uniform(0,vmax,257)
        self.absolute_error=float(np.max(abs(self.cs(test)-compact_fourier(test))))
        if self.absolute_error>1e-10:raise RuntimeError('Fourier interpolation tolerance')
    def moments(self,l,p,k,weights):
        c=self.cfg;r=c['center'];A=2/r-1;L=l*(l+1);P=p*(p+1)
        wi=np.sqrt(k*k+A*L/r**2);wj=np.sqrt(k*k+A*P/r**2)
        prod=self.cs((wi[:,None]+wj[None,:])*c['sigma']/np.sqrt(A))/(2*r*r*np.sqrt(wi[:,None]*wj[None,:]))
        D=-wi[:,None]*wj[None,:]*prod/A
        return joint_moments(D,prod/A,prod/r**2,k,weights,L+P)

def evaluate(v,J):
    a=J*(J+1)/2
    return np.array([[v[0]+2*a*v[1]+a*a*v[2],v[4]+a*v[5]],
                     [v[4]+a*v[5],v[3]]])

def run(cfg,stop_after=None):
    out=ROOT/'results'/cfg['name'];out.mkdir(parents=True,exist_ok=True)
    core.ROOT=ROOT
    dep=['joint_inward.py','localized_trace.py','localized_trace_fast.py','prepared_geometry.py','reference_fluctuations.py']
    fp=hashlib.sha256((json.dumps(cfg,sort_keys=True)+''.join(core.sha(ROOT/'code'/p) for p in dep)).encode()).hexdigest()
    core.atomic_json(out/'config.json',cfg)
    bg=core.Background(cfg['rmin']);win=core.window(bg,cfg)
    x,w=leggauss(cfg['nk']);k=(x+1)*cfg['K']/2;kw=w*cfg['K']/2
    b=cfg['bproper']/np.sqrt(2/cfg['center']-1)
    r,A,*_=bg.at(win['t']).T
    weights=[kw[:,None]*kw[None,:]*np.exp(-b*b*(k[:,None]+s*k[None,:])**2/2) for s in (1,-1)]
    static=FrozenReference(cfg)
    total=np.zeros((cfg['Jmax']+1,2,2));ref=np.zeros_like(total);records=[];newlabels=0;skipped=0;made=0
    cache=OrderedDict()
    def fetch(l):
        nonlocal newlabels
        if l not in cache:
            f,d,new=core.modes(bg,win,cfg,k,l);cache[l]=(f,d);newlabels+=int(new)*len(k)
            if len(cache)>cfg['Jmax']+2:cache.popitem(last=False)
        return cache[l]
    def save(complete):
        core.atomic_json(out/'STATUS.json',dict(complete=complete,config=cfg,records=records,new_mode_labels_this_invocation=newlabels,
            skipped_blocks=skipped,fingerprint=fp,support=win['support'],mean_sqrtK=float(win['qw']@(np.sqrt(48)/r**3)),
            fourier_interpolation_absolute_error=static.absolute_error))
        np.savez_compressed(out/('SPECTRA.npz' if complete else 'PARTIAL.npz'),covariance=total,reference=ref,J=np.arange(len(total)))
    for l in range(cfg['L']+1):
        p=out/f'joint_{l:04d}.npz';meta=p.with_suffix('.json')
        if p.exists() and meta.exists():
            rec=json.loads(meta.read_text())
            if rec['fingerprint']!=fp or core.sha(p)!=rec['sha256']:raise RuntimeError('joint block mismatch')
            z=np.load(p);total+=z['covariance'];ref+=z['reference'];records.append(rec);skipped+=1;continue
        t0=time.perf_counter();f,d=fetch(l);spec=np.zeros_like(total);refs=np.zeros_like(total);moms=[];partners=[]
        for pp in range(l,min(cfg['L'],l+cfg['Jmax'])+1):
            ff,dd=fetch(pp);v=physical_moments(f,d,ff,dd,r,A,win['qw'],k,weights,l*(l+1)+pp*(pp+1))
            vs=static.moments(l,pp,k,weights);partners.append(pp);moms.append(np.r_[v,vs])
            for J in range(pp-l,min(l+pp,cfg['Jmax'])+1,2):
                fac=(2 if l!=pp else 1)*(2*l+1)*(2*pp+1)*(2*J+1)*core.wigner000_squared(l,pp,J)/(16*np.pi**4)
                spec[J]+=fac*evaluate(v,J);refs[J]+=fac*evaluate(vs,J)
        tmp=p.with_suffix('.tmp.npz');np.savez_compressed(tmp,covariance=spec,reference=refs,l=l,partners=partners,moments=moms);tmp.replace(p)
        rec=dict(l=l,seconds=time.perf_counter()-t0,sha256=core.sha(p),fingerprint=fp,utc=core.clock(),pairs=len(partners))
        core.atomic_json(meta,rec);records.append(rec);total+=spec;ref+=refs;made+=1;save(False)
        print('SAVED',cfg['name'],l,round(rec['seconds'],3),flush=True)
        if stop_after and made>=stop_after:print('PLANNED STOP',flush=True);return
    save(True);print('COMPLETE',cfg['name'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');p.add_argument('--stop-after',type=int);a=p.parse_args()
    run(json.loads(Path(a.config).read_text()),a.stop_after)
