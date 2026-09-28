"""Fixed-state improved stress covariance from verified same-window histories.
No new field modes evolved by this program. Independent Huu and IBP checks elsewhere.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,hashlib,time
from collections import OrderedDict
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.chebyshev import chebvander
from scipy.special import expit
from threadpoolctl import threadpool_limits
import localized_trace as core
threadpool_limits(limits=1)
ROOT=Path(__file__).resolve().parents[1]

def second_kernel_ratio(z,sigma):
    """g''/g for compact exp(-z^2)*expit(1/u-1/(1-u)); analytic derivatives."""
    z=np.asarray(z);ratio=4*z*z-2;sel=(abs(z)>4)&(abs(z)<5)
    u=abs(z[sel])-4;v=1/u-1/(1-u);C=expit(v)
    vp=-1/u**2-1/(1-u)**2;vpp=2/u**3-2/(1-u)**3
    first=(1-C)*vp;second=(1-C)*vpp-C*(1-C)*vp*vp
    # Even kernel, use positive z for both symmetric tails.
    zz=abs(z[sel]);ratio[sel]=(-2*zz+first)**2-2+second
    return ratio/sigma**2

def load_hist(source,name,l,x,r,rp):
    p=source/'results'/f'{name}_histories'/f'l{l:04d}.npz'
    z=np.load(p);nk=len(z['k']);val=chebvander(x,int(z['degree']))@z['coefficients']
    W,u,ph=np.split(val,3,axis=1)
    if W.min()<=0:raise ValueError('nonpositive interpolated frequency')
    ph+=z['phase_offset'][None,:]
    f=np.exp(-1j*ph)/(r[:,None]*np.sqrt(2*W));d=(u-1j*W-rp[:,None]/r[:,None])*f
    return f,d,z['k']

def gram_pair(D,K,R,H,k,weights,lp,cutmask=None):
    B=D-lp*R/2;kk=k[:,None]*k[None,:];g=np.zeros((4,4))
    for sign,w in zip((1,-1),weights):
        b=B+sign*kk*K;e=D-(H+b)/6
        vec=[b,R,D,e]
        ww=w if cutmask is None else w*cutmask
        for a in range(4):
            for c in range(a,4):
                v=float(np.sum(ww*np.real(vec[a]*vec[c].conj())))
                g[a,c]+=v
                if c!=a:g[c,a]+=v
    return g

def to_basis(g,J):
    a=J*(J+1)/2
    # [T_minimal, P_minimal, S_uu_conformal]
    t=np.array([[1,a,0,0],[0,0,1,0],[0,-a/6,0,1.]])
    return t@g@t.T

def run(config,source,stop=None):
    name=config['name'];prior=config['prior'];out=ROOT/'results'/name;out.mkdir(parents=True,exist_ok=True)
    old=json.loads((source/'inputs'/f'{prior}.json').read_text());cfg={**old,'nt':config.get('nt',old['nt'])}
    L=config.get('L',old['L']);Jmax=config.get('Jmax',old['Jmax'])
    bg=core.Background(old['rmin']);win=core.window(bg,cfg);r,A,rp,rpp,Ap,*_=bg.at(win['t']).T
    x=win['tau']/(5*cfg['sigma']);nodes,ww=leggauss(old['nk']);k=(nodes+1)*old['K']/2;kw=ww*old['K']/2
    b=cfg['bproper']/np.sqrt(2/cfg['center']-1)
    weights=[kw[:,None]*kw[None,:]*np.exp(-b*b*(k[:,None]+s*k[None,:])**2/2) for s in(1,-1)]
    cut=(k[:,None]<=.8*old['K'])&(k[None,:]<=.8*old['K'])
    kern2=second_kernel_ratio(win['tau']/cfg['sigma'],cfg['sigma'])
    qw=win['qw'];q2=qw*kern2
    srcmanifest=json.loads((source/'FILE_HASHES.json').read_text())
    histmanifest=json.loads((source/'results'/f'{prior}_histories'/'MANIFEST.json').read_text())
    hh={rec['l']:rec['sha256'] for rec in histmanifest['records']}
    fp=hashlib.sha256((json.dumps(config,sort_keys=True)+core.sha(Path(__file__))+
       core.sha(source/'FILE_HASHES.json')).encode()).hexdigest()
    cache=OrderedDict();total=np.zeros((Jmax+1,3,3));truncated=np.zeros_like(total);ktail=np.zeros_like(total);records=[];skipped=0;made=0
    def fetch(l):
        if l not in cache:
            hp=source/'results'/f'{prior}_histories'/f'l{l:04d}.npz'
            if core.sha(hp)!=hh[l]:raise ValueError(f'history hash failed l{l}')
            f,d,ks=load_hist(source,prior,l,x,r,rp)
            if not np.allclose(ks,k,rtol=0,atol=0):raise ValueError('momentum grid mismatch')
            cache[l]=(f,d)
            if len(cache)>Jmax+2:cache.popitem(last=False)
        return cache[l]
    def save(done):
        core.atomic_json(out/'STATUS.json',dict(config=config,prior_config=old,complete=done,records=records,
            fingerprint=fp,skipped_blocks=skipped,new_pair_blocks=made,new_field_evolutions=0,
            support=win['support'],mean_sqrtK=float(qw@(np.sqrt(48)/r**3)),history_max_reported_error=histmanifest['max_error']))
        np.savez_compressed(out/('SPECTRA.npz' if done else 'PARTIAL.npz'),covariance=total,
            truncated_angular=truncated,truncated_momentum=ktail,J=np.arange(Jmax+1))
    for l in range(L+1):
        p=out/f'block_{l:04d}.npz';m=p.with_suffix('.json')
        if p.exists() and m.exists():
            rec=json.loads(m.read_text())
            if rec['fingerprint']!=fp or core.sha(p)!=rec['sha256']:raise ValueError('checkpoint hash mismatch')
            z=np.load(p);total+=z['cov'];truncated+=z['lcov'];ktail+=z['kcov'];records.append(rec);skipped+=1;continue
        t0=time.perf_counter();f,d=fetch(l);spec=np.zeros_like(total);lspec=np.zeros_like(total);kspec=np.zeros_like(total);mom=[];partners=[]
        for pp in range(l,min(L,l+Jmax)+1):
            ff,dd=fetch(pp)
            D=d.T@((qw/A)[:,None]*dd)
            K=f.T@((qw/A)[:,None]*ff)
            R=f.T@((qw/r**2)[:,None]*ff)
            H=f.T@(q2[:,None]*ff)
            g=gram_pair(D,K,R,H,k,weights,l*(l+1)+pp*(pp+1))
            gc=gram_pair(D,K,R,H,k,weights,l*(l+1)+pp*(pp+1),cut)
            partners.append(pp);mom.append(g)
            for J in range(pp-l,min(l+pp,Jmax)+1,2):
                fac=(2 if l!=pp else 1)*(2*l+1)*(2*pp+1)*(2*J+1)*core.wigner000_squared(l,pp,J)/(16*np.pi**4)
                spec[J]+=fac*to_basis(g,J);kspec[J]+=fac*to_basis(gc,J)
                if pp<=L-32:lspec[J]+=fac*to_basis(g,J)
        tmp=p.with_suffix('.tmp.npz');np.savez_compressed(tmp,cov=spec,lcov=lspec,kcov=kspec,l=l,partners=partners,moments=mom);tmp.replace(p)
        rec=dict(l=l,sha256=core.sha(p),fingerprint=fp,seconds=time.perf_counter()-t0,utc=core.clock(),pair_count=len(partners))
        core.atomic_json(m,rec);records.append(rec);total+=spec;truncated+=lspec;ktail+=kspec;made+=1;save(False)
        print('SAVED',name,l,round(rec['seconds'],3),flush=True)
        if stop and made>=stop:return
    save(True);print('COMPLETE',name,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('config');p.add_argument('--source',type=Path,required=True);p.add_argument('--stop-after',type=int);a=p.parse_args()
    run(json.loads(Path(a.config).read_text()),a.source,a.stop_after)
