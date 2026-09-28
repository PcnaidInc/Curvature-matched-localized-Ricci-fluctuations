"""New spacetime-smeared stress covariance from positive vacuum pair amplitudes.
Prior source supplies only SPV-1 background/state and a cached mean for comparison.
No original mean subtraction, Jiang data, or alleged lost covariance is loaded.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import json, time, hashlib, argparse, datetime
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import expit
from threadpoolctl import threadpool_limits
from prepared_geometry import PreparedBackground
threadpool_limits(limits=1)
ROOT=Path(__file__).resolve().parents[1]
COMPONENTS=['rho','p_long','p_ang','T_hat0hatx']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clock():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def atomic_json(p,v):
    p=Path(p);q=p.with_suffix('.tmp.json');q.write_text(json.dumps(v,indent=2));q.replace(p)
def radial_proper_primitive(r):
    r=np.asarray(r);return 2*np.arcsin(np.sqrt(r/2))-np.sqrt(r*(2-r))
def cutoff(z):
    z=np.abs(np.asarray(z)); out=np.ones_like(z)
    out[z>=5]=0; sel=(z>4)&(z<5);u=z[sel]-4
    out[sel]=expit(1/u-1/(1-u));return out
class Background:
    def __init__(self):
        self.base=PreparedBackground(.4)
        tt=np.linspace(0,self.base.end,16001)
        vals=np.array([self.base.at(t) for t in tt]); self.cs=CubicSpline(tt,vals,axis=0)
    def at(self,t):return self.cs(t)
    def rhs(self,t,y,k,L):
        r,A,rp,rpp,Ap,m2,m2p=self.at(t)
        q=k*k+A*(L/r**2+m2)
        qp=Ap*(L/r**2+m2)+A*(-2*L*rp/r**3+m2p)
        w=np.sqrt(q);wp=qp/(2*w);u,d,phase=y.reshape(3,-1)
        return np.array([2*w*d+d*d+rpp/r-u*u,-2*u*(w+d)-wp,w+d]).ravel()
def windows(bg,nt,centers,sigmas):
    x,weights=leggauss(nt); ans=[]
    for rc in centers:
        for sig in sigmas:
            z=5*x; dt=sig*z
            if radial_proper_primitive(rc)-dt.min()>radial_proper_primitive(1.5) or radial_proper_primitive(rc)-dt.max()<radial_proper_primitive(.4):
                raise ValueError('measurement support exceeds common final-state band')
            rr=np.array([brentq(lambda r:radial_proper_primitive(r)-(radial_proper_primitive(rc)-v),.4,1.5,xtol=8e-15) for v in dt])
            ts=bg.base.times(rr)
            qw=weights*5*np.exp(-z*z)*cutoff(z);qw/=qw.sum()
            ans.append(dict(center=rc,sigma=sig,r=rr,t=ts,weight=qw))
    alltimes=np.unique(np.concatenate([v['t'] for v in ans]))
    for win in ans:win['indices']=np.searchsorted(alltimes,win['t'])
    return ans,alltimes

def amplitudes(f,d,rr,A,qw,k,L):
    """Time integrals of bilinear pair terms, never square an averaged field.
    Arrays f,d have (time,k); d=partial_eta f. No conjugation in aa amplitudes.
    """
    weight=qw/A
    Jd=d.T@(weight[:,None]*d)
    Jk=f.T@(weight[:,None]*f)
    Ja=f.T@((qw*L/rr**2)[:,None]*f)
    Jdf=d.T@(weight[:,None]*f)
    ans=[]
    for sign in (1,-1):
        kk=sign*k[:,None]*k[None,:]
        R=Jd-kk*Jk+Ja
        P=Jd-kk*Jk-Ja
        Pa=Jd+kk*Jk
        Flux=1j*(sign*k[None,:]*Jdf + k[:,None]*Jdf.T)
        ans.append(np.array([R,P,Pa,Flux]))
    return ans

def covariance_of_pairs(pair,k,kw,b):
    cov=np.zeros((4,4))
    for sign,B in zip((1,-1),pair):
        window=np.exp(-b*b*(k[:,None]+sign*k[None,:])**2/2)
        measure=kw[:,None]*kw[None,:]*window
        V=(B*np.sqrt(measure)[None,:,:]).reshape(4,-1)
        cov+=np.real(V@V.conj().T)
    # Simultaneously reversing both momenta flips the flux and nothing else.
    # This makes flux-scalar symmetrized covariances exactly zero for an even window.
    cov[:3,3]=cov[3,:3]=0
    return cov/(64*np.pi**4)

def run(cfg,stop_after=None):
    started=time.perf_counter();bg=Background()
    wins,tt=windows(bg,cfg['nt'],cfg['centers'],cfg['sigmas'])
    x,w=leggauss(cfg['nk']);k=(x+1)*cfg['K']/2;kw=w*cfg['K']/2
    out=ROOT/'results'/cfg['name'];out.mkdir(parents=True,exist_ok=True)
    fingerprint=hashlib.sha256((json.dumps(cfg,sort_keys=True)+sha(__file__)+sha(ROOT/'code/prepared_geometry.py')).encode()).hexdigest()
    jobs=[];new=0;skipped=0
    for lo in range(0,cfg['L']+1,cfg['block']):
        hi=min(lo+cfg['block'],cfg['L']+1);p=out/f'block_{lo:04d}_{hi-1:04d}.npz';meta=p.with_suffix('.json')
        if p.exists() and meta.exists():
            m=json.loads(meta.read_text())
            if m['fingerprint']!=fingerprint or sha(p)!=m['sha256']:raise RuntimeError('CHECKPOINT HASH/CONFIG MISMATCH '+str(p))
            skipped+=1;jobs.append(m);print('SKIP',cfg['name'],lo,flush=True);continue
        begin=time.perf_counter();ell=np.arange(lo,hi)
        kk=np.broadcast_to(k,(hi-lo,len(k))).copy().ravel();LL=np.broadcast_to((ell*(ell+1))[:,None],(hi-lo,len(k))).copy().ravel()
        sol=solve_ivp(lambda t,y:bg.rhs(t,y,kk,LL),(0,float(tt[-1])),np.zeros(3*len(kk)),t_eval=tt,
            rtol=cfg['rtol'],atol=cfg['rtol']/100,method='DOP853',max_step=.02)
        if not sol.success:raise RuntimeError(sol.message)
        state=sol.y.reshape(3,hi-lo,len(k),len(tt))
        cov=np.zeros((hi-lo,len(wins),len(cfg['bproper']),4,4))
        norms=[]
        for ii,l in enumerate(ell):
            for iw,win in enumerate(wins):
                r,A,rp,rpp,Ap,m,mp=bg.at(win['t']).T
                u,delta,phase=state[:,ii,:,:][:,:,win['indices']]
                w0=np.sqrt(k[:,None]**2+A[None,:]*(l*(l+1)/r[None,:]**2+m[None,:]))
                W=w0+delta
                if np.min(W)<=0:raise FloatingPointError('nonpositive mode width')
                f=(np.exp(-1j*phase)/(np.sqrt(2*W)*r[None,:])).T
                d=((u-1j*W-rp[None,:]/r[None,:]).T)*f
                pair=amplitudes(f,d,r,A,win['weight'],k,l*(l+1))
                for ib,bp in enumerate(cfg['bproper']):
                    b=bp/np.sqrt(2/win['center']-1)
                    cov[ii,iw,ib]=(2*l+1)*covariance_of_pairs(pair,k,kw,b)
                norms.append(float(np.max(np.abs(r[:,None]**2*(f*d.conj()-d*f.conj())-1j))))
        tmp=p.with_suffix('.tmp.npz')
        np.savez_compressed(tmp,covariance=cov,ell=ell,k=k,kw=kw,
            centers=np.array([v['center'] for v in wins]),sigmas=np.array([v['sigma'] for v in wins]),
            bproper=np.array(cfg['bproper']),last_time=tt[-1],final_state=state[:,:,:,-1])
        tmp.replace(p)
        m=dict(fingerprint=fingerprint,sha256=sha(p),lo=lo,hi=hi,labels=len(kk),seconds=time.perf_counter()-begin,
               nfev=sol.nfev,normalization_error=max(norms),utc=clock())
        atomic_json(meta,m);jobs.append(m);new+=1
        atomic_json(out/'STATUS.json',dict(config=cfg,records=jobs,new_blocks=new,skipped=skipped,complete=False))
        print('SAVED',cfg['name'],lo,hi,round(m['seconds'],2),'s',flush=True)
        if stop_after and new>=stop_after:print('PLANNED STOP',flush=True);return
    atomic_json(out/'STATUS.json',dict(config=cfg,records=jobs,new_blocks=new,skipped=skipped,complete=True,seconds=time.perf_counter()-started))
    print('COMPLETE',cfg['name'],round(time.perf_counter()-started,2),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('config');a.add_argument('--stop-after',type=int);v=a.parse_args();run(json.loads(Path(v.config).read_text()),v.stop_after)
