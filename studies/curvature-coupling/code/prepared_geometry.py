"""Equation-derived first metric perturbation, with documented auxiliary past.

No Jiang input. Numerical tensor is an archived SPV-1 source, NOT recomputed here.
Geometric mass M=1. m=1+e*mu1+...; psi=e*psi1+... .
"""
from __future__ import annotations
import numpy as np
from numpy.polynomial import Chebyshev
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit
from pathlib import Path
import json, hashlib
ROOT=Path(__file__).resolve().parents[1]

def smooth(t,T=2.):
    z=t/T
    if z<=0: return 0.,0.,0.
    if z>=1: return 1.,0.,0.
    f=-1/z+1/(1-z); fp=1/z**2+1/(1-z)**2
    fpp=-2/z**3+2/(1-z)**3
    s=float(expit(f)); v=s*(1-s)
    return s,v*fp/T,v*((1-2*s)*fp*fp+fpp)/(T*T)

def tortoise(r): return r+2*np.log1p(-np.asarray(r)/2)

class PreparedBackground:
    """The prior SPV-1 auxiliary preparation, reproduced explicitly."""
    def __init__(self,rmin=.4):
        self.T=2.; self.rmatch=1.5
        self.rprep=brentq(lambda r:tortoise(r)-(tortoise(1.5)-1),1.5,1.999999,xtol=1e-14)
        self.end=float(2+tortoise(rmin)-tortoise(1.5))
        self.sol=solve_ivp(lambda t,y:[-(2/y[0]-1)*smooth(t)[0]],(0,self.end),[self.rprep],rtol=2.3e-14,atol=1e-15,method='DOP853',max_step=.015,dense_output=True)
        if not self.sol.success: raise RuntimeError(self.sol.message)
    def at(self,t):
        r=float(self.sol.sol(t)[0]); h=2/r-1; s,sd,_=smooth(t)
        rp=-h*s; rpp=-2*h*s*s/r**2-h*sd
        A=1-s+s*h; Ap=sd*(h-1)+2*h*s*s/r**2
        return r,A,rp,rpp,Ap,1-s,-sd
    def times(self,r): return 2+tortoise(r)-tortoise(1.5)

class FirstGeometry:
    def __init__(self,degree=32,fit_top=1.4):
        self.degree=degree
        p=ROOT/'inputs/spv1_cached_vacuum.csv'
        self.input_sha=hashlib.sha256(p.read_bytes()).hexdigest()
        d=np.genfromtxt(p,names=True,delimiter=',')
        mask=(d['r_over_M']<=fit_top+1e-12)&(d['r_over_M']>=.4-1e-12)
        r=d['r_over_M'][mask]
        self.domain=(.4,fit_top)
        self.fits=[]
        for col in ('rho_M4_over_hbar','p_long_M4_over_hbar','p_ang_M4_over_hbar'):
            self.fits.append(Chebyshev.fit(r,r**6*d[col][mask],degree,domain=self.domain))
        def rhs(r,y):
            rho,px,_=self.stress(r)
            return [-4*np.pi*r*r*px,4*np.pi*r/(2/r-1)*(rho+px)]
        kw=dict(method='DOP853',rtol=2e-12,atol=2e-14,max_step=.004,dense_output=True)
        self.inner=solve_ivp(rhs,(1,.4),[0.,0.],**kw)
        self.outer=solve_ivp(rhs,(1,fit_top),[0.,0.],**kw)
        if not (self.inner.success and self.outer.success): raise RuntimeError('geometry integration failed')
    def stress(self,r,der=0):
        r=np.asarray(r)
        if der==0:return np.asarray([f(r)/r**6 for f in self.fits])
        if der==1:return np.asarray([f.deriv()(r)/r**6-6*f(r)/r**7 for f in self.fits])
        if der==2:return np.asarray([f.deriv(2)(r)/r**6-12*f.deriv()(r)/r**7+42*f(r)/r**8 for f in self.fits])
        raise ValueError('unsupported derivative')
    def values(self,r):
        r=float(r)
        mu,psi=(self.inner if r<=1 else self.outer).sol(r)
        rho,px,_=self.stress(r); rhop,pxp,_=self.stress(r,1)
        h=2/r-1; hp=-2/r**2
        mup=-4*np.pi*r*r*px; mupp=-4*np.pi*(2*r*px+r*r*pxp)
        D=4*np.pi*r/h; Dp=4*np.pi*(1/h-r*hp/h**2)
        psip=D*(rho+px); psipp=Dp*(rho+px)+D*(rhop+pxp)
        q=2-r
        alpha=2*mu/q; alphap=2*mup/q+2*mu/q**2
        alphapp=2*mupp/q+4*mup/q**2+4*mu/q**3
        return np.array([mu,psi,mup,psip,mupp,psipp,alpha,alphap,alphapp])
    def deformations(self,r):
        """a1=delta log(A), c1=delta log(C), and radial derivatives.
        ds2=A(-C2 deta2+dx2)+r2dOmega2 with fixed background r(eta).
        Auxiliary switch is not an Einstein solution outside r<=1.
        """
        r=float(r)
        if r>=1.4:return np.zeros(6)
        mu,psi,mp,pp,mpp,ppp,alpha,alp,alpp=self.values(r)
        ar=np.array([2*psi+alpha,2*pp+alp,2*ppp+alpp])
        cr=np.array([-psi-alpha,-pp-alp,-ppp-alpp])
        S,St,Stt=smooth(1.4-r,.3); Sr=-St; Srr=Stt
        def prod(v):return np.array([S*v[0],Sr*v[0]+S*v[1],Srr*v[0]+2*Sr*v[1]+S*v[2]])
        return np.r_[prod(ar),prod(cr)]

def run_geometry(degree=32):
    import time
    t=time.perf_counter(); fg=FirstGeometry(degree)
    r=np.linspace(1.4,.4,501)
    stress=np.array([fg.stress(x) for x in r]); vals=np.array([fg.values(x) for x in r])
    defs=np.array([fg.deformations(x) for x in r])
    F=1-2/r; Fp=2/r**2
    mu,psi,mp,pp,mpp,ppp=vals[:,:6].T
    f1=-2*mu/r; f1p=-2*mp/r+2*mu/r**2;f1pp=-2*mpp/r+4*mp/r**2-4*mu/r**3
    Gtheta=F*ppp+(1.5*Fp+F/r)*pp+.5*f1pp+f1p/r
    scale=8*np.pi*(np.abs(stress[:,0])+np.abs(stress[:,1])+2*np.abs(stress[:,2]))
    res=(Gtheta-8*np.pi*stress[:,2])/scale
    rp=np.array([fg.stress(x,1)[0] for x in r]); rho,px,pa=stress.T
    cons=rp+(-1/(r*(2-r)))*(rho+px)+2/r*(rho+pa)
    cscale=(np.abs(rp)+np.abs((rho+px)/(r*(2-r)))+np.abs(2/r*(rho+pa)))
    csv=np.c_[r,stress,vals,defs,res,cons/cscale]
    header='r,rho0,px0,pa0,mu1,psi1,mu1_prime,psi1_prime,mu1_second,psi1_second,alpha1,alpha1_prime,alpha1_second,a1,a1_prime,a1_second,c1,c1_prime,c1_second,angular_residual_normalized,conservation_normalized'
    out=ROOT/'results'/f'geometry_degree{degree}.csv'
    np.savetxt(out,csv,delimiter=',',header=header,comments='')
    inside=(r>=.45)&(r<=.95)
    meta={'degree':degree,'input_sha256':fg.input_sha,'source':'cached SPV-1; not newly evolved vacuum','matching_r':1.0,'auxiliary_switch':[1.4,1.1], 'Einstein_reporting_domain':[.4,1.0], 'check_domain':[.45,.95],
          'max_angular_residual':float(np.max(np.abs(res[inside]))),'max_source_conservation':float(np.max(np.abs(cons[inside]/cscale[inside]))),'max_first_metric_coefficient':float(np.max(np.maximum(np.abs(vals[r<=1,6]),2*np.abs(vals[r<=1,1])))), 'seconds':time.perf_counter()-t,'r04_values':fg.values(.4).tolist()}
    (ROOT/'results'/f'geometry_degree{degree}.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2))
if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--degree',type=int,default=32);args=a.parse_args();run_geometry(args.degree)
