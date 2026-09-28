"""Equivalent compact-kernel Fourier reference, accelerated by integrating only tiny removed tails.
Physical modes and contractions use the unchanged localized_trace.py engine.
"""
import json,hashlib,argparse
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
import localized_trace as core

def compact_fourier(v,n=40):
 x,w=leggauss(n);z=np.r_[4+(x+1)/2,5+2*(x+1)];wt=np.r_[w/2,2*w]
 rem=2/np.sqrt(np.pi)*wt*np.exp(-z*z)*(1-core.cutoff(z))
 v=np.asarray(v);out=np.exp(-v*v/4)
 for zz,ww in zip(z,rem):out-=ww*np.cos(v*zz)
 return out/(1-rem.sum())

def static_fast(l,p,k,kw,b,cfg):
 r=cfg['center'];A=2/r-1;L=l*(l+1);P=p*(p+1)
 wi=np.sqrt(k*k+A*L/r**2);wj=np.sqrt(k*k+A*P/r**2)
 ft=compact_fourier((wi[:,None]+wj[None,:])*cfg['sigma']/np.sqrt(A))
 prod=ft/(2*r*r*np.sqrt(wi[:,None]*wj[None,:]));R=prod/r**2
 M=(-wi[:,None]*wj[None,:]/A-.5*(L+P)/r**2)*prod;kk=k[:,None]*k[None,:]/A
 out=np.zeros(3)
 for s in (1,-1):
  q=k[:,None]+s*k[None,:];weight=kw[:,None]*kw[None,:]*np.exp(-b*b*q*q/2);B=M+s*kk*prod
  out+=[np.sum(weight*abs(B)**2),np.sum(weight*B*R),np.sum(weight*R*R)]
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('config');a=ap.parse_args();cfg=json.loads(Path(a.config).read_text())
 if cfg.get('reference_engine_hash')!=core.sha(__file__):raise RuntimeError('fast wrapper hash required')
 core.static_moments=static_fast;core.run(cfg)
