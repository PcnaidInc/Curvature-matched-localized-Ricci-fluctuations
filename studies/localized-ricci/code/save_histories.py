"""Supplementary phase-resolved interpolation caches, not new field evolutions.
Production pair moments stay unchanged. Each fit is independently reconstruction-tested.
"""
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[n]='1'
from pathlib import Path
import argparse,json
import numpy as np
from numpy.polynomial.chebyshev import chebvander
import localized_trace as c
ROOT=Path(__file__).resolve().parents[1]
def run(name):
 cfg=json.loads((ROOT/'inputs'/f'{name}.json').read_text());bg=c.Background(cfg['rmin']);win=c.window(bg,cfg)
 x=win['tau']/(5*cfg['sigma']);r,A,rp,*_=bg.at(win['t']).T;out=ROOT/'results'/f'{name}_histories';out.mkdir(exist_ok=True)
 mats={d:chebvander(x,d) for d in (32,48,64,96)};pseudo={d:np.linalg.pinv(v) for d,v in mats.items()};records=[]
 for l in range(cfg['L']+1):
  src=ROOT/'results'/name/'mode_cache'/f'l{l:04d}.npz';meta=json.loads((ROOT/'results'/name/'mode_records'/f'l{l:04d}.json').read_text());assert c.sha(src)==meta['cache_hash']
  z=np.load(src);f=z['f'];d=z['d'];W=-np.imag(d/f);u=np.real(d/f)+rp[:,None]/r[:,None];ph=-np.unwrap(np.angle(f),axis=0);off=ph[len(x)//2].copy();ph-=off;targets=np.concatenate([W,u,ph],axis=1)
  for deg in mats:
   co=pseudo[deg]@targets;wf,uf,pf=np.split(mats[deg]@co,3,axis=1)
   if wf.min()<=0:continue
   ff=np.exp(-1j*(pf+off[None]))/(r[:,None]*np.sqrt(2*wf));dd=(uf-1j*wf-rp[:,None]/r[:,None])*ff
   err=max(float(np.max(abs(ff-f))/np.max(abs(f))),float(np.max(abs(dd-d))/np.max(abs(d))))
   if err<5e-10:break
  if err>=5e-10:raise RuntimeError(f'cache {name} l{l} error {err}')
  p=out/f'l{l:04d}.npz';kk=np.load(ROOT/'results'/name/'final_states'/f'l{l:04d}.npz')['k']
  np.savez_compressed(p,coefficients=co,phase_offset=off,k=kk,l=l,degree=deg,center=cfg['center'],sigma=cfg['sigma'])
  records.append(dict(l=l,error=err,degree=deg,sha256=c.sha(p),source_cache_hash=meta['cache_hash']));print('HISTORY',name,l,err,flush=True)
 c.atomic_json(out/'MANIFEST.json',dict(config=cfg,records=records,max_error=max(v['error'] for v in records),new_evolutions=0,scope='Approximate interpolation at this same time window only; not a new mode campaign or universal history cache.'))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('name');run(a.parse_args().name)
