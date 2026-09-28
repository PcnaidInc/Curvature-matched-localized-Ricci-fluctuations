#!/usr/bin/env python3
"""Separate hash verification, independent checks, cached assembly and full computation."""
from pathlib import Path
import argparse,json,hashlib,subprocess,sys,tempfile,shutil
import numpy as np
ROOT=Path(__file__).resolve().parent
CAMPAIGNS=['outer','middle','inner','middle_time','outer_grid','inner_grid']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(script,args=(),work=ROOT):
 subprocess.run([sys.executable,str(work/'code'/script),*map(str,args)],check=True,cwd=work)
def blocks(root):
 records=0
 for name in CAMPAIGNS:
  path=root/'results'/name;status=json.loads((path/'STATUS.json').read_text());assert status['complete']
  if len(status['records'])!=status['config']['L']+1:raise ValueError(f'missing blocks in {name}')
  total=None
  for rec in status['records']:
   p=path/f"block_{rec['l']:04d}.npz"
   if sha(p)!=rec['sha256']:raise ValueError(f'corrupt block {p}')
   z=np.load(p);a=z['cov'];total=a.copy()if total is None else total+a;records+=1
  if not np.array_equal(total,np.load(path/'SPECTRA.npz')['covariance']):raise ValueError('sum mismatch')
 return records

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=['verify','checks','assemble','full']);p.add_argument('--workdir',type=Path);a=p.parse_args()
 if a.action=='verify':
  manifest=json.loads((ROOT/'FILE_HASHES.json').read_text())
  for name,h in manifest.items():
   if sha(ROOT/name)!=h:raise ValueError(f'hash mismatch {name}')
  print('Verified',len(manifest),'file hashes and',blocks(ROOT),'saved covariance blocks. No physics rerun.');return
 if a.action=='checks':
  with tempfile.TemporaryDirectory(prefix='bh_xi_checks_') as td:
   work=Path(td);shutil.copytree(ROOT/'code',work/'code');(work/'verification').mkdir();(work/'inputs').mkdir();(work/'inputs'/'prior').symlink_to(ROOT/'inputs/prior',target_is_directory=True)
   run('check_tensor.py',work=work)
   for name in ['check_pairs.py','check_gram.py','check_unreduced_trace.py']:run(name,[work/'inputs/prior'],work)
  print('Independent checks passed in a fresh directory. Six small reference modes, not full production.');return
 if a.workdir is None or a.workdir.exists():raise ValueError('Supply --workdir pointing to a NEW directory')
 work=a.workdir.resolve();work.mkdir(parents=True);shutil.copytree(ROOT/'code',work/'code');shutil.copytree(ROOT/'inputs',work/'inputs');(work/'results').mkdir();(work/'verification').mkdir()
 if a.action=='full':
  for name in CAMPAIGNS:
   script='coupling_fresh.py' if name.endswith('_grid') else 'coupling.py'
   run(script,[work/'inputs'/f'{name}.json','--source',work/'inputs/prior'],work)
 else:
  for name in CAMPAIGNS:
   shutil.copytree(ROOT/'results'/name,work/'results'/name,ignore=shutil.ignore_patterns('mode_cache'))
  blocks(work)
  for name in CAMPAIGNS:
   path=work/'results'/name;total=lc=kc=None
   for f in sorted(path.glob('block_*.npz')):
    z=np.load(f)
    if total is None:total=z['cov'].copy();lc=z['lcov'].copy();kc=z['kcov'].copy()
    else:total+=z['cov'];lc+=z['lcov'];kc+=z['kcov']
   np.savez_compressed(path/'SPECTRA.npz',covariance=total,truncated_angular=lc,truncated_momentum=kc,J=np.arange(len(total)))
 run('check_tensor.py',work=work)
 for name in ['check_pairs.py','check_gram.py','check_unreduced_trace.py']:run(name,[work/'inputs/prior'],work)
 run('analyze.py',[work/'inputs/prior'],work)
 d=json.loads((work/'results/FINAL_VALIDATION.json').read_text());assert d['all_primary_pass']
 print('Completed',a.action,'in',work)
if __name__=='__main__':main()
