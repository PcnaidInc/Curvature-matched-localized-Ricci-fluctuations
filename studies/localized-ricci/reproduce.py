"""verify: hashes only; checks: new small independent tests; full: all new modes."""
from pathlib import Path
import argparse,hashlib,json,subprocess,tempfile,shutil,os
ROOT=Path(__file__).resolve().parent
ENV={**os.environ,**{k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def run(args,cwd):subprocess.run(['python',*args],cwd=cwd,env=ENV,check=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['verify','checks','assemble','full']);p.add_argument('--workdir');a=p.parse_args()
 if a.mode=='verify':
  manifest=json.loads((ROOT/'FILE_HASHES.json').read_text())
  for f,h in manifest.items():
   if sha(ROOT/f)!=h:raise RuntimeError('checksum mismatch: '+f)
  print('VERIFIED',len(manifest),'files; no scientific calculation rerun');return
 if a.mode=='checks':
  with tempfile.TemporaryDirectory() as d:
   d=Path(d);shutil.copytree(ROOT/'code',d/'code');shutil.copytree(ROOT/'inputs',d/'inputs');(d/'verification').mkdir();run(['code/check_independent.py'],d)
   q=json.loads((d/'verification/INDEPENDENT_CHECKS.json').read_text());assert q['all_pass'];print('INDEPENDENT CHECKS PASSED; no production modes rerun')
  return
 if a.mode=='assemble':
  run(['code/analyze.py'],ROOT);return
 if not a.workdir:raise ValueError('full requires --workdir NEW_EMPTY_DIRECTORY')
 d=Path(a.workdir).resolve()
 if d==ROOT or d.exists():raise ValueError('choose a new work directory; existing results will not be overwritten')
 d.mkdir();shutil.copytree(ROOT/'code',d/'code');shutil.copytree(ROOT/'inputs',d/'inputs')
 for name in ['results','logs','verification','figures']:(d/name).mkdir()
 run(['code/check_independent.py'],d)
 for name in ['pilot','r060_main','r060_fine','r045_main','r045_fine','r030_main','r030_fine']:
  run(['code/joint_inward.py',f'inputs/{name}.json'],d)
 run(['code/analyze.py'],d)
 print('FULL mode campaigns completed; inspect FINAL_VALIDATION.json rather than assuming all gates passed.')
if __name__=='__main__':main()
