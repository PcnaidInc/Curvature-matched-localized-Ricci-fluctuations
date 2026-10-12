"""Portable evidence checks and explicitly NEW bounded calculations."""
from pathlib import Path
import argparse,json,hashlib,subprocess,sys
R=Path(__file__).resolve().parent

def new_out(p):
 p=Path(p).resolve()
 if p.exists() or p==R or R in p.parents:raise ValueError('Use a NEW directory/file outside the evidence package')
 p.parent.mkdir(parents=True,exist_ok=True);return p

def main():
 ap=argparse.ArgumentParser();sp=ap.add_subparsers(dest='mode',required=True)
 sp.add_parser('verify');sp.add_parser('certificate-check')
 red=sp.add_parser('reduce-saved');red.add_argument('--out',required=True,type=Path)
 ck=sp.add_parser('algebra-checks');ck.add_argument('--out',required=True,type=Path)
 prep=sp.add_parser('preparation');prep.add_argument('--out',required=True,type=Path);prep.add_argument('--time-div',type=int,default=64);prep.add_argument('--directions',type=int,default=8);prep.add_argument('--integral-div',type=int,default=8192);prep.add_argument('--dps',type=int,default=40)
 a=ap.parse_args()
 if a.mode=='verify':
  m=json.loads((R/'MANIFEST.json').read_text())
  for f in m['files']:
   p=(R/f['path']).resolve()
   if R not in p.parents:raise ValueError('unsafe manifest path')
   assert p.stat().st_size==f['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==f['sha256'],f['path']
  print(json.dumps({'payloads_verified':len(m['files']),'physics_reexecuted':False,'physical_scalar_accepted':False}));return
 if a.mode=='reduce-saved':
  out=new_out(a.out);sys.path.insert(0,str(R/'code'));from assemble_results import assemble
  out.write_text(json.dumps(assemble(R),indent=2));return
 if a.mode=='certificate-check':
  sys.path.insert(0,str(R/'code'));from verify_certificate import verify
  verify(R/'results/prep_final70');return
 if a.mode=='algebra-checks':
  out=new_out(a.out);subprocess.run([sys.executable,str(R/'code/check_preparation.py'),str(out)],check=True);return
 if a.mode=='preparation':
  out=new_out(a.out);subprocess.run([sys.executable,str(R/'code/preparation_enclosure.py'),'--out',str(out),'--time-div',str(a.time_div),'--directions',str(a.directions),'--integral-div',str(a.integral_div),'--dps',str(a.dps)],check=True)
if __name__=='__main__':main()
