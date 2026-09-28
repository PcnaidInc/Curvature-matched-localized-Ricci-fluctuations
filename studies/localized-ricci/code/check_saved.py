"""Integrity, deliberate resume, and corruption rejection; not a new mode campaign."""
from pathlib import Path
import hashlib,json,tempfile,shutil,subprocess,os
import joint_inward as j
import localized_trace as c
ROOT=Path(__file__).resolve().parents[1]

def run():
 records=[];num_modes=0;labels=0;pairblocks=0
 for cfgp in sorted((ROOT/'inputs').glob('*.json')):
  cfg=json.loads(cfgp.read_text())
  if not isinstance(cfg,dict) or 'name' not in cfg:continue
  folder=ROOT/'results'/cfg['name'];status=json.loads((folder/'STATUS.json').read_text());assert status['complete']
  for l in range(cfg['L']+1):
   p=folder/f'joint_{l:04d}.npz';rec=json.loads(p.with_suffix('.json').read_text());assert c.sha(p)==rec['sha256'];pairblocks+=1
   m=json.loads((folder/'mode_records'/f'l{l:04d}.json').read_text());assert c.sha(folder/'final_states'/f'l{l:04d}.npz')==m['final_hash'];num_modes+=1;labels+=m['labels']
  records.append(dict(name=cfg['name'],pairs=cfg['L']+1,labels=(cfg['L']+1)*cfg['nk']))
 before={str(p.relative_to(ROOT)):c.sha(p) for p in (ROOT/'results/pilot').glob('joint_*.npz')}
 env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
 q=subprocess.run(['python','code/joint_inward.py','inputs/pilot.json'],cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (ROOT/'logs/resume_test.log').write_text(q.stdout);assert q.returncode==0
 status=json.loads((ROOT/'results/pilot/STATUS.json').read_text());assert status['new_mode_labels_this_invocation']==0 and status['skipped_blocks']==32
 assert before=={str(p.relative_to(ROOT)):c.sha(p) for p in (ROOT/'results/pilot').glob('joint_*.npz')}
 oldj,oldc=j.ROOT,c.ROOT
 try:
  with tempfile.TemporaryDirectory() as temp:
   tmp=Path(temp);shutil.copytree(ROOT/'code',tmp/'code');folder=tmp/'results/pilot';folder.mkdir(parents=True)
   src=ROOT/'results/pilot/joint_0000.npz';shutil.copy(src,folder/src.name);shutil.copy(src.with_suffix('.json'),folder/'joint_0000.json')
   p=folder/src.name;b=bytearray(p.read_bytes());b[len(b)//2]^=1;p.write_bytes(b)
   j.ROOT=tmp;c.ROOT=tmp
   try:j.run(json.loads((ROOT/'inputs/pilot.json').read_text()));rejected=False
   except RuntimeError as e:rejected='block mismatch' in str(e)
   assert rejected
 finally:j.ROOT=oldj;c.ROOT=oldc
 result=dict(all_pass=True,production_mode_labels=labels,field_angular_records=num_modes,verified_pairblocks=pairblocks,resume_new_modes=0,resume_skipped_blocks=32,corruption_rejected=True,campaigns=records)
 (ROOT/'verification/RESTART_AND_HASH_CHECKS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':run()
