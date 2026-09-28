import subprocess,time,json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
queue=['pilot','r060_main','r060_fine','r045_main','r045_fine','r030_main','r030_fine'];jobs={};done=[]
env={**os.environ,**{key:'1' for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}}
while queue or jobs:
 while queue and len(jobs)<2:
  name=queue.pop(0);log=open(ROOT/'logs'/f'{name}.log','a');start=time.time()
  p=subprocess.Popen(['python','code/joint_inward.py',f'inputs/{name}.json'],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
  jobs[name]=(p,log,start);print('LAUNCH',name,p.pid,flush=True)
 for name,(p,log,start) in list(jobs.items()):
  rc=p.poll()
  if rc is not None:
   log.close();done.append(dict(name=name,exitcode=rc,elapsed=time.time()-start));del jobs[name]
   print('EXIT',name,rc,flush=True)
 (ROOT/'logs/QUEUE.json').write_text(json.dumps(dict(done=done,running={n:p.pid for n,(p,l,t) in jobs.items()},waiting=queue),indent=2))
 time.sleep(2)
print('FINISHED',flush=True)
