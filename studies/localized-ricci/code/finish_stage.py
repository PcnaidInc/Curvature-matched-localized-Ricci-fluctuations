from pathlib import Path
import json,time,subprocess,os
R=Path(__file__).resolve().parents[1]
env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
for _ in range(900):
 p=R/'results/r030_fine/STATUS.json'
 if p.exists():
  try:
   if json.loads(p.read_text()).get('complete'):break
  except json.JSONDecodeError:pass
 time.sleep(2)
else:raise TimeoutError('last production campaign not complete')
for args,label in [(['code/analyze.py'],'final_analysis'),(['code/save_histories.py','r030_fine'],'history_r030'),(['code/check_saved.py'],'saved_checks'),(['reproduce.py','checks'],'fresh_checks')]:
 with open(R/'logs'/f'{label}.log','w') as f:
  q=subprocess.run(['python',*args],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT)
 if q.returncode:raise RuntimeError(label+' failed; inspect log')
 print('FINISHED',label,flush=True)
(R/'logs/FINISH_STAGE.json').write_text(json.dumps({'completed':True,'commands':4}))
