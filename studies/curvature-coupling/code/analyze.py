from pathlib import Path
import json,csv,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
XI=[0.,1/12,1/6,1/4]
def projected(cov,xi):
 a=1-6*xi;b=6*xi
 tr=np.array([[a,0,0],[0,a,b]])
 return tr@cov@tr.T

def main(source):
 rows=[];checks={};gates={};allcounts=[]
 pairs={'outer':'outer_grid','middle':'middle_time','inner':'inner_grid'}
 for name in ('outer','middle','inner'):
  path=ROOT/'results'/name
  if not (path/'SPECTRA.npz').exists():continue
  stat=json.loads((path/'STATUS.json').read_text());cfg=stat['prior_config'];z=np.load(path/'SPECTRA.npz');J=z['J']
  prev=np.load(source/'results'/stat['config']['prior']/'SPECTRA.npz')['covariance']
  checks[name]={};testpath=ROOT/'results'/pairs[name]/'SPECTRA.npz';refine=np.load(testpath)['covariance'] if testpath.exists() else None
  for delta in [0]+cfg['angular_deltas']:
   w=(J==0).astype(float) if delta==0 else np.exp(-delta*delta*J*(J+1))
   c=np.einsum('j,jab->ab',w,z['covariance']);lc=np.einsum('j,jab->ab',w,z['truncated_angular']);kc=np.einsum('j,jab->ab',w,z['truncated_momentum']);jc=np.einsum('j,jab->ab',w*(J<=24),z['covariance'])
   old=np.einsum('j,jab->ab',w,prev)
   ov=float(np.max(abs(c[:2,:2]-old)/np.sqrt(np.outer(np.diag(old),np.diag(old)))))
   gates[f'{name}_{delta}_baseline']=ov<1e-6
   # compare 3-source covariance on its own diagonals; tiny eig supported at covariance parent scale
   eig=np.linalg.eigvalsh(c);gates[f'{name}_{delta}_PSD']=eig.min()>=-1e-12*np.diag(c).max()
   for xi in XI:
    q=projected(c,xi);qL=projected(lc,xi);qK=projected(kc,xi);qJ=projected(jc,xi)
    nonzero=np.diag(q)>1e-20*max(np.diag(c));inds=np.where(nonzero)[0]
    denom=np.sqrt(np.outer(np.diag(q)[inds],np.diag(q)[inds]))
    err=lambda x:float(np.max(abs(x[np.ix_(inds,inds)]-q[np.ix_(inds,inds)])/denom))
    lerr,kerr,jerr=err(qL),err(qK),err(qJ)
    ferr=err(projected(np.einsum('j,jab->ab',w,refine),xi)) if refine is not None else None
    key=f'{name}_d{delta}_xi{xi}'
    checks[name][key]=dict(baseline_overlap=ov,angular_covariance=lerr,momentum_tail_covariance=kerr,external_harmonic_covariance=jerr,independent_grid_covariance=ferr)
    gates[key+'_angular']=lerr<.01;gates[key+'_momentum_tail']=kerr<.01;gates[key+'_J']=jerr<.001;gates[key+'_grid']=ferr is not None and ferr<.001
    trsd=np.sqrt(max(q[0,0],0));qsd=np.sqrt(max(q[1,1],0));D=stat['mean_sqrtK'];a=1-6*xi
    rows.append(dict(center=cfg['center'],delta_rad=delta,xi=xi,trace_rms=trsd,observer_source_rms=qsd,
      observer_rms_vs_minimal=qsd/np.sqrt(c[1,1]),trace_rms_vs_minimal=abs(a),source_cross_covariance=q[0,1],
      scalar_curvature_rms_per_epsilon=8*np.pi*trsd/D,observer_Ricci_rms_per_epsilon=8*np.pi*qsd/D,
      scalar_percent_at_epsilon1e9=100*1e-9*8*np.pi*trsd/D,observer_percent_at_epsilon1e9=100*1e-9*8*np.pi*qsd/D,
      mean_background_sqrtK=D,minimal_variance_T=c[0,0],minimal_variance_P=c[1,1],conformal_observer_variance=c[2,2],
      minimal_conformal_covariance=c[1,2],reference_state=cfg['seed']))
 for p in sorted((ROOT/'results').glob('*/STATUS.json')):
  d=json.loads(p.read_text())
  if d.get('complete'):
   rec=d.get('records',[]);allcounts.append(dict(campaign=p.parent.name,pair_blocks=len(rec),pair_evaluations=sum(v['pair_count']for v in rec),summed_block_seconds=sum(v['seconds']for v in rec),new_field_labels=d.get('field_spectral_labels',0),scope='new field/grid'if d.get('field_spectral_labels',0)else'cached phase-history contraction'))
 for fn in ['TENSOR_CHECK.json','PAIR_CHECKS.json','GRAM_CHECK.json','UNREDUCED_TRACE.json']:
  p=ROOT/'verification'/fn
  gates[fn]=p.exists() and json.loads(p.read_text()).get('all_pass',False)
 gates={key:bool(val)for key,val in gates.items()}
 if rows:
  with (ROOT/'results/COUPLING_COMPARISON.csv').open('w')as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 out=dict(all_primary_pass=bool(rows) and len(rows)==36 and all(gates.values()),gates=gates,checks=checks,execution=allcounts,
    accepted_scope='fixed Cauchy state on final Ricci-flat Schwarzschild, trace and one trace-reversed observer projection; not a different collapse preparation or full metric',
    notes='new nontrace contractions; trace dependence exact. Grid comparison is fresh time/momentum/ODE at outer and inner, history time-grid only middle. Cutoff tails are not rigorous infinite-tail bounds.')
 (ROOT/'results/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items()if k not in ['checks','gates']},indent=2));print('FAILED', [k for k,v in gates.items()if not v])
if __name__=='__main__':main(Path(sys.argv[1]))
