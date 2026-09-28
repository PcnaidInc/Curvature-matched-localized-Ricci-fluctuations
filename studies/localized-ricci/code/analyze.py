"""Assemble accepted data without modifying any field component or raw output."""
from pathlib import Path
import json,csv
import numpy as np
import joint_inward as j
import localized_trace as c
ROOT=Path(__file__).resolve().parents[1]
VECTORS={'trace':np.array([1.,0.]),'time_projection':np.array([0.,1.]),'energy_projection':np.array([-.5,1.])}

def spectrum(name,lcut=None,jcut=None):
 folder=ROOT/'results'/name;status=json.loads((folder/'STATUS.json').read_text())
 if not status['complete']:raise ValueError('incomplete '+name)
 cfg=status['config'];L=cfg['L'] if lcut is None else lcut;Jmax=cfg['Jmax'] if jcut is None else jcut
 if lcut is None and jcut is None:
  z=np.load(folder/'SPECTRA.npz');return z['covariance'],z['reference'],status
 out=np.zeros((Jmax+1,2,2));ref=np.zeros_like(out)
 for l in range(L+1):
  z=np.load(folder/f'joint_{l:04d}.npz')
  for p,values in zip(z['partners'],z['moments']):
   p=int(p)
   if p>L:continue
   for J in range(abs(l-p),min(l+p,Jmax)+1,2):
    factor=(2 if p!=l else 1)*(2*l+1)*(2*p+1)*(2*J+1)*c.wigner000_squared(l,p,J)/(16*np.pi**4)
    out[J]+=factor*j.evaluate(values[:6],J);ref[J]+=factor*j.evaluate(values[6:],J)
 return out,ref,status

def weight(J,delta):return (J==0).astype(float) if np.isinf(delta) else np.exp(-delta*delta*J*(J+1))
def projected(S,delta):return np.einsum('j,jab->ab',weight(np.arange(len(S)),delta),S)
def variance(S,delta,vec):return float(vec@projected(S,delta)@vec)

def main():
 checks={};rows=[];spectral=[];candidates=[];gates={};minimum_eig=1.;bg=c.Background(.2)
 old=np.genfromtxt(ROOT/'inputs/ANGULAR_POWER.csv',names=True,delimiter=',')
 for tag in ('r060','r045','r030'):
  name=tag+'_fine';S,F,status=spectrum(name);cfg=status['config'];rc=cfg['center'];ds=[np.inf]+cfg['angular_deltas']
  Cm,Fm,sm=spectrum(tag+'_main');matched,matchedF,_=spectrum(name,sm['config']['L']);
  priorcut=sm['config']['L'];ext24,ref24,_=spectrum(name,jcut=24)
  checks[tag]={};
  for obs,vec in VECTORS.items():
   vf=np.array([variance(S,d,vec) for d in ds]);vmatch=np.array([variance(matched,d,vec) for d in ds]);vm=np.array([variance(Cm,d,vec) for d in ds]);v24=np.array([variance(ext24,d,vec) for d in ds]);
   checks[tag][obs]={'last_angular_fraction':float(np.max(abs(vf/vmatch-1))),
     'matched_L_momentum_time_tolerance_fraction':float(np.max(abs(vmatch/vm-1))),
     'external_J24_to32_fraction':float(np.max(abs(vf/v24-1)))}
   for key,val in checks[tag][obs].items():gates[tag+'_'+obs+'_'+key]=val<(.001 if key.startswith('external') else .01)
   rf=np.array([variance(F,d,vec) for d in ds]);rfmatch=np.array([variance(matchedF,d,vec) for d in ds]);rfmain=np.array([variance(Fm,d,vec) for d in ds])
   checks[tag][obs]['reference_cutoff_fraction']=float(max(np.max(abs(rf/rfmatch-1)),np.max(abs(rfmatch/rfmain-1))))
   gates[tag+'_'+obs+'_reference']=checks[tag][obs]['reference_cutoff_fraction']<.01
   for i,d in enumerate(ds):
    mat=projected(S,d);ev=np.linalg.eigvalsh(mat);minimum_eig=min(minimum_eig,float(ev[0]/ev[-1]));rlo,rhi=status['support'];Ac=2/rc-1
    scales={'time_width':cfg['sigma'],'min_longitudinal_width':cfg['bproper']*np.sqrt((2/rhi-1)/Ac),'min_curvature_length':rlo**1.5/48**.25,'minimum_radius':rlo}
    if np.isfinite(d):scales['min_angular_proper_width']=rlo*d
    minscale=min(scales.values());jj=np.arange(2000);area=1. if np.isinf(d) else 1/np.sum((2*jj+1)*np.exp(-d*d*jj*(jj+1)))
    coeff=8*np.pi*np.sqrt(vf[i])/status['mean_sqrtK'];eps=1e-9
    row=dict(center=rc,observable=obs,delta_rad=0. if np.isinf(d) else d,uniform_sphere=int(np.isinf(d)),effective_sphere_fraction=area,
      sigma_M=cfg['sigma'],longitudinal_M=cfg['bproper'],angular_proper_center_M=0. if np.isinf(d) else rc*d,
      curvature_length_center_M=rc**1.5/48**.25,support_min_r=rlo,support_max_r=rhi,
      variance=vf[i],reference_variance=rf[i],source_RMS=np.sqrt(vf[i]),reference_RMS=np.sqrt(rf[i]),source_to_reference_RMS=np.sqrt(vf[i]/rf[i]),
      mean_background_sqrtK=status['mean_sqrtK'],relative_curvature_RMS_per_epsilon=coeff,epsilon_example=eps,percent_example=100*eps*coeff,
      minimum_screen_length_M=minscale,time_width_Planck_example=cfg['sigma']/np.sqrt(eps),
      minimum_separation_at_unit_ratio=minscale*np.sqrt(coeff),time_width_Planck_at_formal_unit=cfg['sigma']*np.sqrt(coeff))
    rows.append(row)
    for n in (10,30,100):candidates.append(dict(center=rc,observable=obs,delta_rad=row['delta_rad'],N=n,eps_cap=(minscale/n)**2,max_percent=100*(minscale/n)**2*coeff))
  # Cross-covariance full matrix sensitivity normalized by std product rather than a zero element.
  covtests=[]
  for d in ds:
   C=projected(S,d);M=projected(matched,d);N=projected(Cm,d);scale=np.sqrt(np.diag(C)[:,None]*np.diag(C)[None,:]);covtests.append(max(float(np.max(abs(C-M)/scale)),float(np.max(abs(M-N)/scale))))
  checks[tag]['joint_matrix_normalized_fraction']=max(covtests);gates[tag+'_joint_matrix']=max(covtests)<.01
  if tag=='r060':
   for obs,col,idx in [('trace','trace_power',0),('time_projection','time_projection_power',1)]:
    rr=[]
    for d in ds:
     vo=float(weight(old['J'],d)@old[col]);vn=variance(S,d,VECTORS[obs]);rr.append(abs(vn/vo-1))
    checks['overlap_'+obs+'_relative']=max(rr);gates['overlap_'+obs]=max(rr)<.01
  for J in range(len(S)):spectral.append(dict(center=rc,J=J,trace=S[J,0,0],time_projection=S[J,1,1],cross=S[J,0,1],ref_trace=F[J,0,0],ref_time=F[J,1,1],ref_cross=F[J,0,1]))
 ind=json.loads((ROOT/'verification/INDEPENDENT_CHECKS.json').read_text());gates['independent_checks']=ind['all_pass'];gates['positive_covariance']=minimum_eig>=-1e-12
 # Previously rejected pilot retained, no gate requiring acceptance.
 p,pr,ps=spectrum('pilot');so,_,_=spectrum('r060_fine');checks['pilot_variance_error_max']=max(abs(variance(p,d,v)/variance(so,d,v)-1) for d in [.25,.5,np.inf] for v in VECTORS.values())
 for fn,data in [('INWARD_LOCALIZED_RESULTS.csv',rows),('ANGULAR_SPECTRA.csv',spectral),('CONDITIONAL_SCALE_SCREENS.csv',candidates)]:
  with open(ROOT/'results'/fn,'w') as f:w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
 out=dict(status='accepted bounded result' if all(gates.values()) else 'has_failed_gates',all_primary_pass=all(gates.values()),checks=checks,gates=gates,rows=rows,minimum_normalized_covariance_eigenvalue=minimum_eig,scale_screens=candidates,
    scope='leading Ricci/Einstein projections, SPV-1 state, whole sphere and heat kernels with all primary center scales matched; no full metric or endpoint',epsilon_example=1e-9)
 (ROOT/'results/FINAL_VALIDATION.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k not in ('rows','scale_screens')},indent=2));print('ROWS',json.dumps(rows,indent=2))
if __name__=='__main__':main()
