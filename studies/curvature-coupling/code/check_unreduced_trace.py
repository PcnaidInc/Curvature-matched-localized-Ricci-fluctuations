"""Compute conformal trace from the unreduced differential tensor, not the null rule."""
from pathlib import Path
import sys,json
import numpy as np
import coupling as c
ROOT=Path(__file__).resolve().parents[1]
def main(source):
 checks=[]
 for prior in ['r060_fine','r045_fine','r030_fine']:
  cfg=json.loads((source/'inputs'/f'{prior}.json').read_text());bg=c.core.Background(.2);w=c.core.window(bg,cfg)
  r,A,rp,*_=bg.at(w['t']).T;x=w['tau']/(5*cfg['sigma'])
  for l,p,i,j,J in [(0,0,1,2,0),(3,4,13,17,3),(70,73,65,71,5)]:
   f,d,k=c.load_hist(source,prior,l,x,r,rp);ff,dd,kk=c.load_hist(source,prior,p,x,r,rp)
   fi,di,fj,dj=f[:,i],d[:,i],ff[:,j],dd[:,j];a,b=k[i],kk[j];q=fi*fj;dq=di*fj+fi*dj
   ddq=2*di*dj-2*rp/r*dq-(a*a+b*b+A*(l*(l+1)+p*(p+1))/r**2)*q
   for sg in[-1,1]:
    boxq=-(ddq+2*rp/r*dq)/A-((a+sg*b)**2/A+J*(J+1)/r**2)*q
    trace0=(di*dj+sg*a*b*q)/A-(l*(l+1)+p*(p+1)-J*(J+1))*q/(2*r**2)
    raw=trace0+0.5*boxq
    err=float(np.max(abs(raw))/max(np.max(abs(trace0)),1e-100))
    checks.append(dict(prior=prior,l=l,p=p,sign=sg,error=err))
 out=dict(checks=checks,max_relative_parent=max(z['error']for z in checks),all_pass=max(z['error']for z in checks)<1e-8,
          new_mode_evolutions=0,scope='Unreduced field-pair differential trace on actual same-window histories; not a second proof of anomaly renormalization.')
 (ROOT/'verification/UNREDUCED_TRACE.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));
 if not out['all_pass']:raise SystemExit(1)
if __name__=='__main__':main(Path(sys.argv[1]))
