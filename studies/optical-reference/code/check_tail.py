import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import erf
from optical_benchmark import Windows,rule
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a0=p.parse_args();a0.out.mkdir(parents=True,exist_ok=True)
w=Windows(n=256);res={}
for b in [2.,4.,8.]:
    vals=[]
    for lo,hi in [(0,16),(16,32),(32,64),(64,128)]:
        a,q=rule(lo,hi,480);f=w.Wflat(a)
        phase=(2*np.pi)**1.5/b**3*erf(a*b/(2*np.sqrt(2)))
        vals.append(float(np.sum(q*a**4*abs(f)**2*phase)/(4*np.pi**3)))
    res[str(b)]={'bands':[[0,16],[16,32],[32,64],[64,128]],'I_contributions':vals,'tail_to_base':sum(vals[1:])/vals[0]}
(a0.out/'high_frequency_tail.json').write_text(json.dumps(res,indent=2));print(json.dumps(res,indent=2))
# supplemental sensitivity, not an analytic infinite-tail error bound
raise SystemExit(0 if max(v['tail_to_base'] for v in res.values())<1e-6 else 2)
