"""Read-only independent checks of exact coverage, summation and radius brackets.
This verifies the saved interval certificate's accounting; full derivative replay is separate.
"""
from pathlib import Path
from fractions import Fraction as F
import json,sys
from mpmath import mp,iv

def fr(t):
 s,m,e,b=t;return (-1 if s else 1)*F(int(m))*F(2)**int(e)
def ep(a):return tuple(fr(t) for t in a['exact_dyadic'])
def ival(ab):
 a,b=ab
 return iv.mpf([(iv.mpf(a.numerator)/a.denominator).a,(iv.mpf(b.numerator)/b.denominator).b])

def verify(d:Path):
 mp.dps=100;iv.dps=90
 summary=json.loads((d/'summary.json').read_text());g=json.loads((d/'geometry.json').read_text());rows=g['rows']
 nt=197*summary['time_div']//64;nx=summary['directions'];dt=F(1,summary['time_div']);seen=set();sums=[F(0)]*nx
 def tor(r):return r+2*iv.ln(1-r/2)
 base=tor(iv.mpf(3)/2)-1;radius_checks=0
 for i,row in enumerate(rows):
  assert row['index']==i and F(*row['time_fraction'])==i*dt
  lo,hi=ep(row['radius']);integ=ival(ep(row['switch_integral']));target=base+integ
  assert tor(ival((lo,lo)))>target.b and tor(ival((hi,hi)))<target.a
  radius_checks+=1
 assert len(rows)==nt+1 and ep(rows[0]['radius'])[1]<F(7,4)
 for line in (d/'cells.jsonl').read_text().splitlines():
  r=json.loads(line);ij=(r['it'],r['ix']);assert ij not in seen;seen.add(ij)
  assert 0<=ij[0]<nt and 0<=ij[1]<nx
  assert ep(r['w'])[0]>0
  up=ep(r['rate_scaled'])[1];assert up>=0
  sums[ij[1]]+=dt*up
 assert len(seen)==nt*nx
 upper=max(sums)+F(2,2048)
 C=summary['C_integer'];assert C>=upper
 # The serialized interval-summed endpoint may round above the exact Fraction sum.
 assert ep(summary['C_interval'])[1]>=upper
 assert C==(ep(summary['C_interval'])[1].numerator+ep(summary['C_interval'])[1].denominator-1)//ep(summary['C_interval'])[1].denominator
 result={'certificate_directory':d.name,'interval_cells':len(seen),'certified_radius_endpoints_rechecked':radius_checks,
  'exact_direction_sums_and_ceiling_pass':True,'positive_widths':True,'gap_free_coverage':True,
  'C_integer':C,'derivative_cells_reexecuted':False,'science_field_evolutions':0}
 print(json.dumps(result,indent=2));return result
if __name__=='__main__':
 result=verify(Path(sys.argv[1]));Path(sys.argv[2]).write_text(json.dumps(result,indent=2))
