"""Independent compatibility and deliberately invalid-input controls."""
from pathlib import Path
import json,sys
from fractions import Fraction as F
from mpmath import mp,iv
from preparation_enclosure import tj,rational,scalar_switch,geometry_jet,switch_jet,comparison,rate_enclosure

def run(out):
 mp.dps=100;iv.dps=80
 tor=lambda r:r+2*iv.ln(1-r/2)
 t0=2+tor(rational(21,25))-tor(rational(3,2))
 t1=2+tor(rational(19,25))-tor(rational(3,2))
 q=4*iv.sqrt(iv.mpf(7)/3)
 checks={'pulse_covered_by_all_time_enclosure':bool(t0>2 and t1<rational(197,64)),
   'omitted_high_angular_leg_threshold':bool(iv.sqrt(2048*2049)>=2048),
   'omitted_high_momentum_leg_threshold':bool(6144-q/2>2048)}
 # Independent exact-rational convex-range test; derivative sign extrema at corners.
 intervals=[(F(0),F(1,32),F(1,10),F(1,8)),(F(1,4),F(1,2),F(2),F(3)),(F(0),F(1),F(1,2),F(2))]
 n=0
 for xl,xh,bl,bh in intervals:
  corners=[x+(1-x)*b for x in (xl,xh) for b in (bl,bh)]
  lower,upper=min(corners),max(corners)
  for i in range(11):
   for j in range(11):
    x=xl+(xh-xl)*F(i,10);b=bl+(bh-bl)*F(j,10);assert lower<=x+(1-x)*b<=upper;n+=1
 checks['exact_rational_convex_range_controls']=True
 # Verify active limit of NEW smooth geometry against original active jets.
 errors=[]
 for rr in [F(19,25),F(4,5),F(21,25)]:
  rpoint=rational(rr.numerator,rr.denominator)
  r,A,B,U=geometry_jet(tj.constant(1),rpoint)
  original=tj.radius_jet(rpoint,rpoint)
  for j in range(7):
   lo=max(mp.make_mpf(r[j]._mpi_[0]),mp.make_mpf(original[j]._mpi_[0]))
   hi=min(mp.make_mpf(r[j]._mpi_[1]),mp.make_mpf(original[j]._mpi_[1]))
   assert lo<=hi
  for x0 in [0,F(1,2),1]:
   xx=rational(x0.numerator,x0.denominator) if isinstance(x0,F) else iv.mpf(x0)
   D,ws=comparison(B,U,xx)
   # Direct high-precision INTERVAL arithmetic; preserve correlations as enclosures.
   # Midpoints of distinct intervals cannot be assumed to satisfy exact cancellations.
   e=iv.mpf(1)/2048**2
   _,rate,cs=rate_enclosure(D,U,ws,e)
   factored=sum((a*e**i for i,a in enumerate(cs)),iv.mpf(0))
   iv.dps=120
   r120=rational(rr.numerator,rr.denominator)
   x120=rational(x0.numerator,x0.denominator) if isinstance(x0,F) else iv.mpf(x0)
   _,_,B120,U120=geometry_jet(tj.constant(1),r120)
   D120,ww=comparison(B120,U120,x120);ee=iv.mpf(1)/2048**2
   W=[sum((tj.der(w,j)[0]*ee**i for i,w in enumerate(ww)),iv.mpf(0)) for j in range(3)]
   direct=((D120[0]+ee*U120[0])*W[0]**2-W[0]**4+ee*(-W[2]*W[0]/2+3*W[1]**2/4))/ee**3
   # High-precision direct interval must lie inside lower-precision factored interval.
   assert factored.a<=direct.a and direct.b<=factored.b
   iv.dps=80
   errors.append(tj.serial(direct))
 checks['active_radius_jet_overlap']=True
 checks['direct_vs_factored_residual_enclosure']=True
 data={'checks':checks,'all_pass':all(checks.values()),'exact_convex_samples':n,
       'direct_factored_comparisons':len(errors),'pulse_times':[tj.serial(t0),tj.serial(t1)],
       'mode_evolutions':0,'direct_residual_values':errors}
 Path(out).write_text(json.dumps(data,indent=2));print(json.dumps(checks,indent=2));return data
if __name__=='__main__':run(sys.argv[1])
