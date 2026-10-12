"""Independent algebra/endpoint checks; no physical field evolution."""
from pathlib import Path
import json,sys,math
from fractions import Fraction
import sympy as sp
from mpmath import mp,iv
from preparation_enclosure import *

def run(out:Path):
    if out.exists():raise FileExistsError(out)
    mp.dps=100;iv.dps=70
    t,e=sp.symbols('t e');a=sp.Function('a')(t);b=sp.Function('b')(t);cc=sp.Function('c')(t)
    F0=-sp.diff(a,t,2)/(2*a)+3*sp.diff(a,t)**2/(4*a**2)
    F1=-(sp.diff(b,t,2)/a-sp.diff(a,t,2)*b/a**2)/2+sp.Rational(3,2)*(sp.diff(a,t)/a)*(sp.diff(b,t)/a-sp.diff(a,t)*b/a**2)
    U=2*a*b-F0;w=a+e*b+e**2*cc
    N=(a*a+e*U)*w*w-w**4+e*(-sp.diff(w,t,2)*w/2+sp.Rational(3,4)*sp.diff(w,t)**2)
    poly=sp.Poly(sp.expand(N),e)
    checks={}
    for j in range(3):
        x=sp.factor(poly.nth(j).subs(cc,(F1-b*b)/(2*a)))
        checks['residual_coefficient_'+str(j)]=bool(x==0)
    # Canonical vbar has the advertised oscillator defect.
    W=sp.Function('W')(t);vlog=-sp.diff(W,t)/(2*W)-sp.I*W
    actual=sp.diff(vlog,t)+vlog**2+W**2
    expect=-sp.diff(W,t,2)/(2*W)+3*sp.diff(W,t)**2/(4*W**2)
    checks['comparison_mode_defect_identity']=bool(sp.simplify(actual-expect)==0)
    # Verify logistic-jet derivatives against direct high-precision differentiation.
    sm=lambda x:1/(1+mp.exp(2/x-2/(2-x)))
    endpoint=[];interior=[]
    for t0 in [mp.mpf(1)/256,mp.mpf(1)/128,mp.mpf(1)/64]:
        ej=endpoint_jets(False)
        for j in range(9):
            value=mp.diff(sm,t0,j)/mp.factorial(j)
            ok=mp.make_mpf(ej[j]._mpi_[0])<=value<=mp.make_mpf(ej[j]._mpi_[1])
            endpoint.append({'t':str(t0),'order':j,'value':mp.nstr(value,25),'enclosed':bool(ok)})
    for t0 in [mp.mpf(1)/4,mp.mpf(3)/4,mp.mpf(5)/4,mp.mpf(7)/4]:
        tjj=switch_jet(iv.mpf(t0),iv.mpf(t0))
        for j in range(9):
            value=mp.diff(sm,t0,j)/mp.factorial(j)
            lo,hi=tj.bounds(tjj[j]);ok=lo<=value<=hi
            interior.append({'t':str(t0),'order':j,'enclosed':bool(ok)})
    checks['endpoint_derivative_samples_inside_cauchy_bounds']=all(a['enclosed'] for a in endpoint)
    checks['interior_jets_vs_independent_differentiation']=all(a['enclosed'] for a in interior)
    # Analytic Cauchy denominator and monotonic derivative-bound factors, exact rationals.
    delta=Fraction(1,128)
    checks['cauchy_monotonic_for_all_retained_derivatives']=all(delta<Fraction(2,3*j) for j in range(1,9))
    checks['cauchy_disk_inside_switch_domain']=Fraction(3,2)*delta<1
    eta=iv.exp(1/(1-rational(3,2)*rational(1,128))-2/(3*rational(1,128)))
    checks['cauchy_denominator_ge_half']=bool(eta<rational(1,2))
    # Initial mixing measured in an exact constant-frequency basis, not zeroed.
    initial=[]
    for D in [mp.mpf(16)/49,mp.mpf('0.5'),mp.mpf(1)]:
        for k in [2048,4096,8192]:
            eps=mp.mpf(k)**-2
            wa=mp.sqrt(D+eps);wb=mp.sqrt(D)+eps/(2*mp.sqrt(D))-eps**2/(8*D**mp.mpf('1.5'))
            rapidity=abs(mp.log(wa/wb))/2
            initial.append({'D':str(D),'kappa':k,'nonzero_mismatch':mp.nstr(rapidity,30),'below_bound':bool(rapidity<mp.mpf(2)/k**6)})
    checks['initial_sqrt_mismatch_checks']=all(x['below_bound'] for x in initial)
    # Deliberately bad initial-data omission: constant oscillator, different initial width.
    checks['zero_initial_mismatch_control_detected']=all(mp.mpf(x['nonzero_mismatch'])>0 for x in initial)
    result={'checks':checks,'endpoint_samples':endpoint,'interior_samples':interior,'initial_samples':initial,
      'status':'algebra and implementation diagnostics, not independent SPV-1 mode runs',
      'field_mode_evolutions':0,'all_pass':all(checks.values())}
    out.write_text(json.dumps(result,indent=2));print(json.dumps(checks,indent=2));return result

if __name__=='__main__':run(Path(sys.argv[1]))
