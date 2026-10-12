"""Conditional state-width-replacement allowance.
The R13 preparation constant is a reported INPUT, not recertified by this program.
It NEVER returns an accepted physical scalar.
"""
from __future__ import annotations
import argparse,json,time
from pathlib import Path
from mpmath import mp,iv
from active_enclosure import serial,bounds

C_REPORTED=20547271

def pulse_integral(n):
    """Monotone Darboux bounds on integral_0^1 exp(1-1/(1-z^2)) dz."""
    middle=iv.mpf(0)
    for j in range(1,n):
        z=iv.mpf(j)/n
        middle+=iv.exp(1-1/(1-z*z))
    lower=middle/n;upper=(middle+1)/n
    return iv.mpf([lower.a,upper.b])

def spectral_factor(N=2048,V=6144,order=16):
    n0=iv.mpf(N)+iv.mpf(1)/2; L0=n0*n0-iv.mpf(1)/4
    # Exact v integration followed by decreasing first-term + integral bound.
    angular=3/(32*iv.pi)*(n0/L0**(iv.mpf(5)/2)+1/(3*L0**(iv.mpf(3)/2)))
    terms=[]
    for j in range(order+1):
        integer_moment=sum((2*l+1)*(l*(l+1))**j for l in range(N))
        terms.append(iv.mpf((-1)**j*(j+1)*(j+2)*integer_moment)/(2*(5+2*j)*V**(5+2*j))/(4*iv.pi**2))
    assert order%2==0 and iv.mpf(15)/7*((N-1)*N)/V**2<1
    upper=sum(terms,iv.mpf(0));lower=sum(terms[:-1],iv.mpf(0))
    tail=iv.mpf([lower.a,upper.b])
    crude=iv.mpf(N)**2/(20*iv.pi**2*V**5)
    return angular,tail,crude

def run(cert_path:Path,nchi=16384,dps=55):
    start=time.perf_counter();mp.dps=dps+25;iv.dps=dps
    cert=json.loads(cert_path.read_text())
    if not all(cert['checks'].values()):raise ValueError('Active enclosures not accepted')
    p={k:iv.mpf(v) for k,v in cert['rational_promises'].items()}
    N=2048;V=6144;rho0=iv.sqrt(N*(N+1));q=4*iv.sqrt(iv.mpf(7)/3)
    c=1-q/(2*rho0);d=1+q/(2*rho0)
    Econst=2*C_REPORTED/c**5*iv.exp(2*C_REPORTED/(c*rho0)**5)
    b0=2*C_REPORTED/(c*rho0)**5
    ub=p['ubar_abs_upper'];ua=ub+2*C_REPORTED*p['w_upper']*iv.exp(iv.mpf(2*C_REPORTED)/2048**5)/2048**4
    m=p['w_lower']*c;ma=m*iv.exp(-b0);Mhi=p['w_upper']*d
    rlo=iv.mpf(19)/25;rhi=iv.mpf(21)/25
    tor=lambda r:r+2*iv.ln(1-r/2)
    duration=tor(rlo)-tor(rhi);chi_average=pulse_integral(nchi);I=duration*chi_average
    r2max=1/rlo**2;Q=p['H_abs_upper']/3;Qprime=p['Hprime_abs_upper']/3
    P=p['B_upper']+(5*q*q/6+p['B_upper']/3)/rho0**2
    Pprime=p['Bprime_abs_upper']*(1+1/(3*rho0**2))
    Sint=P*I+2*Q/rho0**2
    Ssup=P+6*Q/(duration*rho0**2)
    Svar=2*P+Pprime*I+(400*Q/duration+2*Qprime)/rho0**2
    Ycoef=(Ssup+Svar+4*ub*Sint)/(8*m*m)
    # g positive, Gaussian <=1; its r^-2 factor is not frozen.
    Gint=2*r2max;Gsup=2*r2max/I;Gvar=4*r2max/I+2*p['R2prime_abs_upper']
    Lamcoef=(Gsup+Gvar+4*ua*Gint)/(2*ma*ma)
    parts={'source_change':Econst*Lamcoef*Sint/(4*m),
           'frequency_change':Econst*Lamcoef*duration*2*Mhi*Ycoef,
           'readout_change':Econst*Gint*Ycoef/m}
    K6=sum(parts.values(),iv.mpf(0));ang,tail,crude=spectral_factor(N,V)
    total=K6*(ang+tail)
    old=iv.mpf('0.05835854031693310')
    consts=dict(rho_min=rho0,qmax=q,leg_lower_ratio=c,leg_upper_ratio=d,
      actual_width_relative_coefficient=Econst,max_log_ratio=b0,bar_rate=ub,actual_rate=ua,
      duration=duration,pulse_average=chi_average,pulse_integral=I,bar_width_lower=m,
      actual_width_lower=ma,bar_width_upper=Mhi,P=P,Pprime=Pprime,
      S_integral_scaled=Sint,S_max_scaled=Ssup,S_variation_scaled=Svar,
      normalized_response_bound=Ycoef,normalized_adjoint_rho2_bound=Lamcoef,
      Gint=Gint,Gvar=Gvar,K6=K6,high_angular_factor_upper=ang,
      low_angular_momentum_factor=tail,low_angular_crude_upper=crude,
      omitted_domain_factor_upper=ang+tail,integrated_component_allowance=total,
      recorded_old_allowance=old,improvement_vs_recorded=old/total)
    # Checkability gate, NOT acceptance of physical target.
    return {'run':'R93R13C','dps':dps,'nchi':nchi,'seconds':time.perf_counter()-start,
      'premise':{'C':C_REPORTED,'source':'PCNBR-109 comment13816','original_enclosure_archive_available':False,
      'preparation_certificate_independently_reexecuted':False,
      'mathematical_status':'CONDITIONAL on the reported smooth-state mode bound and its log-derivative corollary'},
      'constants':{k:serial(v) for k,v in consts.items()},'K6_components':{k:serial(v) for k,v in parts.items()},
      'checks':{'active_enclosures_pass':True,'positive_leg_lower':bool(c>0),
       'summable_power_six':True,'tail_bracket_below_crude':bool(tail.b<crude.a),
       'initial_scalar_response_is_zero_by_compact_pulse':True},
      'conditional_component_bound_upper':mp.nstr(bounds(total)[1],35),
      'physical_scalar_accepted':False,'physical_scalar_values':None,'complete_error_bound':None,
      'unbounded_terms':['comparison-width response minus R12 local-order9 trial',
       'original smooth-target preparation proof not revalidated in this release',
       'inherited numerical geometry, integration, finite-sum global errors'],
      'no_new_field_mode_evolution':True}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--certificate',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--nchi',type=int,default=16384);p.add_argument('--dps',type=int,default=55);a=p.parse_args()
 if a.out.exists():raise FileExistsError(a.out)
 r=run(a.certificate,a.nchi,a.dps);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2));print(json.dumps({k:v for k,v in r.items() if k not in ['constants','K6_components']},indent=2));print('ratio',r['constants']['improvement_vs_recorded'])
