"""Read saved certificates/budgets and assemble the bounded release, no field evolution."""
from pathlib import Path
from fractions import Fraction as F
import json,sys
from mpmath import mp

def endpoint(a):
 s,m,e,b=a['exact_dyadic'][1];return (-1 if s else 1)*F(int(m))*F(2)**int(e)

def assemble(root):
 root=Path(root);load=lambda p:json.loads((root/p).read_text())
 campaigns={name:load('results/'+name+'/summary.json') for name in ('prep_coarse','prep_refined','prep_final55','prep_final70')}
 C=campaigns['prep_final70']['C_integer'];budget=load('results/budget_final.json');fine_repeat=load('results/budget_repeat.json')
 old=load('inputs/R13C_budget16384.json');up=endpoint(budget['constants']['integrated_component_allowance'])
 scale=10**9;bound_int=(up.numerator*scale+up.denominator-1)//up.denominator
 conservative=f'{bound_int//scale}.{bound_int%scale:09d}'
 checks0=load('verification/preparation_checks.json');checks1=load('verification/compatibility_checks.json')
 checks2=load('verification/final55_certificate_check.json');checks3=load('verification/final70_certificate_check.json')
 mp.dps=90
 ca=endpoint(campaigns['prep_final55']['C_interval']);cb=endpoint(campaigns['prep_final70']['C_interval'])
 diff_exact=abs(ca-cb)/max(abs(cb),F(1));diff=mp.mpf(diff_exact.numerator)/diff_exact.denominator
 ua_exact=endpoint(budget['constants']['integrated_component_allowance']);ub_exact=endpoint(fine_repeat['constants']['integrated_component_allowance'])
 ua=mp.mpf(ua_exact.numerator)/ua_exact.denominator
 bdiff_exact=abs(ua_exact-ub_exact)/abs(ub_exact);bdiff=mp.mpf(bdiff_exact.numerator)/bdiff_exact.denominator
 gates={'preparation_algebra_and_endpoint_checks':checks0['all_pass'],
    'domain_reference_compatibility_checks':checks1['all_pass'],
    'final55_saved_certificate_checks':all(checks2[k] for k in ('exact_direction_sums_and_ceiling_pass','positive_widths','gap_free_coverage')),
    'final70_saved_certificate_checks':all(checks3[k] for k in ('exact_direction_sums_and_ceiling_pass','positive_widths','gap_free_coverage')),
    'precision_runs_agree_on_integer_C':campaigns['prep_final55']['C_integer']==C,
    'precision_C_relative_change_below_1e_minus_30':bool(diff_exact<F(1,10**30)),
    'propagated_budget_repeat_below_1e_minus_30':bool(bdiff_exact<F(1,10**30)),
    'component_uses_new_C':budget['premise']['C']==C and not budget['premise']['old_R13_C_is_assumed']}
 result={'run':'R93R14-20261012','jira':'PCNBR-111','model':'specified smooth SPV-1 preparation and active Schwarzschild equation; not old spline arrays',
    'preparation_evidence_complete':all(gates.values()),'preparation_C':C,
    'scope':'all t in [0,197/64], all spectral direction x in [0,1], kappa>=2048',
    'relative_width_factor_allowance_at_threshold':campaigns['prep_final70']['width_relative_allowance_at_threshold'],
    'component_upper_outward_decimal':conservative,
    'new_component_allowance':budget['constants']['integrated_component_allowance'],
    'previous_R13C_component_upper':old['conditional_component_bound_upper'],
    'allowance_improvement_factor':str(mp.mpf(old['conditional_component_bound_upper'])/ua),
    'C_precision_relative_difference':str(diff),'component_precision_relative_difference':str(bdiff),
    'acceptance':gates,
    'campaign_summaries':campaigns,
    'propagated_budget':budget,
    'unbounded_errors':['comparison-width retarded response minus R12 local-order9 trial','historical numerical geometry and field integration/roundoff errors','finite-sum global error and remaining full physical acceptance requirements'],
    'physical_scalar_accepted':False,'physical_scalar_values':None,'full_error_bound':None,'feedback_verdict':None,
    'historical_R13_arrays_recovered':False,'historical_R13_C_assumed':False,'historical_field_arrays_certified':False,
    'actual_field_mode_evolutions':0,'new_finite_spectrum_sums':0,'new_full_interval_campaigns':len(campaigns),
    'total_full_interval_cells':sum(d['cells'] for d in campaigns.values()),
    'summed_campaign_seconds_not_turn_elapsed':sum(d['seconds'] for d in campaigns.values()),
    'software_limit':'analytic enclosure under documented interval arithmetic; not external formal verification'}
 return result
if __name__=='__main__':
 root=Path(__file__).resolve().parents[1];d=assemble(root)
 out=Path(sys.argv[1]);out.write_text(json.dumps(d,indent=2));print(json.dumps({k:d[k] for k in ['preparation_C','component_upper_outward_decimal','allowance_improvement_factor','acceptance']},indent=2))
