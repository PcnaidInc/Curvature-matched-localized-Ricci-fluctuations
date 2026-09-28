"""Independent Gaussian time-only null-stress reference.
This auxiliary normalization check does not have the beam benchmark's spatial window.
"""
import json,argparse
from pathlib import Path
import numpy as np
import sympy as s
from scipy.integrate import quad
from numpy.polynomial.legendre import leggauss

def run():
    u=s.symbols('u',real=True)
    # Two-particle phase-space reduction at fixed total frequency, no spatial smearing.
    spatial_poly=s.integrate((1-u*u)*(1-u)**4,(u,-1,1))
    # Independent direct Wick contraction of k.grad phi at a temporal separation:
    # <k.grad phi k.grad phi'>=2/(pi^2(dt-i0)^4), connected square=8/(pi^4(dt-i0)^8).
    wick_rational=s.Rational(8)*s.Rational(48)/s.factorial(7)
    spectral_rational=s.Rational(1,30)*spatial_poly*s.Rational(48,64)
    # independent p,q positive-frequency integrations for the same Wick correlator.
    # angular average <(1-cos(theta))^2>=4/3, and radial p^3 q^3.
    integrand=lambda q,p:p**3*q**3*np.exp(-(p+q)**2/2)
    direct=quad(lambda p:quad(lambda q:integrand(q,p),0,14,epsabs=1e-12,epsrel=1e-12)[0],0,14,epsabs=1e-12,epsrel=1e-12)[0]
    result=2/(np.pi**4*9)*direct
    exact=8/(105*np.pi**4)
    gates={'spatial_polynomial':spatial_poly==s.Rational(64,21),
           'independent_wick_spectral_prefactor':wick_rational==spectral_rational==s.Rational(8,105),
           'direct_positive_frequencies':abs(result/exact-1)<2e-10}
    return {'gates':gates,'all_pass':all(gates.values()),'spatial_polynomial':str(spatial_poly),
      'minimal_time_only_variance_coefficient':exact,'direct_quadrature_value':result,
      'relative_difference':abs(result/exact-1),
      'formula':'Var(T_kk time Gaussian)=8 hbar^2/(105 pi^4 sigma^8) at xi=0',
      'scope':'Independent normalization only: timelike sampling at fixed position with null tensor projection. Not null-line sampling, not the beam window, not a black-hole result.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    data=run();(a.out/'normalization.json').write_text(json.dumps(data,indent=2));print(json.dumps(data,indent=2));raise SystemExit(0 if data['all_pass'] else 2)
