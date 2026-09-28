"""Independent 4D tensor checks, not a quantum-mode simulation."""
from pathlib import Path
import argparse,json,time
import sympy as s

def run():
    began=time.perf_counter()
    r,x,th,ph,xi=s.symbols('r x theta phi xi', real=True)
    coords=[r,x,th,ph];h=2/r-1
    g=s.diag(-1/h,h,r**2,r**2*s.sin(th)**2);ig=g.inv()
    simp=lambda a:s.factor(s.trigsimp(s.cancel(a)))
    G=[[[simp(sum(ig[a,d]*(s.diff(g[d,c],coords[b])+s.diff(g[d,b],coords[c])-s.diff(g[b,c],coords[d])) for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]
    Ric=s.Matrix(4,4,lambda a,b:simp(sum(s.diff(G[c][a][b],coords[c])-s.diff(G[c][a][c],coords[b])+sum(G[c][c][d]*G[d][a][b]-G[c][b][d]*G[d][a][c] for d in range(4)) for c in range(4))))
    k=s.Matrix([-1,1/h,0,0]);kc=g*k
    accel=s.Matrix([simp(sum(k[b]*s.diff(k[a],coords[b])+sum(G[a][b][c]*k[b]*k[c] for c in range(4)) for b in range(4))) for a in range(4)])
    expansion=simp(sum(s.diff(k[a],coords[a])+sum(G[a][a][b]*k[b] for b in range(4)) for a in range(4)))
    # screen B_{AA}=e_A^a e_A^b nabla_b k_a on the reference spheres
    Btheta=simp(-sum(G[c][2][2]*kc[c] for c in range(4))/r**2)
    Bphi=simp(-sum(G[c][3][3]*kc[c] for c in range(4))/(r**2*s.sin(th)**2))
    # Nontrivial exact harmonic solution of Box phi=0; no assumption based on a fitted profile.
    f=x+(r-1)*s.cos(th)
    df=s.Matrix([s.diff(f,c) for c in coords]);F=f**2
    Hess=s.Matrix(4,4,lambda a,b:simp(s.diff(F,coords[a],coords[b])-sum(G[c][a][b]*s.diff(F,coords[c]) for c in range(4))))
    boxF=simp(s.trace(ig*Hess))
    boxphi=simp(sum(ig[a,b]*(s.diff(f,coords[a],coords[b])-sum(G[c][a][b]*df[c] for c in range(4))) for a in range(4) for b in range(4)))
    grad2=simp((df.T*ig*df)[0])
    T=s.simplify(df*df.T-g*grad2/2+xi*(g*boxF-Hess))
    mixed=ig*T
    div=[]
    for b in range(4):
        v=sum(s.diff(mixed[a,b],coords[a])+sum(G[a][a][c]*mixed[c,b]-G[c][a][b]*mixed[a,c] for c in range(4)) for a in range(4))
        div.append(simp(v))
    B=-grad2
    derivative=simp((k.T*df)[0]);second=simp(sum(k[a]*s.diff(derivative,coords[a]) for a in range(4)))
    null=simp((k.T*T*k)[0]);expected=simp((1-2*xi)*derivative**2-2*xi*f*second)
    trace=simp(s.trace(mixed))
    gates={
      'ricci_flat_reference':all(q==0 for q in Ric),
      'radial_tangent_null':simp((k.T*g*k)[0])==0,
      'radial_tangent_affine':all(q==0 for q in accel),
      'expansion_and_zero_shear':simp(expansion+2/r)==0 and simp(Btheta+1/r)==0 and simp(Bphi+1/r)==0,
      'background_raychaudhuri':simp(-s.diff(expansion,r)+expansion**2/2)==0,
      'exact_scalar_field_equation':boxphi==0,
      'all_improved_stress_conservation':all(q==0 for q in div),
      'trace_factor':simp(trace-(1-6*xi)*B)==0,
      'conformal_trace_zero':simp(trace.subs(xi,s.Rational(1,6)))==0,
      'null_hessian_vs_along_geodesic':simp(null-expected)==0,
      'conformal_null_source_not_identically_zero':simp(null.subs(xi,s.Rational(1,6)))!=0,
    }
    return {'M':1,'signature':'-+++','field':'x+(r-1)cos(theta)',
      'k_components':[str(q) for q in k], 'expansion':str(expansion),
      'screen_eigenvalues':[str(Btheta),str(Bphi)],
      'field_first_affine_derivative':str(derivative),'field_second_affine_derivative':str(second),
      'gates':gates,'all_pass':all(gates.values()),'seconds':time.perf_counter()-began,
      'scope':'exact background optical and conserved-source verification; not SPV-1 covariance or a nonlinear Einstein solve'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    data=run();(a.out/'curved_optics.json').write_text(json.dumps(data,indent=2));print(json.dumps(data,indent=2))
    raise SystemExit(0 if data['all_pass'] else 2)
