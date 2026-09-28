"""Independent four-dimensional improved stress construction; no production covariance input."""
import sympy as s
from pathlib import Path
import json,time
ROOT=Path(__file__).resolve().parents[1]
def main():
 t0=time.perf_counter()
 r,x,th,ph,xi=s.symbols('r x th ph xi', real=True);c=[r,x,th,ph];h=2/r-1
 g=s.diag(-1/h,h,r*r,r*r*s.sin(th)**2);inv=g.inv()
 def simp(a):return s.factor(s.trigsimp(s.cancel(a)))
 G=[[[simp(sum(inv[a,d]*(s.diff(g[d,b],c[e])+s.diff(g[d,e],c[b])-s.diff(g[b,e],c[d])) for d in range(4))/2) for e in range(4)] for b in range(4)]for a in range(4)]
 Ric=s.Matrix(4,4,lambda a,b:simp(sum(s.diff(G[d][a][b],c[d])-s.diff(G[d][a][d],c[b])+sum(G[d][d][e]*G[e][a][b]-G[d][b][e]*G[e][a][d] for e in range(4))for d in range(4))))
 phi=x+(r-1)*s.cos(th)
 grad=s.Matrix([s.diff(phi,v)for v in c]);X=simp((grad.T*inv*grad)[0]);q=phi**2
 def Hess(q):return s.Matrix(4,4,lambda a,b:simp(s.diff(q,c[a],c[b])-sum(G[d][a][b]*s.diff(q,c[d])for d in range(4))))
 hf=Hess(phi);hq=Hess(q);boxphi=simp(s.trace(inv*hf));boxq=simp(s.trace(inv*hq))
 T0=grad*grad.T-g*X/2;imp=g*boxq-hq
 T=s.simplify(T0+xi*imp);trace=simp(s.trace(inv*T));Tu=inv*T
 div=[]
 for b in range(4):
  div.append(simp(sum(s.diff(Tu[a,b],c[a])+sum(G[a][a][d]*Tu[d,b]-G[d][a][b]*Tu[a,d]for d in range(4))for a in range(4))))
 # u = -sqrt(h) d/dr (orientation irrelevant in quadratic/Hessian contractions)
 S=T-g*trace/2;P=h*grad[0]**2;B=-X;H=h*hq[0,0]
 projected=simp(h*S[0,0]-(P-xi*(H+B)))
 out=dict(ricci_flat=Ric==s.zeros(4),KG=str(boxphi),box_phi_squared_minus_2gradient=str(simp(boxq-2*X)),
    trace_law=str(simp(trace-(1-6*xi)*B)),conformal_trace=str(simp(trace.subs(xi,s.Rational(1,6)))),
    divergence=[str(z)for z in div],observer_projection=str(projected),seconds=time.perf_counter()-t0,
    test_field='phi=x+(r-1)cos(theta); exact harmonic test with angular and longitudinal dependence',
    general_identity='div(g box - Hess) q = -Ricci . grad(q); zero on Ricci-flat background for any smooth q',
    all_pass=(Ric==s.zeros(4) and boxphi==0 and all(z==0 for z in div) and projected==0 and simp(trace-(1-6*xi)*B)==0))
 (ROOT/'verification/TENSOR_CHECK.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
 if not out['all_pass']:raise SystemExit(1)
if __name__=='__main__':main()
