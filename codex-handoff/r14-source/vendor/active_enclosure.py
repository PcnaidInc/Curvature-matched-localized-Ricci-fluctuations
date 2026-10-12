"""New active-region comparison-width enclosures. NOT the absent R13 preparation proof.
Normalized time jets, mpmath interval arithmetic. Exact rational cell endpoints.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
from mpmath import mp, iv
J=8

def zero(): return iv.mpf(0)
def constant(x): return [iv.mpf(x)]+[zero() for _ in range(J-1)]
def add(a,b): return [x+y for x,y in zip(a,b)]
def neg(a): return [-x for x in a]
def sub(a,b): return add(a,neg(b))
def scale(a,c): return [x*c for x in a]
def mul(a,b): return [sum((a[k]*b[j-k] for k in range(j+1)),zero()) for j in range(J)]
def inv(a):
    b=[1/a[0]]
    for j in range(1,J): b.append(-sum((a[k]*b[j-k] for k in range(1,j+1)),zero())/a[0])
    return b

def div(a,b): return mul(a,inv(b))
def root(a):
    b=[iv.sqrt(a[0])]
    for j in range(1,J): b.append((a[j]-sum((b[k]*b[j-k] for k in range(1,j)),zero()))/(2*b[0]))
    return b

def der(a,n=1):
    for _ in range(n): a=[(j+1)*a[j+1] for j in range(J-1)]+[zero()]
    return a

def radius_jet(a,b):
    r=constant(iv.mpf([a,b]))
    for j in range(J-1): r[j+1]=((1 if j==0 else 0)-2*inv(r)[j])/(j+1)
    return r

def bounds(v): return mp.make_mpf(v._mpi_[0]), mp.make_mpf(v._mpi_[1])
def serial(v):
    a,b=bounds(v)
    return {'lower':mp.nstr(a,32),'upper':mp.nstr(b,32),'exact_dyadic':[list(x) for x in v._mpi_]}

def cell(rlo,rhi,xlo,xhi):
    r=radius_jet(rlo,rhi); ri=inv(r); one=constant(1)
    A=sub(scale(ri,2),one); B=mul(A,mul(ri,ri))
    x=iv.mpf([xlo,xhi]); D=add(one,scale(sub(B,one),1-x))
    U=neg(div(der(r,2),r)); w0=root(D)
    rat=div(der(w0),w0)
    F0=add(scale(div(der(w0,2),w0),-iv.mpf(1)/2),scale(mul(rat,rat),iv.mpf(3)/4))
    w1=div(add(U,F0),scale(w0,2))
    F1=add(scale(sub(div(der(w1,2),w0),div(mul(der(w0,2),w1),mul(w0,w0))),-iv.mpf(1)/2),
           scale(mul(rat,sub(div(der(w1),w0),div(mul(der(w0),w1),mul(w0,w0)))),iv.mpf(3)/2))
    w2=div(sub(F1,mul(w1,w1)),scale(w0,2))
    eps=iv.mpf([0,iv.mpf(1)/2048**2]); w=add(add(w0,scale(w1,eps)),scale(w2,eps**2))
    ub=scale(div(der(w),w),-iv.mpf(1)/2)
    H=add(mul(ri,ri),mul(A,ri))
    R2=mul(ri,ri)
    return {'w':w[0],'ubar':ub[0],'B':B[0],'Bprime':der(B)[0],
            'H':H[0],'Hprime':der(H)[0],'R2prime':der(R2)[0]}

def run(nr=32,nx=32,dps=45):
    mp.dps=dps+20;iv.dps=dps;start=time.perf_counter(); rows=[];lo={};hi={}
    for i in range(nr):
      a=iv.mpf(19)/25+iv.mpf(2)*i/(25*nr);b=iv.mpf(19)/25+iv.mpf(2)*(i+1)/(25*nr)
      for j in range(nx):
        x=iv.mpf(j)/nx;y=iv.mpf(j+1)/nx
        vals=cell(a.a,b.b,x.a,y.b)
        rows.append({'ir':i,'ix':j,'radius':serial(iv.mpf([a.a,b.b])),'direction':serial(iv.mpf([x.a,y.b])),**{k:serial(v) for k,v in vals.items()}})
        for k,v in vals.items():
          l,h=bounds(v);lo[k]=min(lo.get(k,l),l);hi[k]=max(hi.get(k,h),h)
    # Rounded *outward* rational promises. No fitted physical values.
    promises={'w_lower':'0.99','w_upper':'1.70','ubar_abs_upper':'2.10',
              'B_upper':'2.83','Bprime_abs_upper':'22.1','H_abs_upper':'3.89',
              'Hprime_abs_upper':'20','R2prime_abs_upper':'7.5'}
    checks={'positive_comparison_width':lo['w']>mp.mpf(promises['w_lower']),
      'comparison_width_upper':hi['w']<mp.mpf(promises['w_upper']),
      'comparison_rate_upper':max(abs(lo['ubar']),abs(hi['ubar']))<mp.mpf(promises['ubar_abs_upper']),
      'geometry_B':hi['B']<mp.mpf(promises['B_upper']),
      'geometry_Bprime':max(abs(lo['Bprime']),abs(hi['Bprime']))<mp.mpf(promises['Bprime_abs_upper']),
      'geometry_H':max(abs(lo['H']),abs(hi['H']))<mp.mpf(promises['H_abs_upper']),
      'geometry_Hprime':max(abs(lo['Hprime']),abs(hi['Hprime']))<mp.mpf(promises['Hprime_abs_upper']),
      'geometry_R2prime':max(abs(lo['R2prime']),abs(hi['R2prime']))<mp.mpf(promises['R2prime_abs_upper'])}
    return {'nr':nr,'nx':nx,'dps':dps,'cells':len(rows),'seconds':time.perf_counter()-start,
      'domain':{'r':['19/25','21/25'],'x':['0','1'],'epsilon':['0','1/2048^2'],'geometry':'smooth active Schwarzschild only'},
      'extrema':{k:[mp.nstr(lo[k],32),mp.nstr(hi[k],32)] for k in lo},
      'rational_promises':promises,'checks':{k:bool(v) for k,v in checks.items()},'rows':rows,
      'global_preparation_bound_reexecuted':False}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--nr',type=int,default=32);p.add_argument('--nx',type=int,default=32);p.add_argument('--dps',type=int,default=45);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 if a.out.exists():raise FileExistsError(a.out)
 d=run(a.nr,a.nx,a.dps);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(d,indent=2));print(json.dumps({k:v for k,v in d.items() if k!='rows'},indent=2))
