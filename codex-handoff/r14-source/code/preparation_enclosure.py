"""NEW smooth-SPV-1 comparison enclosure; not recovered original R13 code.
All interval endpoints stored exactly as dyadic tuples; display strings are not proof inputs.
No numerical field modes are substituted into the certificate.
"""
from __future__ import annotations
import argparse, json, sys, time, math
from pathlib import Path
from fractions import Fraction
from mpmath import mp, iv
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'vendor'))
import active_enclosure as tj
tj.J=9
cst,add,sub,mul,div,inv,root,der,scale=tj.constant,tj.add,tj.sub,tj.mul,tj.div,tj.inv,tj.root,tj.der,tj.scale
zero=tj.zero

def rational(p,q=1):return iv.mpf(p)/q

def exact_fraction(tup):
    sign,man,exp,bc=tup
    return Fraction((-1 if sign else 1)*int(man))*Fraction(2)**int(exp)

def encode(a):return tj.serial(a)

def point_bound(a,upper=True):return iv.mpf(a.b if upper else a.a)

def tor(r):return r+2*iv.ln(1-r/2)

def jexp(a):
    b=[iv.exp(a[0])]
    for n in range(1,tj.J):b.append(sum((k*a[k]*b[n-k] for k in range(1,n+1)),zero())/n)
    return b

def scalar_switch(t):
    if t<=0:return iv.mpf(0)
    if t>=2:return iv.mpf(1)
    z=t/2
    # Real monotone switch; keep exponent <=0 for numerical conditioning.
    if t<=1:
        E=iv.exp(-1/z+1/(1-z));return E/(1+E)
    return 1-scalar_switch(2-t)

def endpoint_jets(right=False):
    delta=rational(1,128)
    eta=2*iv.exp(1/(1-rational(3,2)*delta)-2/(3*delta))
    vals=[iv.mpf([0,scalar_switch(2*delta).b])]
    for j in range(1,tj.J):
        B=point_bound(eta/delta**j)
        vals.append(iv.mpf([-B.b,B.b]))
    if right:
        return [1-vals[0]]+[(-1)**(j+1)*vals[j] for j in range(1,tj.J)]
    return vals

def switch_jet(lo,hi):
    if hi<=rational(1,64):return endpoint_jets()
    if lo>=2-rational(1,64) and hi<=2:return endpoint_jets(True)
    if lo>=2:return cst(1)
    t=cst(iv.mpf([lo.a,hi.b]));t[1]=iv.mpf(1)
    if hi<=1:
        z=scale(t,rational(1,2));arg=sub(inv(sub(cst(1),z)),inv(z))
        e=jexp(arg);return div(e,add(cst(1),e))
    z=scale(sub(cst(2),t),rational(1,2));e=jexp(sub(inv(sub(cst(1),z)),inv(z)))
    return sub(cst(1),div(e,add(cst(1),e)))

def root_enclosure(target):
    """Proposal via mp bisection, acceptance ONLY by outward iv inequalities."""
    def propose(y):
        lo=mp.mpf('0.5');hi=mp.mpf('1.99')
        for _ in range(190):
            mid=(lo+hi)/2
            if mid+2*mp.log(1-mid/2)>y:lo=mid
            else:hi=mid
        return (lo+hi)/2
    # decreasing tortoise: upper target -> lower radius
    pl=propose(mp.make_mpf(target._mpi_[1]));ph=propose(mp.make_mpf(target._mpi_[0]))
    den=2**100
    lo=rational(int(mp.floor(pl*den))-2,den)
    hi=rational(int(mp.ceil(ph*den))+2,den)
    ok=bool(tor(lo)>target.b and tor(hi)<target.a)
    if not ok:raise ArithmeticError('Radius proposal not certified')
    return iv.mpf([lo.a,hi.b])

def geometry_grid(time_div:int, integral_div:int):
    """Endpoint Darboux integrals. Reflection ensures integral[0,2]s=1 exactly."""
    assert time_div%64==0 and integral_div%time_div==0
    npre=2*time_div;end=197*time_div//64
    accum=iv.mpf(0);integrals={0:iv.mpf(0)}
    stride=integral_div//time_div
    for j in range(1,integral_div+1):
        accum+=scalar_switch(rational(j,integral_div))
        if j%stride==0:
            i=j//stride;t=rational(i,time_div)
            upper=accum/integral_div
            lower=(accum-scalar_switch(t))/integral_div
            integrals[i]=iv.mpf([lower.a,upper.b])
    for i in range(time_div+1,npre+1):
        # integral_0^t s = t-1 + integral_0^(2-t) s
        integrals[i]=rational(i,time_div)-1+integrals[npre-i]
    base=tor(rational(3,2))-1
    rows=[]
    for i in range(end+1):
        I=integrals[i] if i<=npre else rational(i,time_div)-1
        rr=root_enclosure(base+I)
        rows.append({'index':i,'time_fraction':[i,time_div],'switch_integral':encode(I),
                     'radius':encode(rr),'tortoise_inequalities_pass':True})
    return rows

def decode(a):
    aa,bb=a['exact_dyadic'];return iv.mpf([mp.make_mpf(tuple(aa)),mp.make_mpf(tuple(bb))])

def geometry_jet(s,r0):
    r=cst(r0)
    for j in range(tj.J-1):r[j+1]=mul(sub(cst(1),scale(inv(r),2)),s)[j]/(j+1)
    Ai=sub(scale(inv(r),2),cst(1))
    A=add(sub(cst(1),s),mul(s,Ai))
    B=div(A,mul(r,r));U=sub(mul(A,sub(cst(1),s)),div(der(r,2),r))
    return r,A,B,U

def comparison(B,U,x):
    D=add(cst(x),scale(B,1-x))
    # Exact bilinear range for D0 avoids repeated-x dependency; derivatives unchanged.
    candidates=[xx+(1-xx)*bb for xx in (x.a,x.b) for bb in (B[0].a,B[0].b)]
    D[0]=iv.mpf([min((z.a for z in candidates),key=lambda z:mp.make_mpf(z._mpi_[0])),
                  max((z.b for z in candidates),key=lambda z:mp.make_mpf(z._mpi_[1]))])
    w0=root(D)
    F0=add(scale(div(der(w0,2),w0),-rational(1,2)),scale(mul(div(der(w0),w0),div(der(w0),w0)),rational(3,4)))
    w1=div(add(U,F0),scale(w0,2))
    F1=add(scale(sub(div(der(w1,2),w0),div(mul(der(w0,2),w1),mul(w0,w0))),-rational(1,2)),
           scale(mul(div(der(w0),w0),sub(div(der(w1),w0),div(mul(der(w0),w1),mul(w0,w0)))),rational(3,2)))
    w2=div(sub(F1,mul(w1,w1)),scale(w0,2))
    return D,[w0,w1,w2]

def pmul(a,b):
    out=[iv.mpf(0)]*(len(a)+len(b)-1)
    for i in range(len(a)):
        for j in range(len(b)):out[i+j]+=a[i]*b[j]
    return out

def rate_enclosure(D,U,ws,epsmax):
    wv=[w[0] for w in ws];wp=[der(w)[0] for w in ws];wpp=[der(w,2)[0] for w in ws]
    w2=pmul(wv,wv);w4=pmul(w2,w2);duw=pmul([D[0],U[0]],w2)
    fw=pmul(wpp,wv);fd=pmul(wp,wp)
    coeff=[]
    for j in range(3,9):
        z=(duw[j] if j<len(duw) else 0)-w4[j]
        if j>=1 and j-1<len(fw):z+=-fw[j-1]/2+rational(3,4)*fd[j-1]
        coeff.append(z)
    eps=iv.mpf([0,epsmax.b])
    w=wv[0]+eps*wv[1]+eps**2*wv[2]
    if not bool(w>0):raise ArithmeticError('Width positivity failed')
    P=sum((abs(a)*eps**j for j,a in enumerate(coeff)),iv.mpf(0))
    rate=P/(2*w**3)
    return w,rate,coeff

def run(out:Path,time_div=64,directions=8,integral_div=8192,dps=40):
    start=time.perf_counter();mp.dps=dps+30;iv.dps=dps
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True)
    geom=geometry_grid(time_div,integral_div)
    (out/'geometry.json').write_text(json.dumps({'dps':dps,'time_div':time_div,'integral_div':integral_div,'rows':geom},indent=2))
    # Independent outward initial-state inequality. r_initial<7/4 -> D0 >=16/49.
    assert bool(decode(geom[0]['radius'])<rational(7,4))
    Dmin=rational(16,49);e0=rational(1,2048**2)
    initial_coef=1/(32*Dmin**3*(1-e0**3/(16*Dmin**3)))
    assert bool(initial_coef<2)
    sums=[iv.mpf(0) for _ in range(directions)];count=0;wmin=None
    f=(out/'cells.jsonl').open('w')
    for i in range(len(geom)-1):
        lo=rational(i,time_div);hi=rational(i+1,time_div)
        s=switch_jet(lo,hi)
        radius=iv.mpf([decode(geom[i+1]['radius']).a,decode(geom[i]['radius']).b])
        r,A,B,U=geometry_jet(s,radius)
        for j in range(directions):
            x=iv.mpf([rational(j,directions).a,rational(j+1,directions).b])
            D,ws=comparison(B,U,x);w,rate,cs=rate_enclosure(D,U,ws,e0)
            bound=point_bound(rate)
            sums[j]+=bound/time_div
            count+=1
            wmin=w.a if wmin is None or bool(w.a<wmin) else wmin
            f.write(json.dumps({'it':i,'ix':j,'w':encode(w),'rate_scaled':encode(rate),
              'coefficients_e3_to_e8':[encode(z) for z in cs]})+'\n')
        f.flush()
        if i%32==0:print('rows',i,'cells',count,'seconds',time.perf_counter()-start,flush=True)
    f.close()
    Cint=max(sums,key=lambda a:mp.make_mpf(a._mpi_[1]))+rational(2,2048)
    endpoint=exact_fraction(Cint._mpi_[1]);Cceil=(endpoint.numerator+endpoint.denominator-1)//endpoint.denominator
    log_width=rational(2*Cceil,2048**5)
    result={'run':'R93R14','D0_enclosure':'bilinear corner range v2','dps':dps,'time_div':time_div,'directions':directions,'integral_div':integral_div,
      'domain':{'time':['0','197/64'],'direction':['0','1'],'kappa_min':2048},'cells':count,
      'C_by_direction':[encode(c) for c in sums],'initial_coefficient_upper':encode(initial_coef),
      'initial_asinh_beta_bound':'2/kappa^6','C_interval':encode(Cint),'C_integer':Cceil,
      'width_relative_allowance_at_threshold':encode(iv.exp(log_width)-1),
      'positive_width_lower':encode(iv.mpf(wmin)),
      'checks':{'all_widths_positive':True,'initial_mismatch_retained':True,'r_initial_below_7_over_4':True,
                'all_radius_brackets_certified':True,'time_direction_coverage_by_construction':count==(197*time_div//64)*directions},
      'seconds':time.perf_counter()-start,'physical_scalar_accepted':False,
      'old_R13_C_is_input':False,'old_field_arrays_certified':False,'field_mode_evolutions':0}
    (out/'summary.json').write_text(json.dumps(result,indent=2))
    print('DONE',json.dumps({k:result[k] for k in ('C_integer','cells','seconds','checks')}),flush=True)
    return result

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,required=True);a.add_argument('--time-div',type=int,default=64);a.add_argument('--directions',type=int,default=8);a.add_argument('--integral-div',type=int,default=8192);a.add_argument('--dps',type=int,default=40);v=a.parse_args()
    run(v.out,v.time_div,v.directions,v.integral_div,v.dps)
