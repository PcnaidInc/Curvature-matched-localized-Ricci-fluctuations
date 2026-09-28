"""Causal null-beam response and independently normalized flat-vacuum benchmark.
This does NOT calculate the SPV-1 Schwarzschild quantum variance.
Units c=hbar=1 for numerical references; G factors removed from coefficient.
"""
from __future__ import annotations
import argparse, json, hashlib, time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp, quad
from scipy.special import erf, expit


def plain(x):
    if isinstance(x, dict): return {str(k): plain(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)): return [plain(v) for v in x]
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    return x


def rule(a: float,b: float,n: int):
    z,w=leggauss(n);return (a+b)/2+(b-a)*z/2,w*(b-a)/2


def bump(z):
    z=np.asarray(z,dtype=float);a=np.abs(z)
    ans=np.exp(-z*z);ans=np.where(a>=5,0.,ans)
    use=(a>4)&(a<5)
    u=a[use]-4
    ans[use]*=expit(1/u-1/(1-u))
    return ans


def piece_nodes(bounds,n):
    v=[rule(a,b,n) for a,b in zip(bounds[:-1],bounds[1:]) if b>a]
    return np.concatenate([i[0] for i in v]),np.concatenate([i[1] for i in v])


class Windows:
    def __init__(self,start=0.,end=12.,sigma=1.,n=64):
        if end-start <= 10*sigma: raise ValueError('start/end compact supports must be disjoint')
        self.start,self.end,self.sigma,self.n=start,end,sigma,n
        z,w=piece_nodes([-5,-4,0,4,5],n)
        self.z=z;self.qw=w*bump(z);self.qw/=self.qw.sum()
        self.ts=start+sigma*z;self.tf=end+sigma*z
        self.lo=start-5*sigma;self.hi=end+5*sigma
    def C(self,t):
        """CDF via independent adaptive quadrature, including known endpoints."""
        z=(t-self.start)/self.sigma
        if z<=-5:return 0.
        if z>=5:return 1.
        norm=quad(lambda u:float(bump(np.array([u]))[0]),-5,5,epsabs=1e-13,points=[-4,0,4])[0]
        pts=[u for u in [-4,0,4] if -5<u<z]
        return quad(lambda u:float(bump(np.array([u]))[0]),-5,z,epsabs=1e-13,points=pts)[0]/norm
    def kernel(self,t,radial=False,R=30.):
        """w=D0(t) C_start(t) integral_end g(f)(f-t)/D0(f) df."""
        if t<=self.lo or t>=self.hi:return 0.
        C=self.C(t)
        upper=self.end+5*self.sigma;lower=max(t,self.end-5*self.sigma)
        bounds=sorted(set([lower,upper]+[v for v in [self.end-4*self.sigma,self.end,self.end+4*self.sigma] if lower<v<upper]))
        q,w=piece_nodes(bounds,self.n)
        norm=self.sigma*sum(piece_nodes([-5,-4,0,4,5],self.n)[1]*bump(piece_nodes([-5,-4,0,4,5],self.n)[0]))
        g=bump((q-self.end)/self.sigma)/norm
        d0=(R-t) if radial else 1.
        df=(R-q) if radial else np.ones_like(q)
        if min(df)<=0 or d0<=0:raise ValueError('beam crosses caustic/background r=0')
        return float(d0*C*np.sum(w*g*(q-t)/df))
    def Wflat(self,a):
        """Exact compact-window Fourier identity, disjoint equal-shaped endpoints.
        W=[(1-i a T-exp(-i a T))*G - a G']/a^2, start=0.
        """
        a=np.atleast_1d(a).astype(float);u=a*self.sigma
        ph=np.exp(-1j*u[:,None]*self.z[None,:])
        G=ph@self.qw;Gp=ph@(-1j*self.sigma*self.z*self.qw)
        d=self.end-self.start;x=a*d
        A=-np.expm1(-1j*x)-1j*x
        small=np.abs(x)<1e-3
        if small.any():
            xs=x[small];A[small]=sum(-(-1j*xs)**j/float(np.math.factorial(j)) for j in range(2,9))
        out=np.empty_like(G)
        nonzero=a!=0
        out[nonzero]=(A[nonzero]*G[nonzero]-a[nonzero]*Gp[nonzero])/a[nonzero]**2
        out[~nonzero]=d*d/2+self.sigma**2*np.sum(self.qw*self.z**2)
        return out*np.exp(-1j*a*self.start)


def flat_variance(w:Windows,b=4.,n=240,amax=12.,direct=False):
    """Var(y)=G^2*hbar^2*c_xi*I; c_xi=xi^2-xi/3+1/30.
    Fourier convention W(a)=int w(t) exp(-iat)dt. k=(1,1,0,0).
    Longitudinal AND transverse Gaussian exponent widths equal b.
    """
    a,wa=rule(0.,amax/w.sigma,n)
    W=w.Wflat(a)
    if not direct:
        J=np.sqrt(2*np.pi)/b*erf(a*b/(2*np.sqrt(2)))
    else:
        j=[]
        for x in a:
            k,wk=rule(-x/2,12/b,192)
            j.append(np.sum(wk*np.exp(-b*b*k*k/2)*(-np.expm1(-b*b*(x*x+2*x*k)/2))))
        J=np.array(j)
    phase=2*np.pi/b**2*J
    I=float(np.sum(wa*a**4*np.abs(W)**2*phase)/(4*np.pi**3))
    return I


def source_pair_test():
    """Conserved on-shell annihilation-pair tensor, independent 4-vector check."""
    rng=np.random.default_rng(81247);g=np.diag([-1.,1.,1.,1.]);n=np.array([1.,1.,0.,0.])
    errors=[];nullvals=[]
    for j in range(80):
        p3=rng.normal(size=3);q3=rng.normal(size=3)
        p=np.r_[np.linalg.norm(p3),p3];q=np.r_[np.linalg.norm(q3),q3]
        pc=g@p;qc=g@q;P=p+q;Pc=g@P;dot=p@g@q;P2=P@g@P
        for xi in [0.,1/12,1/6,1/4]:
            T=-(np.outer(pc,qc)+np.outer(qc,pc))/2+g*dot/2+xi*(np.outer(Pc,Pc)-g*P2)
            A=float(n@T@n);ap=p[0]-p[1];aq=q[0]-q[1]
            scale=1+np.max(abs(T))+abs(dot)
            errors += [float(np.max(abs(P@T))/((1+np.max(abs(P)))*scale)),abs(np.trace(g@T)-(1-6*xi)*dot)/scale,abs(A-(-ap*aq+xi*(ap+aq)**2))/scale]
            if xi==1/6:nullvals.append(abs(A))
    return dict(max_normalized_identity_residual=max(errors),conformal_null_projection_nonzero=max(nullvals)>1e-3)


def quantum_checks():
    w=Windows(n=80);xs=np.array([0.,1/12,1/6,1/4]);poly=xs**2-xs/3+1/30
    z,qw=rule(-1,1,160)
    angular=np.array([np.sum(qw*(x-.25+.25*z*z)**2)/2 for x in xs])
    values=[]
    for b in [2.,4.,8.]:
        lo=flat_variance(w,b,160,10);mid=flat_variance(w,b,240,12);fine=flat_variance(w,b,360,16)
        d=flat_variance(w,b,240,12,True)
        values.append(dict(sigma=1.,delay=12.,spatial_width=b,I=fine,variance_coefficients=(poly*fine).tolist(),rms_coefficients=np.sqrt(poly*fine).tolist(),cutoff_quadrature_change=abs(mid/fine-1),coarse_change=abs(lo/fine-1),direct_phase_space_change=abs(d/mid-1)))
    # Exact polynomial covariance for null source at two different couplings.
    cross=np.array([[x*y-(x+y)/6+1/30 for y in xs] for x in xs])
    # Independent dimensional rescaling of ALL windows and spatial sizes.
    w2=Windows(start=0,end=24.,sigma=2.,n=80)
    sc=flat_variance(w2,8.,240,12)/flat_variance(w,4.,240,12)
    # Independent time-domain Fourier transform using quadrature and CDF kernel.
    tt,wt=piece_nodes([-5,-4,0,4,5,7,8,12,16,17],44)
    kval=np.array([w.kernel(t) for t in tt])
    freq=np.array([0.,.05,.2,.8,1.4,2.5,4.])
    direct=np.exp(-1j*freq[:,None]*tt)@(wt*kval)
    spectral=w.Wflat(freq)
    return dict(couplings=xs.tolist(),angular_coefficients=poly.tolist(),independent_angular_error=float(np.max(abs(poly-angular))),variance_runs=values,coupling_covariance_eigenvalues=np.linalg.eigvalsh(cross).tolist(),conformal_to_minimal_variance_ratio=float(poly[2]/poly[0]),conformal_to_minimal_rms_ratio=float(np.sqrt(poly[2]/poly[0])),length_doubling_variance_ratio=sc,expected_scaling=1/16,kernel_fourier_relative_error=float(np.max(abs(direct-spectral)/(1+abs(spectral)))))


def optical_checks():
    # Units M=1; compact optical start/end windows remain inside r>0.
    w=Windows(start=.02,end=.14,sigma=.006,n=60);R=.8
    # s in [-.01,.05], final in [.11,.17]; every early window precedes final.
    S=lambda t: np.exp(-((t-.08)/.028)**2)*(1+.6*np.sin(41*t))
    # integrate linear D''=-D0*S/2 from before all starts, then reset using homogeneous data.
    low=w.lo;high=w.hi
    sol=solve_ivp(lambda t,y:[y[1],-(R-t)*S(t)/2],(low,high),[0.,0.],method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
    ys=sol.sol(w.ts);yf=sol.sol(w.tf)
    response=(yf[0][None,:]-ys[0][:,None]-(w.tf[None,:]-w.ts[:,None])*ys[1][:,None])/(R-w.tf)[None,:]
    averaged=float(w.qw@response@w.qw)
    t,qt=piece_nodes([low,w.start-4*w.sigma,w.start,w.start+4*w.sigma,w.start+5*w.sigma,w.end-5*w.sigma,w.end-4*w.sigma,w.end,w.end+4*w.sigma,high],40)
    kernel=np.array([w.kernel(u,True,R) for u in t])
    swapped=float(-np.sum(qt*kernel*np.array([S(u) for u in t]))/2)
    # third route: one direct ODE per sampled start, no anti-derivative reuse
    ss,ww=rule(w.start-5*w.sigma,w.start+5*w.sigma,64);gs=bump((ss-w.start)/w.sigma);ww*=gs;ww/=ww.sum()
    direct=0.
    for s,weight in zip(ss,ww):
        ode=solve_ivp(lambda u,y:[y[1],-(R-u)*S(u)/2],(s,high),[0.,0.],method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
        direct+=weight*float(w.qw@(ode.sol(w.tf)[0]/(R-w.tf)))
    # Null/source point outside either end is exactly excluded by kernel.
    outside=max(abs(w.kernel(low-.01,True,R)),abs(w.kernel(high+.01,True,R)))
    # Changing linear integration constants does not change reset physical data.
    p,q=2.3,-.7
    shifted=((yf[0]+p+q*w.tf)[None,:]-(ys[0]+p+q*w.ts)[:,None]-(w.tf[None,:]-w.ts[:,None])*(ys[1]+q)[:,None])/(R-w.tf)[None,:]
    homogeneous=float(w.qw@shifted@w.qw)
    # Nonlinear optical ODE only, fixed prescribed Ricci, not nonlinear Einstein.
    a=.025;b=.145
    first=solve_ivp(lambda t,y:[y[1],-(R-t)*S(t)/2],(a,b),[0.,0.],rtol=1e-12,atol=1e-14,method='DOP853').y[0,-1]/(R-b)
    errs=[]
    for alpha in [.04,.02,.01,.005]:
        soln=solve_ivp(lambda t,y:[y[1],-alpha*S(t)*y[0]/2],(a,b),[R-a,-1.],rtol=1e-13,atol=1e-14,method='DOP853')
        effect=(soln.y[0,-1]-(R-b))/(R-b)
        errs.append(dict(alpha=alpha,nonlinear_optical_effect=float(effect),linear_effect=float(alpha*first),remainder=float(effect-alpha*first)))
    ratios=[abs(errs[i]['remainder']/errs[i+1]['remainder']) for i in range(len(errs)-1)]
    np.savetxt(OUT/'optical_kernel.csv',np.c_[t,kernel,[S(u) for u in t]],delimiter=',',header='affine_lambda,kernel,manufactured_Rkk_per_amplitude',comments='')
    return dict(M=1.,reference_r_intercept=R,start_center=w.start,end_center=w.end,width=w.sigma,full_support=[low,high],minimum_radius=R-high,quadrature_retarded_value=averaged,swapped_source_value=swapped,independent_initial_value_value=direct,retarded_vs_swapped_relative=abs(averaged/swapped-1),initial_value_vs_swapped_relative=abs(direct/swapped-1),outside_support_value=outside,homogeneous_data_cancellation_absolute=abs(homogeneous-averaged),nonlinear_optical_only=errs,quadratic_halving_ratios=ratios)


def main():
    t=time.perf_counter()
    result={'scope':'definition and reference milestone only; no SPV-1 interior quantum covariance generated','source_pair':source_pair_test(),'quantum_flat_reference':quantum_checks(),'schwarzschild_optical_reference':optical_checks()}
    q=result['quantum_flat_reference'];o=result['schwarzschild_optical_reference']
    gates={
      'flat_conserved_pair':result['source_pair']['max_normalized_identity_residual']<2e-13,
      'conformal_null_projection_present':result['source_pair']['conformal_null_projection_nonzero'],
      'angular_polynomial':q['independent_angular_error']<2e-13,
      'positive_coupling_covariance':min(q['coupling_covariance_eigenvalues'])>-1e-13,
      'flat_direct_phase_space':max(v['direct_phase_space_change'] for v in q['variance_runs'])<2e-7,
      'flat_cutoff_refinement':max(v['cutoff_quadrature_change'] for v in q['variance_runs'])<1e-6,
      'kernel_time_fourier':q['kernel_fourier_relative_error']<2e-8,
      'dimensional_rescaling':abs(q['length_doubling_variance_ratio']-1/16)<1e-12,
      'optical_retarded_integral':o['retarded_vs_swapped_relative']<2e-8,
      'optical_independent_ivp':o['initial_value_vs_swapped_relative']<2e-8,
      'optical_causal_support':o['outside_support_value']==0.,
      'optical_initial_data_accounting':o['homogeneous_data_cancellation_absolute']<2e-10,
      'optical_second_order_remainder':all(3.8<r<4.2 for r in o['quadratic_halving_ratios'])}
    result['gates']=gates;result['all_pass']=all(gates.values());result['seconds']=time.perf_counter()-t
    (OUT/'benchmark.json').write_text(json.dumps(plain(result),indent=2))
    print(json.dumps({'all_pass':result['all_pass'],'gates':gates,'seconds':result['seconds']},indent=2))
    return 0 if result['all_pass'] else 2

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args();OUT=args.out;OUT.mkdir(parents=True,exist_ok=True)
    # Python versions without np.math.
    import math
    np.math=math
    raise SystemExit(main())
