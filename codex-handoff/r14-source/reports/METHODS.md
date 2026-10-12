# R93R14 — independently checkable smooth-state preparation enclosure

Execution date: 12 October 2026 UTC / 11 October America/Chicago. Research owner: Abdul Badran. New implementation and mathematical specialization are assistant-supplied. Ledger: PCNBR-111; feedback parent PCNBR-81. The current request is to continue the next batch and assess closeness to the final goal.

## 1. Result and strict scope

This release provides a NEW independently executable enclosure of the preparation comparison error for the explicitly defined **smooth SPV-1 equation**. It does not recover the missing original R13 source or arrays. It does not take the historically reported constant 20,547,271 as an input. A newly enclosed constant, saved in `results/prep_final70/summary.json`, is propagated through the available R13C amplitude-normalized scalar error-transfer formula.

The preparation-evidence dependency is therefore removed for this new state-width-replacement component, subject to the specified smooth equation, analytical derivation, implementation, and interval-library arithmetic. It is not a formally verified theorem about arbitrary software. In particular, it does NOT certify the historical cubic-spline geometry/field arrays, their quadrature, or the complete continuum scalar. The comparison-width retarded response is not the R12 local-order-nine trial. Their difference remains unbounded, as do the inherited numerical geometry/field/global errors.

The central question remains whether quantum physics becomes essential near the classical-singularity regime while relativity can remain useful macroscopically. This calculation is not a quantum-generated geometry, a feedback verdict, an inward scan, repulsion, a singularity endpoint, a Hawking-origin source, a Big Bang equivalence, or a ring-superposition/bypass test.

## 2. Fixed smooth target and unchanged measurement

Set T=2, M=1 and signature (-,+,+,+). For 0<t<2, z=t/2 and

    s(t) = 1/(1 + exp(1/z - 1/(1-z))),

extended flat by 0 and 1 outside the preparation interval. The radius obeys

    r' = -(2/r-1)*s,
    r(2)=3/2,
    A = 1-s+s*(2/r-1),
    m_squared = 1-s.

The canonical mode psi=r*f obeys

    psi'' + [k^2 + A*(L/r^2 + m_squared) - r''/r]*psi = 0,
    L = ell*(ell+1).

The initial state is normalized positive-frequency data with width
`sqrt(k^2+L/r_initial^2+1)` and zero amplitude rate. This is the specified prepared state, not an asserted astrophysical collapse vacuum.

For t>=2 the equation is the original massless Schwarzschild equation. The anisotropic pulse remains supported on r=0.84 to 0.76, with the original q*b=0 or 1, b=0.25*sqrt(3/7), and xi=0 or 1/6. The internal chi-weighted scalar and original SL0-1 optical detector remain distinct and unchanged. No angular or momentum sum is recalculated in this batch. Both sides of the current omitted spectrum satisfy the comparison threshold: ell>=2048 implies each leg's kappa>=sqrt(2048*2049)>2048; at lower ell, v>=6144 implies |v±q/2|>2048.

## 3. Radius enclosures without treating numerical roots as proof

Let tau(r)=r+2*log(1-r/2). Its derivative is -r/(2-r)<0 for 0<r<2. Since tau'(r)*r'=s,

    tau(r(t)) = tau(3/2)-1 + integral_0^t s(u) du.

The equality integral_0^2 s=1 follows from s(2-t)=1-s(t). For 0<=t<=1, monotone left/right Darboux sums enclose the integral. Reflection provides the rest of preparation exactly in terms of those enclosures; after t=2 the integral equals t-1.

High-precision bisection proposes a radius bracket. Its endpoints are then rounded outwards to rational multiples of 2^-100. Only two outward interval inequalities accept the bracket:

    tau(r_lower) > target_upper,
    tau(r_upper) < target_lower.

The independent saved-certificate checker re-evaluates those inequalities at 90 decimal digits and checks every exact time fraction. For each whole time cell the decreasing radius lies between the lower endpoint enclosure at the cell's right edge and the upper endpoint enclosure at its left edge. No dense output or numerical root convergence tolerance is substituted for this bound.

## 4. Switching endpoints and time jets

Normalized jets store f^(j)/j!, j=0..8. Products, inverse, square root and time differentiation follow coefficient recurrences. The available R13C jet primitives are copied unchanged in `vendor/active_enclosure.py`; their dimension is explicitly set to nine in the preparation process. New full-preparation logic generates the switch, radius, A, B and U jets. The comparisons below require derivatives only through the retained valid orders; padded upper jet coefficients are not used as certified high derivatives.

Interior switch derivatives are evaluated from the analytic exponential/logistic representation. Close to either flat endpoint, singular interval inversion is avoided. For real 0<z0<=delta=1/128, use the complex disk |z-z0|<=z0/2. On it,

    Re(1/z) >= 2/(3*z0),
    |1/(1-z)| <= 1/(1-3*delta/2).

The modulus of E=exp(-1/z+1/(1-z)) is less than 1/2. Hence |E/(1+E)|<=2|E|, and Cauchy's derivative estimate, including t=2z, gives

    |s^(j)(t)|/j! <=
       2*exp[1/(1-3*delta/2)-2/(3*delta)]*delta^(-j),
       j=1..8.

The function exp[-2/(3z)]/z^j increases on this short interval because delta<2/(3j). Reflection bounds the other endpoint. At t=0 and2, flatness supplies the limiting derivatives. The actual switch-value intervals use monotonicity. The analytical endpoint bounds are accompanied by independently differentiated point checks, not inferred from those samples.

## 5. Positive comparison width and factored oscillator defect

For kappa^2=k^2+L, x=k^2/kappa^2, epsilon=kappa^-2, define

    B=A/r^2,
    U=A*m_squared-r''/r,
    D=x+(1-x)*B,
    F(w)=-w''/(2w)+3*w'^2/(4*w^2).

The comparison is the same mathematical fourth-order polynomial used by R13C:

    w0=sqrt(D),
    w1=(U+F(w0))/(2*w0),
    F1=-(w1''/w0-w0''*w1/w0^2)/2
        +(3/2)*(w0'/w0)*(w1'/w0-w0'*w1/w0^2),
    w2=(F1-w1^2)/(2*w0),
    w=w0+epsilon*w1+epsilon^2*w2,
    Wbar=kappa*w.

Every whole cell includes epsilon in [0,1/2048^2], not a single high-frequency point. Width positivity is checked throughout. In the final implementation, D's zeroth coefficient is enclosed by the exact extrema of its bilinear function on the four x/B corners, avoiding an avoidable repeated-x dependency. Higher derivative coefficients remain (1-x)*B^(j)/j!. The old algebraically equivalent enclosure is preserved in the stage archive and source snapshot; this is not a physical model change.

For vbar=exp(-i integral Wbar)/sqrt(2 Wbar), its oscillator defect is

    d=omega^2-Wbar^2-Wbar''/(2Wbar)+3*Wbar'^2/(4Wbar^2).

Avoid evaluating a tiny difference of large quantities. Instead define

    N=(D+epsilon*U)*w^2-w^4
        +epsilon*(-w''*w/2+3*w'^2/4).

Independent SymPy algebra proves the epsilon^0, epsilon^1 and epsilon^2 coefficients vanish identically under the definitions above. The remaining polynomial is

    N=epsilon^3*sum_(j=0..5) N_(j+3)*epsilon^j.

Since d=kappa^2*N/w^2,

    |d|/(2Wbar) <= kappa^-5
       * sum_(j=0..5) |N_(j+3)|*epsilon^j/(2*w^3).

The code encloses the right-hand coefficient on every time/direction cell. It sums upper coefficients along each FIXED direction interval, then takes the maximum of those integrated bounds. This does not incorrectly allow the direction to vary with time. Exact dyadic interval endpoints, not display strings, determine the final conservative integer ceiling.

## 6. Initial data and mode-comparison theorem

Define instantaneous comparison coefficients by psi=alpha*vbar+beta*conj(vbar), with the variation-of-constants constraint alpha'*vbar+beta'*conj(vbar)=0. The normalized Wronskian gives |alpha|^2-|beta|^2=1. The off-diagonal mixing rate is |d|/(2Wbar), so

    asinh|beta(t)| <= asinh|beta(0)| + integral_0^t |d|/(2Wbar) du.

This is a comparison coefficient, not a claim about Hawking quanta. The source paper on Bogoliubov bounds supplies the general oscillator/bounding framework; the specific smooth geometry, enclosure and constants here are project calculations.

The initial radius enclosure proves r_initial<7/4, so D0>=16/49. At t=0 the derivatives of the switch vanish and

    Wactual/kappa=sqrt(D0+epsilon),
    Wbar/kappa=sqrt(D0)+epsilon/(2sqrt(D0))
               -epsilon^2/(8*D0^(3/2)).

Taylor's theorem bounds the positive square-root remainder by epsilon^3/(16*D0^(5/2)). Using -log(1-y)<=y/(1-y) gives

    asinh|beta(0)| = (1/2)*|log(Wactual/Wbar)| <= 2/kappa^6.

It is not set to zero. The final uniform constant includes 2/2048, since kappa>=2048. The script also stores a sharper initial coefficient and verifies it is below2, but the shared proof conservatively uses2.

Consequently, for the NEW integer C in the result,

    asinh|beta(t)| <= C/kappa^5,
    exp(-2C/kappa^5) <= Wactual/Wbar <= exp(2C/kappa^5),
    |zactual-zbar| <= Wbar*(exp(2C/kappa^5)-1),
    zactual=psi'/psi, zbar=-Wbar'/(2Wbar)-i*Wbar.

The last line follows from zactual-zbar=2i*Wbar*beta*conj(vbar)/psi and the Wronskian relation. The mode cannot vanish because its nonzero Wronskian is conserved. These are uniform statements on the specified smooth-target time/direction/high-frequency domain. They do not establish a rigorously accurate historical numerical trajectory.

## 7. Propagation through the existing measurement

This batch uses the available, unchanged R13C normalized error-transfer code, explicitly replacing its previously reported C input by the NEW enclosed C. The active-region input enclosure is included and hashed. This does not change the physical state or recalculate old finite field sums.

For M=sqrt(Wk*Wp), Y=M*Z obeys Y'=-i*Omega*Y+i*S/(4M). The same measurement weight requires g/M in the readout. For actual and comparison widths, the exact source, frequency, initial and changed-readout error terms are retained. R13C's coefficient K6, recomputed with the new C, yields the state-width-replacement envelope K6/(v^2+L)^3.

The nonnegative spectral measure is (2ell+1)/(4pi^2) dv. The omitted domain is partitioned without overlap: ell>=2048 at all v, and ell=0..2047 with v>=6144. The high-angular integral/monotone majorant and finite-angular alternating-series bracket are reused in the source-preserving transfer function. The resulting single-component allowance is saved with exact outward interval endpoints.

The new allowance is an upper bound on this defined smooth-state/comparison component under the analytical and software assumptions. It is not observed error, full stress, a gravitational percentage, an optical variance, or a feedback gain. The comparison-to-R12-local-trial and numerical geometry/field/global errors remain unavailable. The guard still rejects a completed physical scalar.

## 8. Verification and limitations

The completed campaigns and constants are in `results/FINAL_RESULTS.json`. Independent checks cover exact residual cancellations, the comparison-mode defect identity, nonzero initial mismatch, endpoint derivative implementation, rational convex-range samples, active-limit agreement, exact cell coverage, positive widths, exact summation/ceiling, certified radius endpoints and omitted-domain threshold compatibility. Two complete final evaluations use different arithmetic precision; they are numerical enclosure repeats, not independent physical experiments.

An initial auxiliary check substituted separately rounded interval midpoints into identities that require exact relationships, and failed an over-strong containment test. That check is preserved. It was corrected to evaluate the direct residual using 120-digit interval arithmetic and compare with the 80-digit factored interval. A nested high-precision derivative alternative timed out before saving results and was abandoned; no execution counts or successful tests are invented for it. These check-only problems did not enter or change the production enclosure calculations. A foreground transport timeout occurred after both final campaigns had launched; their live processes and saved logs were checked, so no duplicate campaigns were started.

mpmath's interval documentation describes rigorous enclosure semantics but also labels interval support experimental. This release uses supported arithmetic, square root, exponential and logarithm operations with outward endpoints and independent accounting checks. It is not external formal verification of mpmath, Python, the code, or the physics. Supporting software review and alternate implementations remain valuable. Earlier coarse bounds are valid conservative bounds, not numerical instabilities; finer cells reduce their looseness rather than relax any acceptance criterion.

## 9. Relation to the final goal

The advance is a checkable replacement for a previously inherited preparation premise, plus its explicit propagation into one error component. This moves the project closer to a controlled calculation. It does not supply the complete continuum scalar, the physical tensor/conservation calculation, equation-determined mean/initial geometry, shared-source optical covariance or the first feedback verdict. A separate comparison against classical gravity and a validity-controlled inward extension are still needed for the broader near-singularity claim. Task fractions and tighter uncertainty allowances are not percentages of proof.

## References

1. P. Boonserm and M. Visser, *Bounding the Bogoliubov coefficients*, arXiv:0801.0610v1, Annals of Physics 323 (2008). https://arxiv.org/html/0801.0610v1 — general oscillator and comparison-bound framework; not these state values.
2. mpmath documentation, *Contexts*, interval-arithmetic section. https://mpmath.readthedocs.io/en/latest/contexts.html — enclosure semantics and experimental-support caveat; installed version is recorded separately.
3. P. R. Anderson, C. Molina-Paris and E. Mottola, *Linear Response, Validity of Semi-Classical Gravity, and the Stability of Flat Space*, arXiv:gr-qc/0209075. https://arxiv.org/abs/gr-qc/0209075 — semiclassical response/validity framework; its explicit flat-space result is not the project interior result.
4. PCNBR-109/13816 and PCNBR-110/13923: historical R13 and delivered R13C specifications. Their original numbers are retained as prior reported/delivered results, not silently recovered arrays.

Public sources were inspected during this turn. No external paper, library test suite, manuscript approval or publication is invented. Jira binaries remain REFERENCE ONLY—NOT UPLOADED, and the old R4 catalog exception is unchanged.
