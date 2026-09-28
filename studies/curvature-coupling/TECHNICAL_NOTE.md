# Fixed-state curvature-coupling control

Date: 27 September 2026. Ledger: [internal record]; task: [internal record].
This is a new calculation supplement to the author-review manuscript, not a revised
claim of a singularity solution or a submitted publication. The original report
and numerical release are preserved unchanged.

## 1. Question and controlled variation
Does the reported localized Ricci-fluctuation conclusion depend on the minimally
coupled scalar stress tensor? We vary the dimensionless curvature coupling xi,
not Newton's constant and not the black-hole mass. Values 0, 1/12, 1/6, 1/4 were
specified before the new covariance calculation. The final Schwarzschild geometry,
Cauchy two-point function, observer and measurement kernels are identical.

We use signature (-,+,+,+), R_ab=0 and the previously recorded curvature convention.
With action -1/2 integral sqrt(-g)[(grad phi)^2+xi R phi^2], the massless field
obeys (box-xi R)phi=0. Consequently the propagation on this FINAL Ricci-flat patch
is identical for every xi. The same state can be used on the common field algebra.
This does NOT mean that changing xi throughout the auxiliary state preparation
would create the same state: that earlier geometry is not generally Ricci flat.
We do not perform that second experiment here.

## 2. Improved operator, trace and gravity
The full stress is T_ab^(xi)=T_ab^(0)+xi[G_ab phi^2+(g_ab box-nabla_a nabla_b)phi^2].
On the reference, G_ab=0. Set B=T^(0)=-(grad phi)^2, P=(u.grad phi)^2, and
H=u^a u^b nabla_a nabla_b(phi^2). The on-shell identity box(phi^2)=-2B gives

    T^(xi)=(1-6xi) B,
    Q_xi=(T_ab^(xi)-g_ab T^(xi)/2)u^a u^b=P-xi(H+B).

The observer u=A^(-1/2) partial_eta is geodesic on the final homogeneous metric
A(-deta^2+dx^2)+r^2 dOmega^2; hence H=d^2(phi^2)/d tau^2. Define
E=Q_(1/6)=P-(H+B)/6. Then Q_xi=(1-6xi)P+6xi E. The connected covariance of
(T^(xi),Q_xi) is an exact linear transformation of the computed covariance of
(B,P,E), with no independent-noise assumption.

At leading Einstein order only,

    delta R = -8 pi G delta T^(xi),
    delta(R_ab u^a u^b)=8 pi G delta Q_xi.

These algebraic projections do not solve a retarded metric, the Weyl tensor,
intrinsic gravitational noise, or higher-order polarization. The reference Ricci
tensor is zero, so all displayed gravity ratios use the NONZERO sampled
Schwarzschild sqrt(K), never zero background Ricci. The c-number anomaly can give
a nonzero conformal mean trace; it does not generate a connected trace fluctuation
on the fixed background. The known conformal trace null is a benchmark, not a new
physical discovery or a conclusion that all stress fluctuations vanish.

## 3. Sampling and state
Use the inherited SPV-1 prepared scalar state. Its artificial massive ultrastatic
past is a definition of Cauchy data, not a model of astrophysical collapse. At
r_c/M=0.6,0.45,0.3, let s=(r_c/0.6)^(3/2). Proper-time Gaussian exponent width
sigma/M=0.02 s; longitudinal proper width at the center b/M=0.25 s; angular
heat-kernel diffusion widths delta=0.5 sqrt(r_c/0.6) and 0.25 sqrt(r_c/0.6).
The full-sphere control has J=0 only. Kernels are formed AFTER the quadratic
operator, not by squaring a smeared linear field. Time weighting has smooth
compact taper from four to five exponent widths; angular/longitudinal weights
have global smooth tails. Center proper angular length is r_c delta. Their scale
conventions match the earlier report; this is not a new width sweep.

Stress units are hbar/M^4. If V is a source variance, the displayed dimensionless
curvature RMS is epsilon*8pi*sqrt(V)/D, D=integral g(tau)sqrt(K)dtau and
epsilon=G hbar/M^2. The table example uses epsilon=1e-9, as the manuscript did.
It is an illustrative choice, not an observational estimate. A small value is not
proof of semiclassical validity, and no Planck-limit extrapolation is performed.

## 4. New computation versus reused histories
The three primary calculations use the previously generated, separately checked
Chebyshev histories for f_lk and its eta derivative on the exact SAME intervals.
Every reused history file is hash checked. The histories are a numerical input,
not newly evolved modes and not arbitrary-interval propagators. The improved
operator requires NEW pair/time integrals and cross-covariances absent from the
old two-operator tables. Two independent momentum/time-grid refinements evolve
fresh phase-resolved modes from the same initial state; those runs are counted
separately. A middle-window time refinement uses the saved histories.

For each l,p, time integrals are D=<d_l d_p/A>, K=<f_l f_p/A>,
R=<f_l f_p/r^2>, H=<d_tau^2(f_l f_p)>.
Because the normalized coordinate-space angular/longitudinal windows are time
independent, integration by parts gives H=integral g''(tau) f_l f_p dtau, with no
boundary contribution. Analytic derivatives include the smooth compact taper;
we do not differentiate noisy covariance data or omit the taper.

With sign s_k=+/- for longitudinal momentum, and a=J(J+1)/2,

    B_J=D+s_k k k' K-[l(l+1)+p(p+1)]R/2+aR,
    E_J=D-(H+B_J)/6.

There is no complex conjugation inside these annihilation-pair amplitudes.
Conjugation appears in their Gram matrix. Positive pair norms, exact magnetic
Gaunt sums, parity/triangle restrictions, and the same 1/(16pi^4) normalization
as the checked baseline produce the three-source covariance. For numerical
conditioning E is formed at the amplitude level, not by subtracting large
completed variances. Stored per-pair four-moment matrices recover its J dependence;
an independent direct-E sum tests that reduction on E's own variance scale.

## 5. Independent checks
The full 4D improved tensor is constructed independently for the exact on-shell
field phi=x+(r-1)cos(theta). All four conservation components, the trace law and
observer projection reduce symbolically to zero; Ricci flatness is checked from
the connection. The general improvement divergence is -R_ba grad^a(phi^2), zero
here. This analytic identity is stronger than merely enforcing a numerical check.
Direct complex-mode integrations test the saved histories at off-fit points.
Direct mode Hessians and the sampling-derivative route are compared independently.

A separate Minkowski time-only Gaussian benchmark yields the covariance of (B,P,E)
as 1/(pi^4 sigma^8) times

    [[2/35, 3/70, 0],
     [3/70, 3/70, 1/420],
     [0,    1/420, 1/420]].

This follows both from exact polynomial/angular integrals and an independent
positive-frequency double quadrature. It is not the spatially smeared interior
observable and is not used to replace the interior result. In particular, the
old frozen spherical product has nonzero Ricci scalar and cannot be reused as an
on-shell nonminimal control without changing its field equation and state. No
such incorrectly frozen reference is included.

Acceptance is specified in VALIDATION_PLAN.md. Cutoff tails and quadrature changes
are numerical sensitivities, not rigorous infinite-tail or total physical error
bounds. Null-trace rows are checked algebraically, not divided by zero. All failed
or incomplete-development gates are retained separately from final acceptance.

## 6. Interpretation
The trace diagnostic is field-coupling dependent: its standard deviation changes
by exactly |1-6xi| and vanishes at the conformal value. This mathematical result is
known and does not itself test a new black-hole mechanism. The new numerical
nontrace result is needed to decide whether a zero trace hides a nonzero selected
gravitational source. Nonzero nontrace covariance does not specify every stress
component or the full metric. Suppression under one change of field theory is
not dynamic stabilization within a fixed theory.

This control constrains how broadly the manuscript can be interpreted. It does
not refute Abdul's quantum-importance question; it rules out using one
minimal-field trace curve as a universal diagnostic for all quantum matter.
The next distinct objective is a defined causal, gauge-aware nontrace gravitational
response with conserved source correlations, not another trace-only window sweep.
A coupling-dependent preparation and a realistic collapse state remain separate
extensions. Neither a Kerr ring nor its coherence/transport is represented.

## References actually consulted
1. N. G. Phillips and B. L. Hu, Phys. Rev. D 63, 104001 (2001),
   DOI 10.1103/PhysRevD.63.104001, arXiv:gr-qc/0010019v2, particularly Sec. III.3.
2. P. J. Brown, C. J. Fewster and E.-A. Kontou, Gen. Relativ. Gravit. 50, 121 (2018),
   DOI 10.1007/s10714-018-2446-5, Secs. 2–3. Their (+---) convention is not copied
   blindly; the Ricci-flat tensor identities above were independently derived.
3. B. L. Hu and E. Verdaguer, Living Rev. Relativity 11, 3 (2008),
   DOI 10.12942/lrr-2008-3, arXiv:0802.0658. Metadata/abstract consulted for scope.

## 7. Final accepted numerical summary
All 166 recorded numerical/identity gates pass across 36 joint cases. This count includes separate coupling/window checks, not independent physical experiments. New covariance blocks: 960; fresh production field labels: 90112; prior phase histories reused and checked: 480. Summed block wall time 2007.045895 s is not turn elapsed time.
Conformal observer-source RMS/minimal RMS in the narrower windows at r/M=.6,.45,.3: 0.0172379426522, 0.0162270012124, 0.0159617053631.
Largest dimensionless covariance comparisons: {"baseline_overlap": 1.532660428155706e-13, "angular_covariance": 6.163186921135034e-08, "momentum_tail_covariance": 5.572696761839128e-05, "external_harmonic_covariance": 2.8880818545946805e-08, "independent_grid_covariance": 2.3236094649478753e-07}. Multiply these comparisons by100 for percentages. These normalizations are not total error bars.
The independent unreduced conformal trace check has maximum residual 1.3181e-15 on its parent amplitude scale. The flat exact/numerical reference agrees within4.4650e-15, and representative direct-Hessian versus integration-by-parts comparisons within4.8743e-10 on the parent scale. The latter is not a bound on every mode or cancellation condition.
A fresh-directory assembly exactly reproduces COUPLING_COMPARISON.csv; resume skips128 saved outer blocks with no new production or pair evaluations and unchanged hashes; a deliberately corrupted disposable block is rejected. Initial/reference checks and each fresh validation pass evolve six small complex reference modes, counted separately from the two production refinements.
