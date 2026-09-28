# Curvature-matched localized Ricci fluctuations: continuation after interruption

## Scientific question and scope
This is a new numerical completion of [internal record]/70, using the intact [internal record]
release. It asks whether local matter-induced curvature fluctuations increase
relative to a nonzero classical curvature scale when the temporal, longitudinal
and angular resolution conventions all follow the local curvature length inward.
It is not an absolute-vacuum recalculation, nonlinear geometry evolution, physical
collapse preparation, Kerr calculation, coherent-ring calculation, or novelty claim.

The original Schwarzschild field geometry is used, with geometric M=1 and
hbar omitted in numerical field amplitudes. The prepared scalar is free, real,
massless and minimally coupled on the final measurement region. Its ultrastatic
massive past and smooth transition are the inherited SPV-1 state definition, not
an assertion about how astrophysical black holes form. No Jiang numerical inputs
enter this work. The previous corrected mean geometry is NOT extended inward or
silently used as the background. A different field coupling or state may change
these results; in particular conformally coupled massless trace fluctuations are
not represented by the minimal-field trace studied here.

The previous interrupted [internal record] narrative records an unavailable intermediate
archive. These are newly executed, individually saved runs. The old intermediate
is not declared never to have existed and is not used as numerical evidence.

## Observable and units
In ds^2=A(-deta^2+dx^2)+r^2 dOmega^2, A=2/r-1, use u=partial_eta/sqrt(A).
The local quadratic operators are formed BEFORE averaging:

    T = (u.grad phi)^2 - |spatial grad phi|^2,
    P = (u.grad phi)^2,
    rho = P - T/2.

At leading Einstein order about a Ricci-flat reference:

    delta R = -8 pi G delta T,
    delta(R_ab u^a u^b) = 8 pi G delta P,
    delta(G_ab u^a u^b) = 8 pi G delta rho.

These relations do not require a spherical metric perturbation. They are NOT
sufficient to reconstruct localized Weyl curvature, a retarded metric, intrinsic
gravitational fluctuations, or higher-order polarization. Background Ricci and
Einstein tensors vanish; perturbations of these projected quantities are linearly
gauge invariant with the stated background frame. We never divide by zero Ricci.

A proper-time normalized Gaussian-shaped kernel exp[-tau^2/sigma^2] is smoothly
tapered from four to five exponent widths. The longitudinal kernel is a normalized
Gaussian in x; its proper exponent width is b at the CENTER. It evolves across
the time support. The angular kernel is the normalized positive heat kernel

    K_delta(theta) = sum_J (2J+1)/(4 pi) exp[-delta^2 J(J+1)/2] P_J(cos theta).

It has tails over the full sphere, not a sharp angular cap. The effective area
fraction is 1/sum_J (2J+1)exp[-delta^2 J(J+1)]. Its diffusion width delta is a
convention, not a standard deviation or literal patch radius. The full-sphere
kernel retains only J=0. Uniform sampling is a control, not an angularly matched
localized window. No early-to-late difference/measurement interval is introduced
in this single-time averaged observable.

The comparison denominator is D=<sqrt(K_Schwarzschild)>_time,
where sqrt(K)=sqrt(48)/r^3. Stress RMS units are hbar/M^4; variance units hbar^2/M^8.
The reported dimensionless induced-curvature RMS is

    ratio = epsilon * 8 pi * sqrt(Var(Q)) / D, epsilon=G hbar/M^2.

For energy Q=rho this is an Einstein, not Ricci, projection. The source covariance
matrix is ordered [T,P]. To form the Ricci covariance multiply by (8 pi G)^2 and
apply diag(-1,1) on both sides; the sign of the cross entry changes. The source
cross entry must not be mislabeled as a Ricci cross entry. Standard deviations
are not maximum excursions or Gaussian-tail probabilities.

## Matching every chosen center scale
Centers are r/M=.6,.45,.3. Let s=(r/.6)^(3/2). The proper-time exponent width is
.02 s, longitudinal width .25 s, and angular widths are delta=.25 sqrt(r/.6)
and .5 sqrt(r/.6). Thus sigma, b and r*delta are fixed fractions of K^(-1/4).
The full kernel support is checked using the exact proper-time primitive, not
merely the first/last quadrature nodes. It lies strictly above r/M=.2.
The underlying analytic field geometry is independently checked in this region.

Matching these center conventions does not make the whole spacetime geometry or
all directional expansion rates identical. It also does not prove full EFT
validity. The model still has smooth but globally extended angular/longitudinal
weights. Longitudinal homogeneity of the background makes that window definable;
it does not imply independent stochastic cells.

## Mode and pair computation
For chi=r f_lk, the inherited state evolves

    chi'' + [k^2 + A(l(l+1)/r^2 + m_aux^2) - r''/r] chi = 0.

Production uses the Gaussian width and phase, with a separately integrated complex
field equation as an independent check. Both start from the same specified ground
state in the ultrastatic past. Final states and provenance fingerprints are saved
for each angular order; transient time histories are represented by separately
checked interpolation caches in the final release.

Define time-integrated matrices for two field angular orders l,p:

    D_l p = integral f_time d_l d_p/A,
    K_l p = integral f_time f_l f_p/A,
    R_l p = integral f_time f_l f_p/r^2,
    B_J^s = D + s k k' K - [l(l+1)+p(p+1)-J(J+1)]R/2.

There is no complex conjugation inside the vacuum aa amplitude. Conjugation is
applied when forming its two-particle norm or cross product. The trace amplitude
is B, time projection is D, and energy is D-B/2. All three use the SAME mode pairs.
The longitudinal sign sectors carry weights exp[-b_coord^2(k+s k')^2/2]. Exact
magnetic-index sums yield a factor

    (2l+1)(2p+1)(2J+1) (3j(l,p,J;0,0,0))^2/(16 pi^4).

Only l<=p is explicitly evaluated; unequal pairs have multiplicity two. Gaunt
triangle/parity conditions are imposed. The stored six real moments recover the
full 2x2 covariance at each J, not just independent variances. Angular localization
then weights it by exp[-delta^2 J(J+1)]. Momentum and angular cutoffs are explicit;
no fitted ultraviolet tail or adjusted pressure enforces the result.

Connected covariances are unchanged by additive local c-number renormalizations
of the mean operator. That fact is NOT permission to use an unsmeared divergent
coincidence limit. The spacetime kernels, high-frequency convergence and field
normalization remain essential. Finite gravitational couplings and polarization
that change the gravitational response remain outside this leading-order map.

## Independent reference
A frozen ultrastatic product at the same radius and longitudinal scale is used as
a derivative-stress control with the same time and angular kernels. Its mode
frequency and time integrals are analytic positive-frequency expressions; no
coefficient is fitted to the interior data. It has a different state/history and
geometry, so its agreement is not an exact decomposition of noise into independent
flat and curved parts. Its gravitational field is not being solved. The massless
undifferentiated product-space zero-mode infrared issue is avoided only for the
finite derivative observables; no general phi^2 vacuum claim is made.

This continuation accelerates the compact-kernel Fourier evaluation using a
checked dense interpolation. Pointwise comparisons against the original direct
Fourier function and independent time quadrature are saved. The previous source
code, numerical normalization and reference formulas are preserved and credited.

## Validity and stop decision
The accepted domains and all numerical tolerances are in VALIDATION_PLAN.md and
results/FINAL_VALIDATION.json. Spectral comparisons are sensitivity tests, not
rigorous infinite-tail bounds or physical confidence intervals. The input state,
leading-order gravity and nonzero background normalization are separate assumptions.

The illustrative epsilon is 1e-9. Algebraic screens requiring the minimum of the
specified temporal, angular, longitudinal, radius and curvature scales to exceed
10/30/100 Planck units are also provided. These factors are chosen conventions,
not sufficient semiclassical-validity theorems. Rescaling to order-one curvature
noise outside that hierarchy is not a new simulation or a singularity solution.

The endpoint for this bounded project is a state/resolution-controlled result,
including a negative result. Reference agreement does NOT refute quantum importance:
ordinary quantum behavior could in principle become important. The issue here is
whether its computed gravitational magnitude is comparable while the calculation
remains within its declared regime. A stronger claim requires dynamics beyond
this fixed-background leading-order study, not merely a smaller sampling window.
