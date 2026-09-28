# Angular localization of two sampled leading Ricci projections

This calculation is independently implemented for a specified scalar state. It does
not claim novel physics, a final interior solution, or a Kerr-ring result.

## 1. Motivation and scope
The prior work averaged over a whole sphere. The current experiment asks how much
that averaging suppresses a well-defined local fluctuation. It does NOT insert a
patch into the spherical mass or Weyl formula. Instead it uses two general leading
Einstein relations: delta R=-8 pi G delta T and delta(R_ab u^a u^b)=8 pi G delta(pi^2)
for the chosen minimally coupled massless scalar. Here pi=u.grad phi. The first is
a scalar and the second uses the explicitly specified background interior observer.
Because the background Ricci tensor vanishes, their first perturbations have the
usual linear gauge invariance. Higher-order polarization, intrinsic gravitational
fluctuations, and the full Weyl/metric response are not determined by these identities.

We retain the original four-dimensional nonrotating Schwarzschild geometry and
SPV-1 smooth in-vacuum preparation. The auxiliary massive ultrastatic past and
mass-removal history define a state, not a realistic stellar-collapse solution.
We reuse that background/state code with attribution to the previous release, but
new phase-resolved modes generate the new angular covariance. Old whole-sphere
values are used only as regression data. No Jiang coefficients/tables enter it.

Units: geometric M=1, hbar omitted from mode amplitudes. Stress standard deviations
are in hbar/M^4; variances in hbar^2/M^8. The illustrative gravitational coupling is
epsilon=G hbar/M^2. Field is free, real, massless and MINIMALLY coupled on the final
patch, not a general matter sector. In particular, this fluctuating trace is not
the conformally coupled massless trace anomaly; an additive anomaly is a c-number
and does not fluctuate. A trace measurement can miss traceless fluctuations.

## 2. Measured operator and angular window
The operator is quadratic locally BEFORE averaging:

    Q_delta = integral d tau f_sigma(tau) integral dx g_b(x)
              integral dOmega K_delta(Omega) T(tau,x,Omega).

The time weight is normalized exp[-(tau/sigma)^2] times a C-infinity cutoff equal
one up to four exponent widths and zero at five. The longitudinal Gaussian is
normalized exp[-x^2/b^2]/(sqrt(pi)b). The proper longitudinal exponent width equals
.25M at the center; it evolves across the time support because the coordinate
window is fixed. The angular heat kernel, centered on an arbitrary axis, is

    K_delta(theta) = sum_J (2J+1)/(4pi)
                    exp[-delta^2 J(J+1)/2] P_J(cos theta).

It is normalized and positive. It has smooth tails across the sphere, not a hard
cap boundary. Delta is its diffusion-width parameter, not an exact angular radius.
The uniform limit is K=1/(4pi). Its effective solid angle is

    Omega_eff = 1 / integral K_delta^2 dOmega,
    Omega_eff/(4pi) = 1 / sum_J (2J+1) exp[-delta^2 J(J+1)].

Effective area is a weighting measure, not the fraction containing all support.
Primary location r/M=.6, sigma/M=.02. The exact full support is recorded in JSON.
The initial angular cases were delta=1,.5,.25 and uniform. Delta=.125 is an added
narrow-window check with separately extended external harmonics; its acceptance
uses the same frozen numerical tolerances.

## 3. State and field modes
The metric is ds^2=A(-deta^2+dx^2)+r^2 dOmega^2, with A=2/r-1 and r'=-A after
preparation. Expand phi in exp(ikx)/sqrt(2pi) times normalized Y_lm. For chi=r f,

    chi'' + [k^2 + A(l(l+1)/r^2 + m^2) - r''/r] chi=0.

The main solver evolves the Gaussian width and phase, starting in the ultrastatic
positive-frequency ground state. Direct complex solutions are independent checks.
Normalization is r^2(f d* - d f*)=i, with d=partial_eta f. The Gaussian evolution
builds in normalization, so the direct solver, not merely that identity, is needed
as an independent test. Every completed mode integration has a final state/hash.

## 4. Full angular coupling, not the l=0 field
For a trace pair with angular orders l,p and external multipole J, put
L=l(l+1), P=p(p+1), g=(L+P-J(J+1))/2. Integration by parts on the sphere gives
integral K gradY_l.gradY_p through this g times the corresponding Gaunt coefficient.
The two positive-longitudinal-momentum sign sectors have pair amplitudes

    B_J^s = integral d tau f_sigma [d_l d_p/A + s k k' f_l f_p/A
                                   - g f_l f_p/r^2].

There is NO complex conjugation in these vacuum aa pairs until taking their norms.
The Gaunt/Wigner orthogonality sums all magnetic indices exactly. The angular
power contributing to a unit-centered kernel is

    P_J = (2J+1)/(16pi^4) sum_(l,p) (2l+1)(2p+1) (3j(l,p,J;0,0,0))^2
          integral_0^infty dk dk' sum_s exp[-b^2(k+s k')^2/2] |B_J^s|^2.

The two-particle Fock norm and both signs of longitudinal momentum are included.
The implementation uses l<=p and doubles unequal pairs. Triangle/parity rules are
imposed. At J=0, l=p and this reduces to the old whole-sphere trace covariance.
At J>0, distinct field angular orders mix; no spherical gravitational formula is
assumed. For the time-direction projection, replace B by integral f_sigma d_l d_p/A.

The sampled trace variance is sum_J exp[-delta^2 J(J+1)] P_J. Cross-covariance
between two widths on the same axis uses exp[-(delta1^2+delta2^2)J(J+1)/2]. Its
positivity follows from the underlying pair norms and is checked numerically.
Both derivative field terms and angular gradient terms are retained.

## 5. Why some amplification is guaranteed, not a discovered mechanism
For a rotationally invariant state and fixed common time/longitudinal sampling,
the covariance is diagonal in angular multipoles with nonnegative powers. Every
normalized angular window has the same monopole contribution. The whole-sphere
average contains only that monopole; a nonuniform window adds nonnegative terms.
Thus it cannot have smaller variance than the sphere average in this setting.
The numerical question is HOW MUCH it grows and whether that exceeds the vacuum
control, not whether angular averaging suppresses any noise at all.

If the low-J spectrum is nearly white, P_J approximately (2J+1)P_0, the variance
increase is the inverse effective-area fraction. This is a reference approximation,
not an assumption that physical patches are statistically independent. The exact
calculation retains non-white spectra and all allowed couplings.

## 6. Independent references and gravitational normalization
The comparison is a frozen ultrastatic R^(1,1)xS^2 derivative-field vacuum at the
same sphere radius and local longitudinal scale, with the SAME time, longitudinal
and angular kernels. Its time integral follows from positive-frequency modes and
the Fourier transform of the compact Gaussian. A fast but equivalent evaluation
uses the analytic Gaussian transform minus the removed small tail. It is checked
against direct time quadrature; no fitting to interior results occurs.

This frozen geometry has different dynamics and state/history. Its Ricci tensor
is not the Schwarzschild tensor, and it is NOT being solved as an Einstein system.
The comparison is between stress projections, not a decomposition of independent
curved and flat random noises. Agreement cannot be converted into an exact
percentage of noise attributed to one cause.

An additional independent Minkowski point-worldline pure-time reference gives
Var(Q_trace)=2/(35 pi^4 sigma^8) and Var(Q_pi^2)=3/(70 pi^4 sigma^8). Direct double
positive-frequency integrals and the Fourier transform of the separated Wightman
contractions verify the normalization. These point-worldline controls have DIFFERENT
spatial sampling from the sphere experiment; they are not its numerical substitute.
The massless static sphere's undifferentiated zero mode has an infrared caveat;
only derivative observables with finite integrals are used here.

A nonzero classical comparison is the same-time sampled Schwarzschild sqrt(K):

    D = integral f_sigma sqrt(48)/r^3 d tau,
    RMS(Q_deltaR)/D = epsilon * 8pi sqrt(Var(Q_trace))/D.

The time-projection uses Var(Q_pi^2) similarly. We never divide by the zero
background Ricci scalar or by a nearly zero quantum mean. These are not variances
of K or percentages of all gravitational degrees of freedom. A standard deviation
is not a maximum event or a Gaussian-tail probability.

## 7. Validation and reusable histories
Cutoff, time/quadrature/tolerance, external angular-harmonic, spherical overlap,
Gaunt/gradient, coordinate rotation, independent mode, positivity and flat/static
reference checks are recorded separately. The rejected pilot is retained. Larger
runs use the same state and physical measurement, not retuned physics.

Final spectra are constructed from direct time integrals and explicit spectral
cutoffs. No fitted ultraviolet tail or pressure adjustment is inserted. Coupling
rescaling is algebraic, not a new nonlinear simulation. Planck-scale sampling
screens are declared conventions, not sufficient validity theorems.

To avoid repeatedly recomputing field phases for every future operator, validated
Chebyshev histories of u,W,phase are included for the final campaign. They are a
supplementary approximate cache, NOT newly evolved modes or replacements for the
accepted contraction data. Reconstruction tests and hashes are explicit. Full
reproduction regenerates transient dense histories; verification and assembly do
not. The cache is valid only within this measurement's time interval.

## 8. Relation to Abdul's hypothesis
This is a controlled test of an averaging loophole in the previous small spherical
results. It does not establish or refute quantum dominance near a singularity,
quantize the ring, predict a bounce, or test hit/bypass probabilities. It addresses
one free field and two leading Ricci projections in a finite nonrotating interior.
Localized Weyl/tidal response, full conserved tensor metric covariance, intrinsic
gravity, nonlinear effects and realistic formation/preparation remain separate.

Primary methodological sources: Phillips & Hu gr-qc/0010019 (noise kernels and
trace distinction); Hu & Verdaguer0802.0658 (induced vs intrinsic metric correlations);
NIST DLMF34.3.19–22 (Gaunt/Wigner sums); Wu, Ford & Schiappacasse2104.04446
(finite spacetime sampling). These are established frameworks, not our discoveries.
