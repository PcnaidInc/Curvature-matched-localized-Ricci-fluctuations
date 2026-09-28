# Fixed-final-state curvature-coupling control — 2026-09-27

Freeze before new covariance results: xi = 0, 1/12, 1/6, 1/4; same scalar Cauchy
state, final Ricci-flat Schwarzschild geometry, observer and existing normalized
spacetime kernels. The conformal trace null is established physics, not a fitted
result; all nontrace results require separate numerical checks.

Use cached, hashed same-window Chebyshev mode histories of [internal record], and clearly
count new pair integrations separately from new mode evolutions (none planned in
production). Use all three final windows if their individual gates pass.

Primary gates:
- Independently constructed improved tensor: on-shell conservation, trace and
  trace-reversed observer identity; symbolic zero or representative residual <1e-8.
- Minimal covariance overlap <1e-6 relative for nonzero entries normalized by the
  diagonal geometric mean; new conformal trace implemented from the exact identity,
  with an independent unreduced derivation, never cancellation of noisy variances.
- Independent direct Hessian versus integration-by-parts amplitude <1e-6 relative
  on the dominant parent-amplitude scale. Document near-null cancellations.
- Field-angular extension (last 32 orders) <1% variance; last 20% explicit momentum
  contribution <1%; external-harmonic J24->32 <0.1% for primary localized windows.
- Time-quadrature refinement <0.1% and direct complex-mode history spot checks <1e-8.
- Covariance positivity >=-1e-12 relative to maximum diagonal. Near-zero eigenvalues
  and conformal trace are not treated with ill-conditioned relative errors.
- Direct conformal pair-sum versus transformed covariance tested with its own scale;
  preserve any failed settings and do not relax gates after results.
- Archive/manifest and resumable blocks verified; unchanged original manuscript/data.

No production of a coupling-dependent preparation, full backreacted metric, mean
renormalized stress, global collapse state, Kerr-ring physics or endpoint. Flat
plane-wave reference checks are not substituted for the interior. The frozen
non-Ricci-flat product used in the previous report is NOT an on-shell nonminimal
reference without changing its wave equation, and is not used to validate xi!=0.

## Refinement addition before final acceptance
Two fresh field/grid campaigns (outer_grid and inner_grid) and one cached-history
middle-time quadrature refinement were added after the initial implementation
checks to protect the cancellation-sensitive nontrace calculation. These are
additional independent checks, not relaxation of the gates. The outer/inner grids
raise the momentum range and quadrature/time sampling; the ODE tolerance stays
2e-11, matching the original fine runs. Their number of new field-mode labels is
recorded separately. Missing/in-progress comparisons were not accepted as passing.
