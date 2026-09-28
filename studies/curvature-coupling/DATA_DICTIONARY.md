# Data dictionary

- `center`: areal radius divided by geometric black-hole mass M.
- `delta_rad`: heat-kernel angular diffusion-width convention. Zero denotes the
  uniform-sphere CONTROL; it does not mean a pointlike angular measurement.
- `xi`: dimensionless curvature coupling. It is not epsilon or a black-hole mass.
- `trace_rms`, `observer_source_rms`: connected source standard deviations,
  numerical units hbar/M^4, not maxima or energy changes.
- `trace_rms_vs_minimal`, `observer_rms_vs_minimal`: ratios at IDENTICAL state and
  window to xi=0, not percentages of all gravity. Trace ratio is |1-6xi| exactly.
- `source_cross_covariance`: symmetrized connected covariance of (T_xi,Q_xi),
  units hbar²/M^8. The Ricci-scalar/time-projection cross entry changes sign because
  delta R=-8piG T_xi while delta R_uu=+8piG Q_xi.
- `*_rms_per_epsilon`: 8pi times source RMS divided by the NONZERO sampled sqrt(K).
- `*_percent_at_epsilon1e9`: percentage for epsilon=10^-9 (the name's '1e9' is the
  column identifier; the exponent is NEGATIVE). Formula =100*10^-9*coefficient.
- `mean_background_sqrtK`: properly time-averaged sqrt(48)/r³, units M^-2.
- `minimal_variance_T`, `minimal_variance_P`, `conformal_observer_variance`,
  `minimal_conformal_covariance`: basis entries before transformations.
- `reference_state`: inherited auxiliary seed mass-squared setting (=1); not xi.

`SPECTRA.npz` covariance basis order is (minimal trace B, minimal observer source P,
conformal observer source E), resolved by external angular order J. Heat-kernel
weights are exp[-delta² J(J+1)]. Each raw block stores the l contribution and
underlying (B0,R,D,E0) per-partner Gram matrix. The finite xi results are exact
linear transformations of this SINGLE correlated dataset at each window.

`truncated_angular` omits the last 32 field orders. `truncated_momentum` removes
nodes above 0.8 of the explicit K while keeping original quadrature weights: it
is a tail test, not a separately optimized quadrature. The two `_grid` campaigns
supply independent momentum/time grids and fresh field histories at the outer
and inner locations. `middle_time` changes time quadrature only.
