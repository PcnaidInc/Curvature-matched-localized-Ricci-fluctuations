# Result and evidence files
- `results/FINAL_VALIDATION.json`: authoritative pass/fail status, measured numerical
  sensitivities, all result rows and conditional screens. Never infer acceptance
  from a successful program exit alone.
- `results/INWARD_LOCALIZED_RESULTS.csv`: center, observable, angular convention,
  sampling support, source/reference variance and RMS, and gravitational normalization.
  `trace` maps to delta Ricci scalar; `time_projection` maps to delta R_uu;
  `energy_projection` maps to delta G_uu. The latter is supplementary.
- `uniform_sphere=1, delta_rad=0` denotes the UNIFORM kernel, not a zero-size detector.
  Nonuniform delta values are heat-kernel diffusion widths in radians.
- `relative_curvature_RMS_per_epsilon` is 8pi sqrt(variance)/sampled sqrt(K).
  Multiply by epsilon=G hbar/M^2 for a dimensionless ratio; multiply by 100 for percent.
  It is not a fluctuation of K or a percentage of all gravity.
- `source_to_reference_RMS` compares different state/geometry calculations using
  matched filters. It is not a partition into additive independent noise sources.
- `ANGULAR_SPECTRA.csv`: ordered [T,pi^2] source spectral covariance and reference.
  The Ricci-scalar/time-projection cross covariance has the opposite sign to the
  stored source cross covariance. The overall gravitational factor is (8piG)^2.
- `CONDITIONAL_SCALE_SCREENS.csv`: algebraic rescalings, not new mode runs.
  N=10,30,100 are chosen diagnostic separation factors, not proven validity bounds.
- `joint_####.npz/.json`: exact saved angular pair moments, covariances and hash.
  Six moment coefficients describe trace/time cross covariance. Next six are frozen
  reference moments. Angular partners and multiplicity are retained.
- `final_states/`, `mode_records/`: exact end state, integration count/time and source
  fingerprint. `*_histories/` stores separately validated approximate phase-resolved
  histories on the SAME time interval. They are not freshly evolved states.
- The original dense `mode_cache/` arrays are transient and omitted from the final ZIP
  to avoid a multi-gigabyte duplicate. The accepted covariance, exact final states,
  checked portable histories, configuration and full reproduction code are retained.
- `pilot` is a deliberately coarse calculation and may be rejected. It is never
  substituted for the finest result. `main`/`fine` field-angular truncations are
  compared using both l and p cutoffs, not only the outer loop index.
