# varKgen 0.2.0 (update relative to 0.1.0)

## New: full covariance specification with `var_type`

`calibrate_varK()` gains the arguments `var1`, `var2`, `cov12` (defaults 1,
1, 0) and `var_type = c("process", "process_error")`:

- `var_type = "process"` (new default): `var1`, `var2`, `cov12` fix the
  stationary process covariance (top-left 2 x 2 block of `P`). The required
  process error covariance is derived exactly in every loss evaluation by a
  new internal linear solver (`solve_process_error()`, in
  `R/matrix_utils.R`), which also handles unequal variances. With
  `var1 = var2 = 1`, standardized and unstandardized cross-lagged
  coefficients coincide, so targets can be given in the standardized metric.
  Inadmissible targets (derived process error covariance not positive
  semi-definite) receive a graded penalty during optimization; if the final
  solution is inadmissible or unstable, a warning is issued.
- `var_type = "process_error"`: the three values fix the process error
  covariance directly (previous behaviour, minus the free correlation).

## Removed

- The free process error correlation parameter (`rho_raw` / `tanh`) is gone
  in both modes; the full contemporaneous covariance is always set by the
  user. The parameter vector `theta` now has length `4K - 2` with layout
  `(ax[2..K], ay[2..K], b21[1..K], b12[1..K])`. `default_theta0()` returns
  this length.
- The arguments `var_u`, `var_e`, `cov_ue` and the term "innovation" are
  removed everywhere; the terminology is "process error" throughout.
  `make_mats_varK()` now takes `pe_var1`, `pe_var2`, `pe_cov12`.

## Return object additions

- `process_cov`: achieved stationary process covariance (2 x 2).
- `params` now contains `pe_var1`, `pe_var2`, `pe_cov12` (the process error
  covariance actually used, derived in "process" mode).
- `settings` records `var1`, `var2`, `cov12`, `var_type`.
- The print method shows the mode and both covariances.

## Files changed in this update

R/matrix_utils.R (adds `solve_process_error()`, wording),
R/companion.R (adds internal `companion_F()`, `pe_*` arguments),
R/calibrate.R (full rewrite of the calibration layer),
R/simulate.R (wording only),
man/make_mats_varK.Rd, man/calibrate_varK.Rd, man/default_theta0.Rd,
man/stationary_cov.Rd, man/simulate_panel_stationary.Rd,
tests/testthat/test-companion.R, test-stationary.R, test-implied.R,
test-simulate.R, test-calibrate.R,
inst/scripts/reproduce_original_analysis.R.

Unchanged: NAMESPACE (no export changes), R/implied.R, R/clpm.R,
R/verify.R, remaining man pages, tests/testthat/test-clpm.R.

Please bump `Version:` in DESCRIPTION to 0.2.0 and run
`devtools::test()` / `R CMD check` once locally.
