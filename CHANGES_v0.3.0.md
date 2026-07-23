# varKgen 0.3.0 (update relative to 0.2.0)

## Breaking: no fixed lag-1 autoregressive VAR coefficients

The arguments `alpha11_1` and `alpha22_1` of `calibrate_varK()` have been
removed. They fixed the two diagonal elements of the first VAR coefficient
matrix. Every element of every VAR coefficient matrix is now freely
optimized. Counting only VAR coefficients, the number of free parameters
rises from `4K - 2` to `4K`. Counting the full optimization dimension, i.e.
including the free process error correlation that version 0.1.0 still had
and 0.2.0 removed, it rises from `4K - 1` (0.1.0) to `4K`. The internal
`theta` layout becomes `(ax[1..K], ay[1..K], b21[1..K], b12[1..K])`.

The removed arguments are not accepted as aliases for the new target
arguments, because they denoted a conceptually different quantity. Passing
them raises an informative error.

## New: target autoregressive projection effects at `Delta_1`

`calibrate_varK()` gains the required arguments `target_b11_1` and
`target_b22_1`. They specify the population-level autoregressive projection
effects at the first target interval `Delta_1 = 1` time unit, that is the
diagonal elements of the population projection matrix `B(Delta_1)`, each
controlling for the respective other process.

## Notation of the target intervals

The target intervals are the fixed grid
`Delta = (Delta_1, ..., Delta_delta_max) = (1, ..., delta_max)`, so
`Delta_j = j` and in particular `Delta_1 = 1` time unit. The subscript
denotes the position in this fixed grid; a freely chosen interval vector is
not supported. `targets_b21[j]` is the target for `b21(Delta_j)`, and
`target_b11_1`, `target_b22_1` are the targets for `b11(Delta_1)` and
`b22(Delta_1)`.

For `K = 1`, `B(Delta_1)` equals the VAR(1) coefficient matrix. For
`K > 1` the projection effects generally differ from the corresponding
elements of the first VAR coefficient matrix because they depend on the
complete VAR(K) process; equality can still occur in special cases, for
instance when all higher VAR lag matrices are zero. The projection effects
come from the single canonical routine `implied_clpm_profile()`, whose first
column is `Delta_1`; no second projection routine exists, so the loss, the
result object and `implied_clpm_profile()` cannot diverge.

Each individual target coefficient contributes one squared deviation with
unit weight; no separate block weighting is applied. The loss is therefore a
sum of `2 * delta_max + 2` squared deviations.

## `verify_varK()` reports the autoregressive comparison at `Delta_1`

The returned list gains the element `ar1`, a one-row data frame with
`target_b11_1`, `implied_b11_1`, `estimated_b11_1`, the corresponding `b22`
triple and their differences. This separates calibration error (implied vs.
target) from simulation and estimation error (estimated vs. implied). The
interval of one time unit is always evaluated, even if `1` is not contained
in `deltas`. The main table gains implied and estimated autoregressive
columns for every requested interval.

## Optimizer diagnostics are fully retained

`restart_log` now carries `convergence`, `fn_evals` and `message` per
restart, and the new matrix `restart_theta` stores the full pre-polish
Nelder-Mead endpoint of every restart. Validity, spectral radius, smallest
eigenvalue of the process error covariance, implied autoregressive patterns
and all near-best variability measures can therefore be reconstructed after
the fact.

The optimizer stages are kept apart: `theta_nm` / `loss_nm` hold the best
pre-polish Nelder-Mead solution, on which hits and near-best variability are
defined, while `theta_polished`, `loss_polished`, `polish_convergence` and
`polish_improved` describe the BFGS step. `theta` and `loss` remain the
solution finally used.

## Fewer eigendecompositions per loss evaluation

`stationary_cov()` and `solve_process_error()` gain `check_stability`
(default `TRUE`). `loss_varK()` establishes stability itself and calls them
with `check_stability = FALSE`, removing one `2K x 2K` eigendecomposition
per loss evaluation.

## Stricter input validation

`restarts`, `maxit`, `jitter_sd`, `early_tol`, `seed`, `stab_thresh`,
`polish`, `parallel`, `verbose` and `cores` are checked for length, type,
finiteness and, where applicable, integrality up front, so e.g.
`early_tol = NA` fails immediately instead of during optimization.
`default_theta0()`, `stationary_cov()`, `solve_process_error()`,
`simulate_panel_stationary()`, `verify_varK()` and the print method validate
their own arguments as well. `settings` additionally records the start
vector `theta0`, `cores_requested`, `cores_used`, `package_version`,
`R_version` and `RNGkind()`. This covers the run configuration; bit-identical
reproduction may additionally depend on the BLAS/LAPACK implementation, so
study code should still store a `sessionInfo()` protocol.

## Exact covariance factor for simulation

The internal Cholesky-with-jitter helper is replaced by `cov_factor_psd()`,
an exact symmetric factor based on the eigendecomposition. `chol()` fails
for exactly singular positive semi-definite matrices, and the previous
diagonal jitter silently simulated from `M + epsilon * I` instead of `M`.
Singular covariance matrices are reachable: `make_mats_varK()` admits
`abs(pe_cov12) = sqrt(pe_var1 * pe_var2)`, and calibrated process error
covariances can be numerically near-singular. Eigenvalues within a relative
tolerance of zero are now set to zero exactly and clearly negative
eigenvalues raise an error. For a given seed the simulated realizations
therefore differ from earlier versions, while the distribution is exactly
the specified one.

## Restart failures no longer abort a calibration

`stats::optim()` inside a restart is wrapped in `tryCatch()`. A failed
restart is logged with `status = "error"`, `value = Inf`, `NA` convergence
diagnostics, the condition message, and an all-`NA` row in `restart_theta`,
so it can be excluded from hit and near-best analyses. This matters for the
planned optimizer study with more than 100,000 restarts.

## Failure handling and diagnostics

A restart that raises an R error, or that returns a non-finite or malformed
result, is logged with `status = "error"` or `"invalid"`, `value = Inf` and
an all-`NA` row in `restart_theta`. The best solution is selected explicitly
from the usable restarts rather than by incremental update, and a run in
which no restart is usable now stops with a clear error instead of silently
returning an object built from `theta0`. The final warning also triggers for
a non-finite best loss. The single restart is implemented in the internal,
testable helper `run_restart()`.

The polish stage reports `polish_attempted`, `polish_status` (`"not
requested"`, `"skipped (penalty solution)"`, `"completed"`, `"error"`,
`"invalid"`), `polish_message`, `polish_fn_evals` and `polish_seconds`, so a
failed BFGS call is distinguishable from one that was never attempted.

`early_tol` must stay below the penalty sentinel and early stopping only
triggers on genuine fit values. The sum of squared targets must stay below
1e6 so that a valid fit cannot reach the sentinel scale. Numeric control
arguments are type-checked with `is.numeric()` and rejected when logical.
`parallel::detectCores()` returning `NA` falls back to one worker, `cores` is
capped at the number of restarts, the master's `.libPaths()` are propagated
to the PSOCK workers so that the package is loadable in project, check or
user libraries, and the cluster is closed directly after `parLapply()`
rather than at function exit.

## Simulation is explicit about what it simulates

`simulate_panel_stationary()` now requires `mats$Q` to be zero outside the
top-left 2 x 2 block instead of silently ignoring the rest, checks that the
covariance matrices are matrices and symmetric and that the state dimension
is even, and returns `psd_correction` for `P_stationary` and the process
error block. `solve_process_error()` additionally verifies the full Lyapunov
residual and computes the smallest eigenvalue with the symmetric eigensolver
instead of an overflow-prone closed form. Relative numerical tolerances no
longer use a fixed floor of 1, so they remain meaningful for covariance
scales far below 1. `fit_clpm_ols()` requires finite numeric columns and
rejects rank-deficient designs, and `make_delta_pairs_overlapping()` returns
`id`, `p` and `q` so that cluster-robust standard errors are actually
usable.

## Naming of the coefficient vectors

`theta`, `theta_nm` and `theta_polished` carry the canonical names
`ax1..axK`, `ay1..ayK`, `b211..b21K`, `b121..b12K`, matching the column
names of `restart_theta`. The decomposed blocks in `params` (`ax`, `ay`,
`b21`, `b12`) are returned unnamed, since they are indexed by lag.

## Robustness

All calls into `stats` now use the `stats::` prefix, including `rnorm()`, so
regenerating `NAMESPACE` with `devtools::document()` cannot break the
package. The roxygen comments additionally carry `@importFrom stats rnorm`,
so a regenerated `NAMESPACE` reproduces the current one.

## Other changes

- `default_theta0()` gains the arguments `target_b11_1` and `target_b22_1`,
  which are used as start values for `ax[1]` and `ay[1]`, and returns a
  vector of length `4K`.
- `calibrate_varK()` returns `target_b11_1`, `target_b22_1`,
  `implied_b11_1`, `implied_b22_1` and the deviations `error_b11_1`,
  `error_b22_1` (implied minus target); `settings` reports the new targets
  instead of the removed `alpha11_1` / `alpha22_1`.
- The print method reports target, implied and error for `b11(Delta_1)`
  and `b22(Delta_1)`.
