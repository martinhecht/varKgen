# varKgen

Calibrated VAR(K) data-generating processes for interval-specific
cross-lagged effects ("varK generator").

The package refactors `simulation_script_20260222.R` into a documented,
tested R package. It constructs bivariate VAR(K) processes whose implied
population CLPM coefficients `b21(delta)` and `b12(delta)` across intervals
`delta = 1, ..., delta_max` match user-specified target profiles, with the
lag-1 autoregressions of the DGP held fixed.

## Installation

From the unzipped source directory:

```r
install.packages("path/to/varKgen", repos = NULL, type = "source")
# or
devtools::install_local("path/to/varKgen")
```

Then run the unit tests once locally:

```r
devtools::test("path/to/varKgen")   # or: testthat::test_local()
```

Note: this package was assembled in a sandbox without an R installation,
so the code was reviewed statically but not executed there. The test suite
in `tests/testthat/` covers the mathematical identities (Lyapunov residual,
K = 1 / delta = 1 identity, agreement of the exact and iterative solvers,
`lm()` equivalence, calibration smoke test); please run it once before
production use. Also adjust the placeholder e-mail address in
`DESCRIPTION`.

## Quickstart

```r
library(varKgen)

targets_b21 <- c(0.30, 0.25, 0.15, 0.12, 0.10, 0.08, 0.06, 0.05, 0.04, 0.03)
targets_b12 <- c(0.15, 0.20, 0.10, 0.10, 0.08, 0.07, 0.05, 0.04, 0.03, 0.02)

sol <- calibrate_varK(targets_b21, targets_b12, K = 5,
                      alpha11_1 = 0.5, alpha22_1 = 0.5,
                      restarts = 120, seed = 1)
print(sol)

ver <- verify_varK(sol, N = 50000, T_obs = 30, seed = 1)
round(ver$table, 4)
```

`inst/scripts/reproduce_original_analysis.R` reproduces the full original
workflow with the original settings.

## Mapping: original script -> package

| Original                          | Package                                  |
|-----------------------------------|------------------------------------------|
| `make_mats_varK(ax, ay, b, c, ...)` | `make_mats_varK(ax, ay, b21, b12, ...)` (args renamed; `K = 1` fixed) |
| `solve_lyapunov(F, Q)`            | `stationary_cov(Fmat, Qmat)` (exact Kronecker solver; `method = "iterate"` keeps the old behaviour) |
| `coeffs_delta(F, P, delta)`       | `coeffs_delta(Fmat, P, delta)` (same output) plus `implied_clpm_profile()` for all intervals at once |
| `unpack_theta_fixedAR1()`         | internal `unpack_theta()` (returns `b21`, `b12` instead of `b`, `c`) |
| `loss_fixedAR1()`                 | internal `loss_varK()` (identical at stable points, graded penalty in the unstable region) |
| `calibrate_fixedAR1_fast(theta0, K, delta_max, ...)` | `calibrate_varK(targets_b21, targets_b12, K, ...)` (`theta0` optional via `default_theta0()`; extras: `polish`, `parallel`, `cores`) |
| `simulate_panel_stationary(N, T, mats, P)` | `simulate_panel_stationary(N, T_obs, mats, P_stationary, seed = NULL)` (Cholesky draws, no MASS) |
| `make_delta_pairs_overlapping()`  | same name (validates `delta < T_obs`)    |
| `fit_clpm_ols()`                  | same name and output (one QR via `lm.fit` instead of two `lm()` calls) |
| manual verification block         | `verify_varK()` (adds the estimate-vs-implied comparison) |

## Bug fixes

1. `make_mats_varK()` crashed for `K = 1` (`for (k in 2:K)` evaluates to
   `c(2, 1)` in R); now guarded.
2. `make_delta_pairs_overlapping()` silently produced invalid indices for
   `delta >= T`; now an informative error.
3. `solve(VarNow, ...)` in the implied-coefficient computation was
   unguarded and could abort the optimizer; now wrapped, degenerate points
   receive the penalty value.
4. The stationary covariance is symmetrized before Cholesky/simulation
   (numerical asymmetry could otherwise propagate).
5. Shadowing-prone names removed (`T` -> `T_obs`, function argument `F` ->
   `Fmat`, cross-lag argument `c` -> `b12`).
6. `restart_log` is preallocated instead of grown via `rbind()` (O(n)
   instead of O(n^2)).

## Performance changes

- Exact Lyapunov solve via `(I - F kron F) vec(P) = vec(Q)` instead of a
  fixed-point iteration whose iteration count explodes as the spectral
  radius approaches 1 (with `stab_thresh = 0.998`, several thousand
  iterations per loss evaluation were possible). Roughly one to two orders
  of magnitude faster per loss evaluation for persistent systems is
  plausible; the exact factor depends on the spectral radius and was not
  benchmarked here (no R in the sandbox). Results carry a residual check;
  the old iteration remains available as `method = "iterate"`.
- Implied profiles use cumulative 2-row blocks of `F^delta` instead of full
  matrix powers per interval (O(delta_max) instead of O(delta_max^2) matrix
  products).
- Graded instability penalty `1e9 + 1e7 * (rho - stab_thresh)` gives
  Nelder-Mead a direction signal in the infeasible region; at feasible
  points the loss is unchanged.
- Simulation draws use precomputed Cholesky factors instead of
  `MASS::mvrnorm` per time step (removes the MASS dependency and one
  eigendecomposition per step).
- `fit_clpm_ols()` uses one QR with a two-column response.
- Optional `polish = TRUE` (BFGS from the best Nelder-Mead solution, kept
  only if it improves the loss) and `parallel = TRUE` (PSOCK cluster over
  restarts; requires the installed package; no early stopping).

## Reproducibility notes

- For a given `seed`, the jittered start values are identical to the
  original script (same `rnorm` sequence; Nelder-Mead is deterministic).
  Final estimates are numerically equivalent but not bit-identical because
  the Lyapunov solution is now exact rather than iterative.
- The simulation RNG stream differs from the original for the same seed
  (Cholesky instead of eigendecomposition-based draws); the distribution of
  the simulated data is identical.
- The loss is a sum of `2 * delta_max` squared deviations; the default
  `early_tol = 1e-10` corresponds to a root mean squared calibration error
  of about 2.2e-6 for 20 targets. Set `early_tol = 0` to disable early
  stopping (useful for optimizer studies on `restart_log`).
- If `2 * delta_max > 4K - 1`, a message notes that the targets can in
  general only be matched approximately.
