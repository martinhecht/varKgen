# varKgen

Calibrated VAR(K) data-generating processes for interval-specific
cross-lagged effects ("varK generator").

The package refactors `simulation_script_20260222.R` into a documented,
tested R package. It constructs bivariate VAR(K) processes whose implied
population CLPM coefficients `b21(delta)` and `b12(delta)` across intervals
`delta = 1, ..., delta_max` match user-specified target profiles, jointly
with the population autoregressive projection effects `b11(Delta_1)` and
`b22(Delta_1)`. The target intervals are the fixed grid `Delta_j = j`, so
`Delta_1 = 1` time unit. All VAR coefficients are freely optimized; none is
held fixed.

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

The test suite in `tests/testthat/` covers the mathematical identities
(Lyapunov residual, the `K = 1` / `Delta_1` identity, agreement of the exact
and iterative solvers, `lm()` equivalence, the loss formula, and end-to-end
calibration).

## Quickstart

```r
library(varKgen)

targets_b21 <- c(0.30, 0.25, 0.15, 0.12, 0.10, 0.08, 0.06, 0.05, 0.04, 0.03)
targets_b12 <- c(0.15, 0.20, 0.10, 0.10, 0.08, 0.07, 0.05, 0.04, 0.03, 0.02)

sol <- calibrate_varK(targets_b21, targets_b12, K = 5,
                      target_b11_1 = 0.5, target_b22_1 = 0.5,
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
| `unpack_theta_fixedAR1()`         | internal `unpack_theta()` (4K layout, no fixed lag-1 diagonal; returns `b21`, `b12` instead of `b`, `c`) |
| `loss_fixedAR1()`                 | internal `loss_varK()` (**not** equivalent: two additional free VAR coefficients and two additional autoregressive targets at `Delta_1` enter the loss; graded penalty in the unstable region) |
| `calibrate_fixedAR1_fast(theta0, K, delta_max, ...)` | `calibrate_varK(targets_b21, targets_b12, K, target_b11_1, target_b22_1, var1, var2, cov12, var_type, theta0, stab_thresh, restarts, maxit, jitter_sd, early_tol, seed, polish, parallel, cores, verbose)` (`theta0` optional via `default_theta0()`; the fixed `alpha11_1` / `alpha22_1` are gone) |
| `simulate_panel_stationary(N, T, mats, P)` | `simulate_panel_stationary(N, T_obs, mats, P_stationary, seed = NULL)` (exact PSD-factor draws, no MASS) |
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
4. The stationary covariance is symmetrized before factorization/simulation
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
  plausible; the exact factor depends on the spectral radius and has not
  been benchmarked systematically. Results carry a residual check;
  the old iteration remains available as `method = "iterate"`.
- Implied profiles use cumulative 2-row blocks of `F^delta` instead of full
  matrix powers per interval (O(delta_max) instead of O(delta_max^2) matrix
  products).
- Graded instability penalty `1e9 + 1e7 * (rho - stab_thresh)` gives
  Nelder-Mead a direction signal in the infeasible region; at feasible
  points the loss is unchanged.
- Simulation draws use precomputed exact PSD factors instead of
  `MASS::mvrnorm` per time step (removes the MASS dependency and one
  eigendecomposition per step).
- `fit_clpm_ols()` uses one QR with a two-column response.
- Optional `polish = TRUE` (BFGS from the best Nelder-Mead solution, kept
  only if it improves the loss) and `parallel = TRUE` (PSOCK cluster over
  restarts; requires the installed package; no early stopping).

## Reproducibility notes

- For a given `seed`, the jittered start values are reproducible within
  version 0.3.0 (the draws are taken upfront; Nelder-Mead is
  deterministic). They are **not** identical to those of the original
  script or of earlier package versions, because the length and layout of
  the parameter vector have changed: the vector now holds `4K` elements
  with layout `(ax[1..K], ay[1..K], b21[1..K], b12[1..K])`, so the same
  `rnorm` sequence is mapped to different coefficients. Results are
  therefore not comparable coefficient by coefficient across versions.
- The simulation RNG stream differs from the original script, and from
  package versions up to 0.2.0, for the same seed: draws now use an exact
  symmetric PSD factor instead of a Cholesky factor with diagonal jitter.
  The distribution of the simulated data is exactly the specified one, also
  for singular positive semi-definite covariance matrices, where the old
  jittered Cholesky silently simulated from `M + epsilon * I`.
- The loss is a sum of `2 * delta_max + 2` squared deviations (the
  cross-lag profiles plus the two autoregressive projection effects at
  the first target interval `Delta_1`); each individual target coefficient
  contributes one squared deviation with unit weight, no block weighting is
  applied. The default
  `early_tol = 1e-10` corresponds to a root mean squared calibration error
  of `sqrt(early_tol / (2 * delta_max + 2))`, i.e. about 2.13e-6 for
  `delta_max = 10` (22 targets: 10 for `b21`, 10 for `b12`, one for
  `b11(Delta_1)`, one for `b22(Delta_1)`). Set `early_tol = 0` to disable
  early
  stopping (useful for optimizer studies on `restart_log`).
- If `2 * delta_max + 2 > 4K`, a message notes that the targets can in
  general only be matched approximately.
- `target_b11_1` and `target_b22_1` are population projection effects, i.e.
  the diagonal of `B(Delta_1)`, not the diagonal of the first VAR matrix
  `A1`. For `K = 1` the two coincide; for `K > 1` they generally differ,
  though equality can still occur in special cases such as all higher VAR
  lag matrices being zero.
  The two coincide only for `K = 1`.
