# ---------------------------------------------------------------------------
# Re-run the analysis of simulation_script_20260222.R with varKgen >= 0.2.0
# ---------------------------------------------------------------------------
# Settings, seeds and workflow follow the original script where possible.
# Two structural differences to the original are unavoidable in >= 0.2.0:
#   - The process error correlation is no longer a free calibration
#     parameter: the full 2 x 2 covariance is always set by the user
#     (var1, var2, cov12) and interpreted according to var_type. Below it is
#     fixed at cov12 = 0 in "process_error" mode, which corresponds to the
#     original problem restricted to uncorrelated process errors. If the
#     original calibration found a non-zero error correlation, results will
#     differ accordingly; set cov12 to that value to match more closely.
#   - The number of free parameters is therefore 4K - 2 (was 4K - 1), and
#     jittered start vectors are shorter, so optimization paths differ.
# In addition, as before: the stationary covariance is solved exactly
# (Kronecker) instead of by fixed-point iteration, and simulation draws use
# Cholesky factors, so results are numerically equivalent but not
# bit-identical to the original script.
# ---------------------------------------------------------------------------

library(varKgen)

# ---- settings (as in the original script) ---------------------------------
delta_max <- 10
K <- 5

alpha11_1 <- 0.50
alpha22_1 <- 0.50

targets_b21 <- c(0.30, 0.25, 0.15, 0.12, 0.10, 0.08, 0.06, 0.05, 0.04, 0.03)
targets_b12 <- c(0.15, 0.20, 0.10, 0.10, 0.08, 0.07, 0.05, 0.04, 0.03, 0.02)

stopifnot(length(targets_b21) == delta_max,
          length(targets_b12) == delta_max)

# ---- calibration -----------------------------------------------------------
sol <- calibrate_varK(
  targets_b21 = targets_b21,
  targets_b12 = targets_b12,
  K = K,
  alpha11_1 = alpha11_1,
  alpha22_1 = alpha22_1,
  var1 = 1, var2 = 1, cov12 = 0,      # process error covariance, see header
  var_type = "process_error",
  stab_thresh = 0.998,
  restarts = 120,
  maxit = 25000,
  jitter_sd = 0.10,
  early_tol = 1e-10,
  seed = 1,
  polish = TRUE,      # set FALSE for the closest match to the original
  parallel = FALSE,   # TRUE distributes restarts over cores (pkg installed)
  verbose = TRUE
)

print(sol)

# Optional: save restart runtimes for reporting
# write.csv(sol$restart_log, "optimizer_restart_log.csv", row.names = FALSE)

# ---- verification (as in the original script: N = 50000, T = 30) ----------
ver <- verify_varK(sol, N = 50000, T_obs = 30, seed = 1)

cat("\nOLS estimates from one simulated dataset:\n")
print(round(ver$estimates, 4))

cat("\nEstimates vs. targets and vs. implied population values:\n")
print(round(ver$table, 6), row.names = FALSE)

cat(sprintf("\nMax abs diff (est - target):  b21 = %.6f, b12 = %.6f\n",
            max(abs(ver$table$diff_b21_vs_target)),
            max(abs(ver$table$diff_b12_vs_target))))
cat(sprintf("Max abs diff (est - implied): b21 = %.6f, b12 = %.6f\n",
            max(abs(ver$table$diff_b21_vs_implied)),
            max(abs(ver$table$diff_b12_vs_implied))))
