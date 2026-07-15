#' Verify a calibrated VAR(K) generator on simulated data
#'
#' Convenience wrapper reproducing the verification step of the original
#' script: simulates one stationary panel from a calibration result, fits the
#' interval-specific CLPM regressions by pooled OLS for each requested
#' interval, and tabulates the estimates against both the target values and
#' the implied population values.
#'
#' @details
#' The comparison against the \emph{implied} population coefficients is the
#' statistically relevant one for verifying the simulation and estimation
#' pipeline (only sampling error should remain). The comparison against the
#' \emph{targets} additionally contains the calibration error, i.e. the
#' deviation of the implied population values from the targets. Keeping the
#' two apart avoids misattributing calibration error to the estimator.
#'
#' @param calibration An object of class \code{"varK_calibration"} from
#'   \code{\link{calibrate_varK}}.
#' @param N Number of persons to simulate.
#' @param T_obs Number of time points to simulate.
#' @param deltas Intervals to evaluate; defaults to
#'   \code{1:delta_max} of the calibration. All values must be smaller than
#'   \code{T_obs} and at most \code{delta_max}.
#' @param seed Optional seed for the simulation.
#' @param return_data If \code{TRUE}, the simulated panel is returned as
#'   element \code{sim} (can be large).
#'
#' @return A list with elements \code{estimates} (4 x length(deltas) matrix
#'   with rows \code{b11_hat}, \code{b21_hat}, \code{b22_hat},
#'   \code{b12_hat}) and \code{table} (data frame with targets, implied
#'   values, estimates and their differences per interval), plus optionally
#'   \code{sim}.
#'
#' @examples
#' \donttest{
#' targets_b21 <- c(0.30, 0.25, 0.15)
#' targets_b12 <- c(0.15, 0.20, 0.10)
#' sol <- calibrate_varK(targets_b21, targets_b12, K = 2,
#'                       alpha11_1 = 0.5, alpha22_1 = 0.5,
#'                       restarts = 10, maxit = 5000, verbose = FALSE)
#' v <- verify_varK(sol, N = 2000, T_obs = 20, seed = 1)
#' v$table
#' }
#'
#' @export
verify_varK <- function(calibration, N, T_obs, deltas = NULL, seed = NULL,
                        return_data = FALSE) {
  stopifnot(inherits(calibration, "varK_calibration"))
  delta_max <- calibration$settings$delta_max
  if (is.null(deltas)) deltas <- seq_len(delta_max)
  stopifnot(all(deltas >= 1), all(deltas == round(deltas)),
            all(deltas <= delta_max), all(deltas < T_obs))

  sim <- simulate_panel_stationary(N, T_obs, calibration$mats,
                                   calibration$P, seed = seed)

  est <- matrix(NA_real_, nrow = 4L, ncol = length(deltas),
                dimnames = list(c("b11_hat", "b21_hat", "b22_hat", "b12_hat"),
                                paste0("delta=", deltas)))
  for (j in seq_along(deltas)) {
    dfD <- make_delta_pairs_overlapping(sim$Y1, sim$Y2, deltas[j])
    co <- fit_clpm_ols(dfD)
    est[, j] <- co[c("b11", "b21", "b22", "b12")]
  }

  idx <- deltas
  tab <- data.frame(
    delta       = deltas,
    target_b21  = calibration$targets_b21[idx],
    implied_b21 = as.numeric(calibration$implied["b21", idx]),
    est_b21     = as.numeric(est["b21_hat", ]),
    target_b12  = calibration$targets_b12[idx],
    implied_b12 = as.numeric(calibration$implied["b12", idx]),
    est_b12     = as.numeric(est["b12_hat", ])
  )
  tab$diff_b21_vs_implied <- tab$est_b21 - tab$implied_b21
  tab$diff_b21_vs_target  <- tab$est_b21 - tab$target_b21
  tab$diff_b12_vs_implied <- tab$est_b12 - tab$implied_b12
  tab$diff_b12_vs_target  <- tab$est_b12 - tab$target_b12

  out <- list(estimates = est, table = tab)
  if (return_data) out$sim <- sim
  out
}
