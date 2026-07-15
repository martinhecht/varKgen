#' Unpack the parameter vector of the calibration problem
#'
#' Internal helper. The parameter vector \code{theta} has length
#' \code{4K - 1} and layout
#' \code{(ax[2..K], ay[2..K], b21[1..K], b12[1..K], rho_raw)}; the lag-1
#' autoregressive coefficients are fixed at \code{alpha11_1} and
#' \code{alpha22_1} and are not free parameters. The innovation covariance is
#' parameterized as \code{cov_ue = tanh(rho_raw) * sqrt(var_u * var_e)}, which
#' keeps the innovation correlation inside (-1, 1).
#'
#' @param theta Numeric parameter vector of length \code{4K - 1}.
#' @param K Number of lags.
#' @param alpha11_1,alpha22_1 Fixed lag-1 autoregressive coefficients.
#' @param var_u,var_e Innovation variances.
#'
#' @return List with elements \code{ax}, \code{ay}, \code{b21}, \code{b12},
#'   \code{cov_ue}.
#' @keywords internal
unpack_theta <- function(theta, K, alpha11_1, alpha22_1, var_u = 1, var_e = 1) {
  stopifnot(length(theta) == 4L * K - 1L)
  idx <- 1L

  ax <- numeric(K)
  ay <- numeric(K)
  ax[1L] <- alpha11_1
  ay[1L] <- alpha22_1

  if (K > 1L) {
    ax[2L:K] <- theta[idx:(idx + K - 2L)]; idx <- idx + (K - 1L)
    ay[2L:K] <- theta[idx:(idx + K - 2L)]; idx <- idx + (K - 1L)
  }

  b21 <- theta[idx:(idx + K - 1L)]; idx <- idx + K
  b12 <- theta[idx:(idx + K - 1L)]; idx <- idx + K

  rho_raw <- theta[idx]
  cov_ue <- tanh(rho_raw) * sqrt(var_u * var_e)

  list(ax = ax, ay = ay, b21 = b21, b12 = b12, cov_ue = cov_ue)
}

#' Calibration loss
#'
#' Sum of squared deviations between the implied population cross-lagged
#' coefficients \code{b21(delta)}, \code{b12(delta)} for
#' \code{delta = 1, ..., delta_max} and the target profiles. Unstable or
#' numerically degenerate parameter values receive a large penalty.
#'
#' @details
#' Differences to the original script: (i) the stationary covariance is
#' obtained from the exact Kronecker solver instead of the fixed-point
#' iteration, which is substantially faster for persistent systems;
#' (ii) the instability penalty is graded
#' (\code{1e9 + 1e7 * (max |eigenvalue| - stab_thresh)}) instead of a flat
#' \code{1e9}, so that Nelder-Mead receives a direction signal back towards
#' the stable region instead of a plateau. At all stable (feasible) points the
#' loss is identical to the original definition.
#'
#' @param theta Parameter vector, see \code{\link{unpack_theta}}.
#' @param K Number of lags.
#' @param delta_max Largest interval; \code{length(targets_b21)} must equal
#'   \code{delta_max}.
#' @param alpha11_1,alpha22_1 Fixed lag-1 autoregressive coefficients.
#' @param targets_b21,targets_b12 Target profiles of length \code{delta_max}.
#' @param var_u,var_e Innovation variances.
#' @param stab_thresh Stability threshold for the spectral radius of the
#'   companion matrix.
#'
#' @return A single non-negative loss value.
#' @keywords internal
loss_varK <- function(theta, K, delta_max, alpha11_1, alpha22_1,
                      targets_b21, targets_b12,
                      var_u = 1, var_e = 1, stab_thresh = 0.998) {
  if (!all(is.finite(theta))) return(1e12)

  par <- unpack_theta(theta, K, alpha11_1, alpha22_1, var_u, var_e)
  mats <- make_mats_varK(par$ax, par$ay, par$b21, par$b12,
                         var_u = var_u, var_e = var_e, cov_ue = par$cov_ue)

  eig <- eigen(mats$F, only.values = TRUE)$values
  if (!all(is.finite(eig))) return(1e12)
  maxmod <- max(Mod(eig))
  if (maxmod >= stab_thresh) return(1e9 + 1e7 * (maxmod - stab_thresh))

  P <- stationary_cov(mats$F, mats$Q, method = "kron")
  if (is.null(P)) return(1e9)

  prof <- implied_clpm_profile(mats$F, P, delta_max)
  if (!all(is.finite(prof))) return(1e9)

  sum((prof["b21", ] - targets_b21)^2) + sum((prof["b12", ] - targets_b12)^2)
}

#' Default start values for the calibration
#'
#' Reproduces the heuristic of the original script: small positive values for
#' the higher-lag autoregressive coefficients, the first cross-lagged target
#' as start value for the lag-1 cross effects, small positive values for the
#' remaining cross-lag coefficients, and zero for the raw innovation
#' correlation parameter.
#'
#' @param K Number of lags.
#' @param targets_b21,targets_b12 Target profiles (only the first element of
#'   each is used).
#'
#' @return Numeric vector of length \code{4K - 1}.
#'
#' @examples
#' default_theta0(K = 5,
#'                targets_b21 = c(0.30, 0.25, 0.15, 0.12, 0.10),
#'                targets_b12 = c(0.15, 0.20, 0.10, 0.10, 0.08))
#'
#' @export
default_theta0 <- function(K, targets_b21, targets_b12) {
  stopifnot(K >= 1L, K == round(K),
            length(targets_b21) >= 1L, length(targets_b12) >= 1L)
  c(
    if (K > 1L) rep(0.05, K - 1L) else numeric(0),
    if (K > 1L) rep(0.05, K - 1L) else numeric(0),
    c(targets_b21[1L], rep(0.03, K - 1L)),
    c(targets_b12[1L], rep(0.03, K - 1L)),
    0
  )
}

#' Calibrate a VAR(K) data-generating process to target cross-lag profiles
#'
#' Searches, via multi-start Nelder-Mead (optionally followed by a BFGS
#' polishing step), for VAR(K) parameters such that the implied population
#' cross-lagged CLPM coefficients \code{b21(delta)} and \code{b12(delta)}
#' across intervals \code{delta = 1, ..., delta_max} match the supplied
#' target profiles. The lag-1 autoregressive coefficients of the
#' data-generating process are held fixed at \code{alpha11_1} and
#' \code{alpha22_1}.
#'
#' @details
#' The number of free parameters is \code{4K - 1} while the number of target
#' values is \code{2 * delta_max}. If \code{2 * delta_max > 4K - 1}, the
#' targets can in general only be matched approximately and the attainable
#' loss has a positive floor; a message points this out.
#'
#' For a given \code{seed}, the jittered start values are identical to those
#' of the original script (they are drawn upfront in the same order;
#' Nelder-Mead itself does not consume random numbers). Loss values at stable
#' points are identical to the original definition, but because the
#' stationary covariance is now computed exactly rather than iteratively,
#' optimization paths and final estimates can differ in the last decimals;
#' results are numerically equivalent, not bit-identical.
#'
#' \code{parallel = TRUE} distributes the restarts over a PSOCK cluster.
#' This requires the package to be installed (workers load it via its
#' namespace) and disables early stopping.
#'
#' @param targets_b21 Target profile for the cross effect process 2 -> 1,
#'   length \code{delta_max}.
#' @param targets_b12 Target profile for the cross effect process 1 -> 2,
#'   same length as \code{targets_b21}.
#' @param K Number of lags of the VAR(K) data-generating process.
#' @param alpha11_1,alpha22_1 Fixed lag-1 autoregressive coefficients of the
#'   data-generating process. Note that these fix the DGP coefficients, not
#'   the implied CLPM autoregressions \code{b11(1)}, \code{b22(1)}, which
#'   will generally differ once higher lags are present.
#' @param var_u,var_e Innovation variances (default 1).
#' @param theta0 Optional start vector of length \code{4K - 1}; defaults to
#'   \code{\link{default_theta0}}.
#' @param stab_thresh Stability threshold for the spectral radius (default
#'   0.998).
#' @param restarts Number of jittered Nelder-Mead restarts (default 120).
#' @param maxit Maximum Nelder-Mead iterations per restart (default 25000).
#' @param jitter_sd Standard deviation of the normal jitter added to
#'   \code{theta0} for each restart (default 0.10).
#' @param early_tol Early-stopping tolerance on the loss (default 1e-10).
#'   Note the loss is a sum of \code{2 * delta_max} squared deviations, so
#'   for 20 targets this corresponds to a root mean squared calibration
#'   error of about 2.2e-6.
#' @param seed Random seed for the restart jitter (default 1).
#' @param polish If \code{TRUE} (default), a BFGS run from the best
#'   Nelder-Mead solution is attempted and kept only if it improves the loss.
#' @param parallel If \code{TRUE}, restarts run in parallel (see Details).
#' @param cores Number of worker processes if \code{parallel = TRUE};
#'   defaults to \code{parallel::detectCores() - 1}.
#' @param verbose Print progress information (default \code{TRUE}).
#'
#' @return An object of class \code{"varK_calibration"}: a list with elements
#'   \code{theta} (best parameter vector), \code{loss}, \code{implied}
#'   (4 x delta_max matrix of implied population coefficients, rows
#'   \code{b11}, \code{b21}, \code{b22}, \code{b12}), \code{params}
#'   (unpacked parameters \code{ax}, \code{ay}, \code{b21}, \code{b12},
#'   \code{cov_ue}), \code{mats} (companion matrices \code{F}, \code{Q}),
#'   \code{P} (stationary covariance), \code{targets_b21},
#'   \code{targets_b12}, \code{settings}, \code{runtime_seconds} and
#'   \code{restart_log}.
#'
#' @examples
#' \donttest{
#' targets_b21 <- c(0.30, 0.25, 0.15, 0.12, 0.10)
#' targets_b12 <- c(0.15, 0.20, 0.10, 0.10, 0.08)
#' sol <- calibrate_varK(targets_b21, targets_b12, K = 3,
#'                       alpha11_1 = 0.5, alpha22_1 = 0.5,
#'                       restarts = 10, maxit = 5000, verbose = FALSE)
#' print(sol)
#' }
#'
#' @export
calibrate_varK <- function(targets_b21, targets_b12, K,
                           alpha11_1, alpha22_1,
                           var_u = 1, var_e = 1,
                           theta0 = NULL,
                           stab_thresh = 0.998,
                           restarts = 120L,
                           maxit = 25000L,
                           jitter_sd = 0.10,
                           early_tol = 1e-10,
                           seed = 1L,
                           polish = TRUE,
                           parallel = FALSE,
                           cores = NULL,
                           verbose = TRUE) {

  stopifnot(length(targets_b21) == length(targets_b12),
            length(targets_b21) >= 1L,
            all(is.finite(targets_b21)), all(is.finite(targets_b12)),
            length(K) == 1L, K >= 1L, K == round(K),
            is.finite(alpha11_1), is.finite(alpha22_1),
            var_u > 0, var_e > 0,
            stab_thresh > 0, stab_thresh < 1,
            restarts >= 1L, maxit >= 1L, jitter_sd >= 0)
  delta_max <- length(targets_b21)
  n_par <- 4L * K - 1L

  if (is.null(theta0)) theta0 <- default_theta0(K, targets_b21, targets_b12)
  stopifnot(length(theta0) == n_par, all(is.finite(theta0)))

  if (2L * delta_max > n_par && verbose) {
    message(sprintf(paste0(
      "Note: %d target values vs. %d free parameters (4K - 1). ",
      "Targets can in general only be matched approximately ",
      "(attainable loss may have a positive floor)."),
      2L * delta_max, n_par))
  }

  loss_args <- list(K = K, delta_max = delta_max,
                    alpha11_1 = alpha11_1, alpha22_1 = alpha22_1,
                    targets_b21 = targets_b21, targets_b12 = targets_b12,
                    var_u = var_u, var_e = var_e, stab_thresh = stab_thresh)

  t_start <- Sys.time()

  # Draw all jittered starts upfront (same rnorm sequence as drawing them
  # one-by-one in a loop, since Nelder-Mead does not use the RNG).
  set.seed(seed)
  starts <- lapply(seq_len(restarts),
                   function(r) theta0 + rnorm(n_par, sd = jitter_sd))

  run_one <- function(th) {
    t0 <- proc.time()[["elapsed"]]
    fit <- do.call(stats::optim,
                   c(list(par = th, fn = loss_varK,
                          method = "Nelder-Mead",
                          control = list(maxit = maxit)),
                     loss_args))
    list(value = fit$value, par = fit$par,
         seconds = proc.time()[["elapsed"]] - t0)
  }

  restart_log <- data.frame(restart = seq_len(restarts),
                            value = NA_real_, seconds = NA_real_)
  best_theta <- theta0
  best_val <- Inf

  if (parallel) {
    if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    res <- parallel::parLapply(cl, starts, run_one)
    for (r in seq_len(restarts)) {
      restart_log$value[r]   <- res[[r]]$value
      restart_log$seconds[r] <- res[[r]]$seconds
      if (res[[r]]$value < best_val) {
        best_val   <- res[[r]]$value
        best_theta <- res[[r]]$par
      }
    }
  } else {
    n_eval <- 0L
    for (r in seq_len(restarts)) {
      res_r <- run_one(starts[[r]])
      restart_log$value[r]   <- res_r$value
      restart_log$seconds[r] <- res_r$seconds
      n_eval <- r
      if (res_r$value < best_val) {
        best_val   <- res_r$value
        best_theta <- res_r$par
        if (verbose) {
          cat(sprintf("New best loss = %.3e (restart %d, %.2fs)\n",
                      best_val, r, res_r$seconds))
        }
        if (best_val < early_tol) break
      } else if (verbose && r %% 10L == 0L) {
        cat(sprintf("Restart %d done (loss = %.3e, %.2fs). Current best = %.3e\n",
                    r, res_r$value, res_r$seconds, best_val))
      }
    }
    restart_log <- restart_log[seq_len(n_eval), , drop = FALSE]
  }

  if (polish && is.finite(best_val) && best_val < 1e9) {
    fit_p <- tryCatch(
      do.call(stats::optim,
              c(list(par = best_theta, fn = loss_varK,
                     method = "BFGS",
                     control = list(maxit = 2000L)),
                loss_args)),
      error = function(e) NULL)
    if (!is.null(fit_p) && is.finite(fit_p$value) && fit_p$value < best_val) {
      if (verbose) {
        cat(sprintf("BFGS polish improved loss: %.3e -> %.3e\n",
                    best_val, fit_p$value))
      }
      best_val   <- fit_p$value
      best_theta <- fit_p$par
    }
  }

  par_best <- unpack_theta(best_theta, K, alpha11_1, alpha22_1, var_u, var_e)
  mats <- make_mats_varK(par_best$ax, par_best$ay, par_best$b21, par_best$b12,
                         var_u = var_u, var_e = var_e, cov_ue = par_best$cov_ue)
  P <- stationary_cov(mats$F, mats$Q, method = "kron")
  implied <- implied_clpm_profile(mats$F, P, delta_max)

  sec_total <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (verbose) {
    cat(sprintf("\nCalibration runtime: %.2f seconds (%.2f minutes)\n",
                sec_total, sec_total / 60))
    cat(sprintf("Restarts evaluated: %d\n", nrow(restart_log)))
  }

  out <- list(
    theta = best_theta,
    loss = best_val,
    implied = implied,
    params = par_best,
    mats = mats,
    P = P,
    targets_b21 = targets_b21,
    targets_b12 = targets_b12,
    settings = list(K = K, delta_max = delta_max,
                    alpha11_1 = alpha11_1, alpha22_1 = alpha22_1,
                    var_u = var_u, var_e = var_e,
                    stab_thresh = stab_thresh, restarts = restarts,
                    maxit = maxit, jitter_sd = jitter_sd,
                    early_tol = early_tol, seed = seed,
                    polish = polish, parallel = parallel),
    runtime_seconds = sec_total,
    restart_log = restart_log
  )
  class(out) <- "varK_calibration"
  out
}

#' Print method for varK calibration objects
#'
#' @param x An object of class \code{"varK_calibration"}.
#' @param digits Number of digits for the implied coefficient table.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @export
print.varK_calibration <- function(x, digits = 4, ...) {
  s <- x$settings
  cat(sprintf("varK calibration (K = %d, delta_max = %d)\n", s$K, s$delta_max))
  cat(sprintf("Fixed lag-1 AR of the DGP: alpha11_1 = %g, alpha22_1 = %g\n",
              s$alpha11_1, s$alpha22_1))
  cat(sprintf("Best loss: %.3e | runtime: %.1f s | restarts evaluated: %d\n",
              x$loss, x$runtime_seconds, nrow(x$restart_log)))
  dev21 <- x$implied["b21", ] - x$targets_b21
  dev12 <- x$implied["b12", ] - x$targets_b12
  cat(sprintf("Max |implied - target|: b21 = %.2e, b12 = %.2e\n",
              max(abs(dev21)), max(abs(dev12))))
  cat("\nImplied population CLPM coefficients (rows b11, b21, b22, b12):\n")
  print(round(x$implied, digits))
  invisible(x)
}
