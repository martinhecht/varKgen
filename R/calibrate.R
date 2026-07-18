#' Unpack the parameter vector of the calibration problem
#'
#' Internal helper. The parameter vector \code{theta} has length
#' \code{4K - 2} and layout
#' \code{(ax[2..K], ay[2..K], b21[1..K], b12[1..K])}; the lag-1
#' autoregressive coefficients are fixed at \code{alpha11_1} and
#' \code{alpha22_1} and are not free parameters. The process error
#' correlation is no longer a free parameter: the full 2 x 2 covariance
#' (process or process error, depending on \code{var_type}) is set by the
#' user in \code{\link{calibrate_varK}}.
#'
#' @param theta Numeric parameter vector of length \code{4K - 2}.
#' @param K Number of lags.
#' @param alpha11_1,alpha22_1 Fixed lag-1 autoregressive coefficients.
#'
#' @return List with elements \code{ax}, \code{ay}, \code{b21}, \code{b12}.
#' @keywords internal
unpack_theta <- function(theta, K, alpha11_1, alpha22_1) {
  stopifnot(length(theta) == 4L * K - 2L)
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
  b12 <- theta[idx:(idx + K - 1L)]

  list(ax = ax, ay = ay, b21 = b21, b12 = b12)
}

#' Calibration loss
#'
#' Sum of squared deviations between the implied population cross-lagged
#' coefficients \code{b21(delta)}, \code{b12(delta)} for
#' \code{delta = 1, ..., delta_max} and the target profiles. Unstable
#' parameter values receive a graded penalty; in \code{var_type = "process"}
#' mode, parameter values for which the derived process error covariance is
#' not positive semi-definite (inadmissible target process covariance)
#' receive a graded penalty as well.
#'
#' @details
#' With \code{var_type = "process"}, the stationary process covariance block
#' is fixed at \code{(var1, var2, cov12)} and the process error covariance is
#' derived in every evaluation via the internal linear solver
#' \code{solve_process_error()} (one LU factorization with three right-hand
#' sides; the stationary covariance follows by linearity without a second
#' Lyapunov solve). With \code{var_type = "process_error"}, \code{(var1,
#' var2, cov12)} specify the process error covariance directly and the
#' stationary covariance is computed with \code{\link{stationary_cov}}.
#'
#' @param theta Parameter vector, see \code{\link{unpack_theta}}.
#' @param K Number of lags.
#' @param delta_max Largest interval; \code{length(targets_b21)} must equal
#'   \code{delta_max}.
#' @param alpha11_1,alpha22_1 Fixed lag-1 autoregressive coefficients.
#' @param targets_b21,targets_b12 Target profiles of length \code{delta_max}.
#' @param var1,var2,cov12 Covariance settings, interpreted according to
#'   \code{var_type}.
#' @param var_type Either \code{"process"} or \code{"process_error"}.
#' @param stab_thresh Stability threshold for the spectral radius of the
#'   companion matrix.
#'
#' @return A single non-negative loss value.
#' @keywords internal
loss_varK <- function(theta, K, delta_max, alpha11_1, alpha22_1,
                      targets_b21, targets_b12,
                      var1 = 1, var2 = 1, cov12 = 0,
                      var_type = "process", stab_thresh = 0.998) {
  if (!all(is.finite(theta))) return(1e12)

  par <- unpack_theta(theta, K, alpha11_1, alpha22_1)
  Fmat <- companion_F(par$ax, par$ay, par$b21, par$b12)

  eig <- eigen(Fmat, only.values = TRUE)$values
  if (!all(is.finite(eig))) return(1e12)
  maxmod <- max(Mod(eig))
  if (maxmod >= stab_thresh) return(1e9 + 1e7 * (maxmod - stab_thresh))

  if (var_type == "process") {
    pe <- solve_process_error(Fmat, var1, var2, cov12)
    if (is.null(pe)) return(1e9)
    tol_pe <- 1e-10 * max(1, abs(pe$Q2[1L, 1L]) + abs(pe$Q2[2L, 2L]))
    if (pe$min_eig < -tol_pe) return(1e9 + 1e7 * (-pe$min_eig))
    P <- pe$P
  } else {
    dimS <- nrow(Fmat)
    Qmat <- matrix(0, dimS, dimS)
    Qmat[1L, 1L] <- var1
    Qmat[2L, 2L] <- var2
    Qmat[1L, 2L] <- Qmat[2L, 1L] <- cov12
    P <- stationary_cov(Fmat, Qmat, method = "kron")
    if (is.null(P)) return(1e9)
  }

  prof <- implied_clpm_profile(Fmat, P, delta_max)
  if (!all(is.finite(prof))) return(1e9)

  sum((prof["b21", ] - targets_b21)^2) + sum((prof["b12", ] - targets_b12)^2)
}

#' Default start values for the calibration
#'
#' Adapts the start-value heuristic of the original script to the
#' \code{4K - 2} parameter layout: small positive values for the higher-lag
#' autoregressive coefficients, the first cross-lagged target as start value
#' for the lag-1 cross effects, and small positive values for the remaining
#' cross-lag coefficients. There is no start value for a process error
#' correlation because that parameter no longer exists (the full covariance
#' is set by the user).
#'
#' @param K Number of lags.
#' @param targets_b21,targets_b12 Target profiles (only the first element of
#'   each is used).
#'
#' @return Numeric vector of length \code{4K - 2} with layout
#'   \code{(ax[2..K], ay[2..K], b21[1..K], b12[1..K])}.
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
    c(targets_b12[1L], rep(0.03, K - 1L))
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
#' The contemporaneous 2 x 2 covariance is always fully specified by
#' \code{var1}, \code{var2} and \code{cov12}; its meaning depends on
#' \code{var_type}:
#' \describe{
#'   \item{\code{"process"} (default)}{\code{var1}, \code{var2},
#'     \code{cov12} are the target stationary process covariance, i.e. the
#'     top-left 2 x 2 block of the stationary state covariance \code{P}. In
#'     every loss evaluation the required process error covariance is derived
#'     exactly by an internal linear solver (also for unequal variances).
#'     With \code{var1 = var2 = 1} the process variances are one, so
#'     standardized and unstandardized cross-lagged coefficients coincide and
#'     target profiles can be specified in the standardized metric. Not every
#'     target process covariance is admissible for given dynamics; parameter
#'     values whose derived process error covariance is not positive
#'     semi-definite receive a graded penalty, and the optimizer avoids such
#'     regions.}
#'   \item{\code{"process_error"}}{\code{var1}, \code{var2}, \code{cov12}
#'     are the process error covariance itself; the stationary process
#'     covariance is a result and is returned as \code{process_cov}.}
#' }
#'
#' Compared to earlier versions, the process error correlation is no longer a
#' free calibration parameter (it is fully determined by the covariance
#' settings), so the number of free parameters is \code{4K - 2}: the
#' higher-lag autoregressions \code{ax[2..K]}, \code{ay[2..K]} and the
#' cross-lag coefficients \code{b21[1..K]}, \code{b12[1..K]}. The number of
#' target values is \code{2 * delta_max}. If \code{2 * delta_max > 4K - 2},
#' the targets can in general only be matched approximately and the
#' attainable loss has a positive floor; a message points this out.
#'
#' For a given \code{seed}, the jittered start values are drawn upfront
#' (Nelder-Mead itself does not consume random numbers). If the best loss is
#' at the penalty level (>= 1e9), no stable and admissible solution was found
#' and a warning is issued.
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
#' @param var1,var2 Variances (both default 1), interpreted according to
#'   \code{var_type}; must be positive.
#' @param cov12 Covariance (default 0), interpreted according to
#'   \code{var_type}; must satisfy \code{cov12^2 < var1 * var2}.
#' @param var_type Either \code{"process"} (default; \code{var1},
#'   \code{var2}, \code{cov12} fix the stationary process covariance) or
#'   \code{"process_error"} (they fix the process error covariance). See
#'   Details.
#' @param theta0 Optional start vector of length \code{4K - 2}; defaults to
#'   \code{\link{default_theta0}}.
#' @param stab_thresh Stability threshold for the spectral radius (default
#'   0.998).
#' @param restarts Number of jittered Nelder-Mead restarts (default 120).
#' @param maxit Maximum Nelder-Mead iterations per restart (default 25000).
#' @param jitter_sd Standard deviation of the normal jitter added to
#'   \code{theta0} for each restart (default 0.10).
#' @param early_tol Early-stopping tolerance on the loss (default 1e-10; set
#'   to 0 to disable early stopping, e.g. for optimizer studies).
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
#'   (\code{ax}, \code{ay}, \code{b21}, \code{b12} plus the process error
#'   covariance actually used: \code{pe_var1}, \code{pe_var2},
#'   \code{pe_cov12}), \code{mats} (companion matrices \code{F}, \code{Q}),
#'   \code{P} (stationary state covariance), \code{process_cov} (achieved
#'   stationary process covariance, the top-left 2 x 2 block of \code{P}),
#'   \code{targets_b21}, \code{targets_b12}, \code{settings},
#'   \code{runtime_seconds} and \code{restart_log}.
#'
#' @examples
#' \donttest{
#' targets_b21 <- c(0.30, 0.25, 0.15, 0.12, 0.10)
#' targets_b12 <- c(0.15, 0.20, 0.10, 0.10, 0.08)
#' sol <- calibrate_varK(targets_b21, targets_b12, K = 3,
#'                       alpha11_1 = 0.5, alpha22_1 = 0.5,
#'                       var1 = 1, var2 = 1, cov12 = 0,
#'                       var_type = "process",
#'                       restarts = 10, maxit = 5000, verbose = FALSE)
#' print(sol)
#' sol$process_cov   # equals diag(2) up to numerical precision
#' }
#'
#' @export
calibrate_varK <- function(targets_b21, targets_b12, K,
                           alpha11_1, alpha22_1,
                           var1 = 1, var2 = 1, cov12 = 0,
                           var_type = c("process", "process_error"),
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

  var_type <- match.arg(var_type)

  stopifnot(length(targets_b21) == length(targets_b12),
            length(targets_b21) >= 1L,
            all(is.finite(targets_b21)), all(is.finite(targets_b12)),
            length(K) == 1L, K >= 1L, K == round(K),
            is.finite(alpha11_1), is.finite(alpha22_1),
            length(var1) == 1L, is.finite(var1), var1 > 0,
            length(var2) == 1L, is.finite(var2), var2 > 0,
            length(cov12) == 1L, is.finite(cov12),
            cov12^2 < var1 * var2,
            stab_thresh > 0, stab_thresh < 1,
            restarts >= 1L, maxit >= 1L, jitter_sd >= 0)
  delta_max <- length(targets_b21)
  n_par <- 4L * K - 2L

  if (is.null(theta0)) theta0 <- default_theta0(K, targets_b21, targets_b12)
  stopifnot(length(theta0) == n_par, all(is.finite(theta0)))

  if (2L * delta_max > n_par && verbose) {
    message(sprintf(paste0(
      "Note: %d target values vs. %d free parameters (4K - 2). ",
      "Targets can in general only be matched approximately ",
      "(attainable loss may have a positive floor)."),
      2L * delta_max, n_par))
  }

  loss_args <- list(K = K, delta_max = delta_max,
                    alpha11_1 = alpha11_1, alpha22_1 = alpha22_1,
                    targets_b21 = targets_b21, targets_b12 = targets_b12,
                    var1 = var1, var2 = var2, cov12 = cov12,
                    var_type = var_type, stab_thresh = stab_thresh)

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

  # ---- reconstruct the solution ------------------------------------------
  par_best <- unpack_theta(best_theta, K, alpha11_1, alpha22_1)
  Fmat <- companion_F(par_best$ax, par_best$ay, par_best$b21, par_best$b12)
  dimS <- nrow(Fmat)

  if (var_type == "process") {
    pe <- solve_process_error(Fmat, var1, var2, cov12)
    if (is.null(pe)) {
      warning(paste0("The returned parameter vector is not stable; the ",
                     "process error covariance and the stationary covariance ",
                     "are unavailable (NA)."))
      Q2 <- matrix(NA_real_, 2L, 2L)
      P <- NULL
    } else {
      Q2 <- pe$Q2
      P <- pe$P
      tol_pe <- 1e-10 * max(1, abs(Q2[1L, 1L]) + abs(Q2[2L, 2L]))
      if (pe$min_eig < -tol_pe) {
        warning(paste0("The derived process error covariance at the returned ",
                       "solution is not positive semi-definite: the target ",
                       "process covariance is not admissible for these ",
                       "dynamics."))
      }
    }
  } else {
    Q2 <- matrix(c(var1, cov12, cov12, var2), 2L, 2L)
    Qmat0 <- matrix(0, dimS, dimS)
    Qmat0[1:2, 1:2] <- Q2
    P <- stationary_cov(Fmat, Qmat0, method = "kron")
  }

  Qmat <- matrix(0, dimS, dimS)
  Qmat[1:2, 1:2] <- Q2
  mats <- list(F = Fmat, Q = Qmat)
  implied <- implied_clpm_profile(Fmat, P, delta_max)
  process_cov <- if (is.null(P)) matrix(NA_real_, 2L, 2L) else P[1:2, 1:2]

  if (is.finite(best_val) && best_val >= 1e9) {
    warning(paste0("No admissible solution found: the best loss is at the ",
                   "penalty level. The returned object should not be used ",
                   "as a data-generating process."))
  }

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
    params = list(ax = par_best$ax, ay = par_best$ay,
                  b21 = par_best$b21, b12 = par_best$b12,
                  pe_var1 = Q2[1L, 1L], pe_var2 = Q2[2L, 2L],
                  pe_cov12 = Q2[1L, 2L]),
    mats = mats,
    P = P,
    process_cov = process_cov,
    targets_b21 = targets_b21,
    targets_b12 = targets_b12,
    settings = list(K = K, delta_max = delta_max,
                    alpha11_1 = alpha11_1, alpha22_1 = alpha22_1,
                    var1 = var1, var2 = var2, cov12 = cov12,
                    var_type = var_type,
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
  cat(sprintf("varK calibration (K = %d, delta_max = %d, var_type = \"%s\")\n",
              s$K, s$delta_max, s$var_type))
  cat(sprintf("Fixed lag-1 AR of the DGP: alpha11_1 = %g, alpha22_1 = %g\n",
              s$alpha11_1, s$alpha22_1))
  if (identical(s$var_type, "process")) {
    cat(sprintf("Target process covariance: var1 = %g, var2 = %g, cov12 = %g\n",
                s$var1, s$var2, s$cov12))
    cat(sprintf("Derived process error covariance: pe_var1 = %.4g, pe_var2 = %.4g, pe_cov12 = %.4g\n",
                x$params$pe_var1, x$params$pe_var2, x$params$pe_cov12))
  } else {
    cat(sprintf("Process error covariance (set): pe_var1 = %g, pe_var2 = %g, pe_cov12 = %g\n",
                s$var1, s$var2, s$cov12))
    cat(sprintf("Resulting process covariance: var1 = %.4g, var2 = %.4g, cov12 = %.4g\n",
                x$process_cov[1L, 1L], x$process_cov[2L, 2L],
                x$process_cov[1L, 2L]))
  }
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
