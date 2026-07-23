#' Unpack the parameter vector of the calibration problem
#'
#' Internal helper. The parameter vector \code{theta} has length \code{4K}
#' and layout \code{(ax[1..K], ay[1..K], b21[1..K], b12[1..K])}. Since
#' version 0.3.0 every element of every VAR coefficient matrix is a free
#' parameter, including the lag-1 autoregressive diagonal elements
#' \code{ax[1]} and \code{ay[1]}. The process error correlation is not a
#' free parameter: the full 2 x 2 covariance (process or process error,
#' depending on \code{var_type}) is set by the user in
#' \code{\link{calibrate_varK}}.
#'
#' @param theta Numeric parameter vector of length \code{4K}.
#' @param K Number of lags.
#'
#' @return List with the unnamed elements \code{ax}, \code{ay},
#'   \code{b21}, \code{b12}, each of length \code{K}.
#' @keywords internal
unpack_theta <- function(theta, K) {
  stopifnot(length(theta) == 4L * K)
  idx <- 1L

  ax  <- theta[idx:(idx + K - 1L)]; idx <- idx + K
  ay  <- theta[idx:(idx + K - 1L)]; idx <- idx + K
  b21 <- theta[idx:(idx + K - 1L)]; idx <- idx + K
  b12 <- theta[idx:(idx + K - 1L)]

  # the coefficient blocks are indexed by lag, so they are returned without
  # names; canonical names live on `theta` itself (and on `restart_theta`)
  list(ax = unname(ax), ay = unname(ay),
       b21 = unname(b21), b12 = unname(b12))
}

#' Calibration loss
#'
#' Sum of squared deviations between the implied population coefficients and
#' the targets. Each individual target coefficient contributes one squared
#' deviation with unit weight; no separate block weighting is applied. The
#' targets are the cross-lagged coefficients \code{b21(Delta_j)} and
#' \code{b12(Delta_j)} for all target intervals \code{Delta_j = j},
#' \code{j = 1, ..., delta_max}, plus the two autoregressive projection
#' effects \code{b11(Delta_1)} and \code{b22(Delta_1)} at the first target
#' interval \code{Delta_1 = 1} time unit.
#' Unstable parameter values receive a graded penalty; in
#' \code{var_type = "process"} mode, parameter values for which the derived
#' process error covariance is not positive semi-definite (inadmissible
#' target process covariance) receive a graded penalty as well.
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
#' @param targets_b21,targets_b12 Target profiles of length \code{delta_max}.
#' @param target_b11_1,target_b22_1 Target population autoregressive
#'   projection effects at the first target interval \code{Delta_1 = 1} time
#'   unit, i.e. targets for \code{b11(Delta_1)} and \code{b22(Delta_1)}.
#' @param var1,var2,cov12 Covariance settings, interpreted according to
#'   \code{var_type}.
#' @param var_type Either \code{"process"} or \code{"process_error"}.
#' @param stab_thresh Stability threshold for the spectral radius of the
#'   companion matrix.
#'
#' @return A single non-negative loss value.
#' @keywords internal
loss_varK <- function(theta, K, delta_max,
                      targets_b21, targets_b12,
                      target_b11_1, target_b22_1,
                      var1 = 1, var2 = 1, cov12 = 0,
                      var_type = "process", stab_thresh = 0.998) {
  if (!all(is.finite(theta))) return(1e12)

  par <- unpack_theta(theta, K)
  Fmat <- companion_F(par$ax, par$ay, par$b21, par$b12)

  eig <- eigen(Fmat, only.values = TRUE)$values
  if (!all(is.finite(eig))) return(1e12)
  maxmod <- max(Mod(eig))
  if (maxmod >= stab_thresh) return(1e9 + 1e7 * (maxmod - stab_thresh))

  if (var_type == "process") {
    # stability was just established above, so the solver must not repeat
    # the eigendecomposition
    pe <- solve_process_error(Fmat, var1, var2, cov12,
                              check_stability = FALSE)
    if (is.null(pe)) return(1e9)
    tol_pe <- 1e-10 * max(.Machine$double.eps,
                          abs(pe$Q2[1L, 1L]) + abs(pe$Q2[2L, 2L]))
    if (pe$min_eig < -tol_pe) return(1e9 + 1e7 * (-pe$min_eig))
    P <- pe$P
  } else {
    dimS <- nrow(Fmat)
    Qmat <- matrix(0, dimS, dimS)
    Qmat[1L, 1L] <- var1
    Qmat[2L, 2L] <- var2
    Qmat[1L, 2L] <- Qmat[2L, 1L] <- cov12
    P <- stationary_cov(Fmat, Qmat, method = "kron",
                        check_stability = FALSE)
    if (is.null(P)) return(1e9)
  }

  # Single canonical computation: implied_clpm_profile() evaluates the fixed
  # target grid Delta_j = j, so column 1 is by construction the first target
  # interval Delta_1 = 1 time unit. No second projection routine is used.
  prof <- implied_clpm_profile(Fmat, P, delta_max)
  if (!all(is.finite(prof))) return(1e9)

  sum((prof["b21", ] - targets_b21)^2) +
    sum((prof["b12", ] - targets_b12)^2) +
    (prof["b11", 1L] - target_b11_1)^2 +
    (prof["b22", 1L] - target_b22_1)^2
}

#' Default start values for the calibration
#'
#' Adapts the start-value heuristic to the \code{4K} parameter layout: the
#' autoregressive target values as start values for the lag-1 autoregressive
#' coefficients, small positive values for the higher-lag autoregressive
#' coefficients, the first cross-lagged target as start value for the lag-1
#' cross effects, and small positive values for the remaining cross-lag
#' coefficients. There is no start value for a process error correlation
#' because that parameter does not exist (the full covariance is set by the
#' user).
#'
#' Using \code{target_b11_1} and \code{target_b22_1} as start values for
#' \code{ax[1]} and \code{ay[1]} is a heuristic: the two quantities coincide
#' only for \code{K = 1}, but for higher \code{K} the target value is still a
#' reasonable starting point.
#'
#' @param K Number of lags.
#' @param targets_b21,targets_b12 Target profiles (only the first element of
#'   each is used).
#' @param target_b11_1,target_b22_1 Targets for \code{b11(Delta_1)} and
#'   \code{b22(Delta_1)} at the first target interval \code{Delta_1 = 1}
#'   time unit.
#'
#' @return Numeric vector of length \code{4K} with layout
#'   \code{(ax[1..K], ay[1..K], b21[1..K], b12[1..K])}.
#'
#' @examples
#' default_theta0(K = 5,
#'                targets_b21 = c(0.30, 0.25, 0.15, 0.12, 0.10),
#'                targets_b12 = c(0.15, 0.20, 0.10, 0.10, 0.08),
#'                target_b11_1 = 0.5, target_b22_1 = 0.5)
#'
#' @export
default_theta0 <- function(K, targets_b21, targets_b12,
                           target_b11_1, target_b22_1) {
  stopifnot(length(K) == 1L, is.finite(K), K >= 1L, K == round(K),
            length(targets_b21) >= 1L, all(is.finite(targets_b21)),
            length(targets_b12) >= 1L, all(is.finite(targets_b12)),
            length(target_b11_1) == 1L, is.finite(target_b11_1),
            length(target_b22_1) == 1L, is.finite(target_b22_1))
  c(
    c(target_b11_1, if (K > 1L) rep(0.05, K - 1L) else numeric(0)),
    c(target_b22_1, if (K > 1L) rep(0.05, K - 1L) else numeric(0)),
    c(targets_b21[1L], rep(0.03, K - 1L)),
    c(targets_b12[1L], rep(0.03, K - 1L))
  )
}

#' Run a single optimizer restart
#'
#' Internal helper. Wraps one \code{stats::optim()} Nelder-Mead call so that
#' neither an R error nor a formally returned but non-finite result can abort
#' a calibration or be mistaken for a completed restart.
#'
#' @param theta_start Start vector for this restart.
#' @param fn Objective function.
#' @param maxit Maximum Nelder-Mead iterations.
#' @param fn_args List of further arguments passed to \code{fn}.
#'
#' @return A list with \code{value}, \code{par}, \code{seconds},
#'   \code{convergence}, \code{fn_evals}, \code{message} and
#'   \code{status}. \code{status} is \code{"completed"} if the optimizer
#'   returned a finite value and a finite parameter vector of the expected
#'   length, \code{"error"} if the call raised an R error, and
#'   \code{"invalid"} if it returned a non-finite or malformed result. For
#'   the latter two, \code{value} is \code{Inf} and \code{par} is all
#'   \code{NA}, so such restarts cannot become the best solution and are
#'   excluded from hit and near-best analyses.
#' @keywords internal
run_restart <- function(theta_start, fn, maxit, fn_args = list()) {
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    do.call(stats::optim,
            c(list(par = theta_start, fn = fn, method = "Nelder-Mead",
                   control = list(maxit = maxit)),
              fn_args)),
    error = function(e) e)
  secs <- proc.time()[["elapsed"]] - t0
  failed <- function(status, msg) {
    list(value = Inf, par = rep(NA_real_, length(theta_start)),
         seconds = secs, convergence = NA_integer_, fn_evals = NA_integer_,
         message = msg, status = status)
  }
  if (inherits(fit, "error")) {
    return(failed("error", conditionMessage(fit)))
  }
  ok <- is.list(fit) &&
    length(fit$value) == 1L && is.finite(fit$value) &&
    length(fit$par) == length(theta_start) && all(is.finite(fit$par))
  if (!ok) {
    return(failed("invalid",
                  "optimizer returned a non-finite or malformed result"))
  }
  list(value = fit$value, par = fit$par, seconds = secs,
       convergence = fit$convergence,
       fn_evals = unname(fit$counts[["function"]]),
       message = if (is.null(fit$message)) NA_character_ else
         as.character(fit$message),
       status = "completed")
}

#' Calibrate a VAR(K) data-generating process to target coefficient profiles
#'
#' Searches, via multi-start Nelder-Mead (optionally followed by a BFGS
#' polishing step), for VAR(K) parameters such that the implied population
#' CLPM coefficients match the supplied targets. Two groups of targets are
#' matched jointly: the cross-lagged coefficients \code{b21(Delta_j)} and
#' \code{b12(Delta_j)} across the fixed target intervals
#' \code{Delta_j = j}, \code{j = 1, ..., delta_max}, and the two
#' autoregressive projection effects \code{b11(Delta_1)} and
#' \code{b22(Delta_1)} at the first target interval.
#'
#' @details
#' Since version 0.3.0 all elements of all VAR coefficient matrices are
#' freely optimized, including the lag-1 autoregressive diagonal elements.
#' The former arguments \code{alpha11_1} and \code{alpha22_1}, which fixed
#' those two DGP coefficients, have been removed and are not accepted as
#' aliases: they denoted a conceptually different quantity. The number of
#' free parameters is \code{4K}: \code{ax[1..K]}, \code{ay[1..K]},
#' \code{b21[1..K]}, \code{b12[1..K]}. The number of target values is
#' \code{2 * delta_max + 2}; each individual target coefficient contributes
#' one squared deviation with unit weight, and no separate block weighting
#' is applied.
#'
#' The target intervals are fixed:
#' \code{Delta = (Delta_1, ..., Delta_delta_max) = (1, ..., delta_max)}, so
#' \code{Delta_j = j} and in particular \code{Delta_1 = 1} time unit. The
#' subscript denotes the position in this fixed grid; a freely chosen
#' interval vector is not supported. Accordingly \code{targets_b21[j]} is
#' the target for \code{b21(Delta_j)}, and \code{target_b11_1},
#' \code{target_b22_1} are the targets for \code{b11(Delta_1)} and
#' \code{b22(Delta_1)}.
#'
#' \code{target_b11_1} and \code{target_b22_1} are population-level
#' projection effects, that is the diagonal elements of the population
#' projection matrix \code{B(Delta_1)}, each controlling for the respective
#' other process. For \code{K = 1}, \code{B(Delta_1)} equals the VAR(1)
#' coefficient matrix. For \code{K > 1}, the projection effects generally
#' differ from the corresponding elements of the first VAR coefficient
#' matrix because they depend on the complete VAR(K) process; equality can
#' still occur in special cases, for instance when all higher VAR lag
#' matrices are zero. They come from the same canonical routine as
#' the rest of the profile, \code{\link{implied_clpm_profile}}, whose first
#' column is \code{Delta_1}.
#'
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
#'     standardized and unstandardized coefficients coincide and targets can
#'     be specified in the standardized metric. Not every target process
#'     covariance is admissible for given dynamics; parameter values whose
#'     derived process error covariance is not positive semi-definite receive
#'     a graded penalty, and the optimizer avoids such regions.}
#'   \item{\code{"process_error"}}{\code{var1}, \code{var2}, \code{cov12}
#'     are the process error covariance itself; the stationary process
#'     covariance is a result and is returned as \code{process_cov}.}
#' }
#'
#' If \code{2 * delta_max + 2 > 4K}, the targets can in general only be
#' matched approximately and the attainable loss has a positive floor; a
#' message points this out.
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
#' @param target_b11_1 Target population autoregressive projection effect of
#'   \code{Y1} at time \code{t} on \code{Y1} at time
#'   \code{t + Delta_1}, where \code{Delta_1 = 1} time unit, controlling
#'   for \code{Y2} at time \code{t}. See Details.
#' @param target_b22_1 Target population autoregressive projection effect of
#'   \code{Y2} at time \code{t} on \code{Y2} at time
#'   \code{t + Delta_1}, where \code{Delta_1 = 1} time unit, controlling
#'   for \code{Y1} at time \code{t}. See Details.
#' @param var1,var2 Variances (both default 1), interpreted according to
#'   \code{var_type}; must be positive.
#' @param cov12 Covariance (default 0), interpreted according to
#'   \code{var_type}; must satisfy \code{cov12^2 < var1 * var2}.
#' @param var_type Either \code{"process"} (default; \code{var1},
#'   \code{var2}, \code{cov12} fix the stationary process covariance) or
#'   \code{"process_error"} (they fix the process error covariance). See
#'   Details.
#' @param theta0 Optional start vector of length \code{4K}; defaults to
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
#' @param ... Must be empty. Any additional argument raises an error; the
#'   removed arguments \code{alpha11_1} and \code{alpha22_1} are
#'   intercepted with a specific message.
#'
#' @importFrom stats rnorm
#'
#' @return An object of class \code{"varK_calibration"}: a list with elements
#'   \code{theta} (best parameter vector), \code{loss}, \code{implied}
#'   (4 x delta_max matrix of implied population coefficients, rows
#'   \code{b11}, \code{b21}, \code{b22}, \code{b12}), \code{implied_b11_1}
#'   and \code{implied_b22_1} (achieved autoregressive projection effects at
#'   \code{Delta_1}), \code{error_b11_1} and \code{error_b22_1} (implied
#'   minus target at \code{Delta_1}), \code{params} (\code{ax}, \code{ay},
#'   \code{b21}, \code{b12} plus the process error covariance actually used:
#'   \code{pe_var1}, \code{pe_var2}, \code{pe_cov12}), \code{mats}
#'   (companion matrices \code{F}, \code{Q}), \code{P} (stationary state
#'   covariance), \code{process_cov} (achieved stationary process covariance,
#'   the top-left 2 x 2 block of \code{P}), \code{targets_b21},
#'   \code{targets_b12}, \code{target_b11_1}, \code{target_b22_1},
#'   \code{settings}, \code{runtime_seconds}, \code{restart_log} and
#'   \code{restart_theta}.
#'
#'   \code{theta} and \code{loss} always describe the solution finally used.
#'   The optimizer stages are kept separately: \code{theta_nm} and
#'   \code{loss_nm} are the best pre-polish Nelder-Mead solution, on which
#'   hits and near-best variability are defined; \code{theta_polished},
#'   \code{loss_polished}, \code{polish_convergence} and
#'   \code{polish_improved} describe the BFGS step (\code{NULL} / \code{NA}
#'   if no polish was run).
#'
#'   \code{restart_log} has one row per evaluated restart with
#'   \code{restart} (the index), \code{value}, \code{seconds},
#'   \code{convergence}, \code{fn_evals}, \code{message} and
#'   \code{status}. \code{status} is \code{"completed"} if the optimizer
#'   returned a finite value and a finite parameter vector, \code{"error"}
#'   if the call raised an R error, and \code{"invalid"} if it returned a
#'   non-finite or malformed result. Note that \code{"completed"} does not
#'   imply convergence: a non-zero \code{convergence} code remains possible
#'   and must be judged separately. For failing restarts \code{value} is
#'   \code{Inf} and the corresponding row of \code{restart_theta} is all
#'   \code{NA}; such restarts are excluded when the best solution is
#'   selected. If no restart is usable, the function stops with an error.
#'   \code{restart_theta} is the matching matrix of full
#'   pre-polish Nelder-Mead endpoints (one row per restart, columns named
#'   \code{ax1..axK}, \code{ay1..ayK}, \code{b211..b21K},
#'   \code{b121..b12K}). Validity, spectral radius, smallest eigenvalue of
#'   the process error covariance, implied autoregressive patterns and all
#'   near-best variability measures can be reconstructed from these
#'   endpoints.
#'
#' @examples
#' \donttest{
#' targets_b21 <- c(0.30, 0.25, 0.15, 0.12, 0.10)
#' targets_b12 <- c(0.15, 0.20, 0.10, 0.10, 0.08)
#' sol <- calibrate_varK(targets_b21, targets_b12, K = 3,
#'                       target_b11_1 = 0.5, target_b22_1 = 0.5,
#'                       var1 = 1, var2 = 1, cov12 = 0,
#'                       var_type = "process",
#'                       restarts = 10, maxit = 5000, verbose = FALSE)
#' print(sol)
#' sol$process_cov   # equals diag(2) up to numerical precision
#' }
#'
#' @export
calibrate_varK <- function(targets_b21, targets_b12, K,
                           target_b11_1, target_b22_1,
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
                           verbose = TRUE,
                           ...) {

  dots <- list(...)
  if (length(dots)) {
    removed <- intersect(names(dots), c("alpha11_1", "alpha22_1"))
    if (length(removed)) {
      stop("`alpha11_1` and `alpha22_1` are no longer supported. The VAR ",
           "coefficients are now freely optimized. Use `target_b11_1` and ",
           "`target_b22_1` to specify the population-level autoregressive ",
           "projection effects at the first target interval Delta_1 = 1 ",
           "time unit.",
           call. = FALSE)
    }
    nm <- names(dots)
    nm <- if (is.null(nm)) character(0) else nm[nzchar(nm)]
    stop("unused argument(s): ",
         if (length(nm)) paste(nm, collapse = ", ") else "unnamed",
         call. = FALSE)
  }

  var_type <- match.arg(var_type)

  stopifnot(length(targets_b21) == length(targets_b12),
            length(targets_b21) >= 1L,
            all(is.finite(targets_b21)), all(is.finite(targets_b12)),
            length(K) == 1L, is.numeric(K), !is.logical(K), is.finite(K),
            K >= 1L, K == round(K), K <= 1000L,
            length(target_b11_1) == 1L, is.finite(target_b11_1),
            length(target_b22_1) == 1L, is.finite(target_b22_1),
            length(var1) == 1L, is.finite(var1), var1 > 0,
            length(var2) == 1L, is.finite(var2), var2 > 0,
            length(cov12) == 1L, is.finite(cov12),
            cov12^2 < var1 * var2)
  num1 <- function(x) length(x) == 1L && is.numeric(x) && !is.logical(x) &&
    is.finite(x)
  int1 <- function(x) num1(x) && x == round(x) && abs(x) <= .Machine$integer.max
  stopifnot(num1(stab_thresh),
            stab_thresh > 0, stab_thresh < 1,
            int1(restarts), restarts >= 1L,
            int1(maxit), maxit >= 1L,
            num1(jitter_sd), jitter_sd >= 0,
            num1(early_tol), early_tol >= 0, early_tol < 1e9,
            int1(seed),
            length(polish) == 1L, is.logical(polish), !is.na(polish),
            length(parallel) == 1L, is.logical(parallel), !is.na(parallel),
            length(verbose) == 1L, is.logical(verbose), !is.na(verbose),
            is.null(cores) || (int1(cores) && cores >= 1L))
  # the values 1e9 / 1e12 are sentinels for invalid parameters; a valid SSE
  # must stay far below them, otherwise valid solutions could be classified
  # as penalty solutions
  target_scale <- sum(targets_b21^2) + sum(targets_b12^2) +
    target_b11_1^2 + target_b22_1^2
  if (target_scale >= 1e6) {
    stop("The targets are too large for the penalty sentinels used here ",
         "(sum of squared targets must stay below 1e6). Rescale the ",
         "problem instead.", call. = FALSE)
  }

  delta_max <- length(targets_b21)
  n_par <- 4L * K
  n_targets <- 2L * delta_max + 2L

  if (is.null(theta0)) {
    theta0 <- default_theta0(K, targets_b21, targets_b12,
                             target_b11_1, target_b22_1)
  }
  stopifnot(length(theta0) == n_par, all(is.finite(theta0)))

  if (n_targets > n_par && verbose) {
    message(sprintf(paste0(
      "Note: %d target values vs. %d free parameters (4K). ",
      "Targets can in general only be matched approximately ",
      "(attainable loss may have a positive floor)."),
      n_targets, n_par))
  }

  loss_args <- list(K = K, delta_max = delta_max,
                    targets_b21 = targets_b21, targets_b12 = targets_b12,
                    target_b11_1 = target_b11_1,
                    target_b22_1 = target_b22_1,
                    var1 = var1, var2 = var2, cov12 = cov12,
                    var_type = var_type, stab_thresh = stab_thresh)

  t_start <- Sys.time()

  # Draw all jittered starts upfront (same rnorm sequence as drawing them
  # one-by-one in a loop, since Nelder-Mead does not use the RNG).
  set.seed(seed)
  starts <- lapply(seq_len(restarts),
                   function(r) theta0 + stats::rnorm(n_par, sd = jitter_sd))

  run_one <- function(th) run_restart(th, loss_varK, maxit, loss_args)

  theta_names <- c(paste0("ax", seq_len(K)), paste0("ay", seq_len(K)),
                   paste0("b21", seq_len(K)), paste0("b12", seq_len(K)))
  restart_log <- data.frame(restart = seq_len(restarts),
                            value = NA_real_, seconds = NA_real_,
                            convergence = NA_integer_,
                            fn_evals = NA_integer_,
                            message = NA_character_,
                            status = NA_character_,
                            stringsAsFactors = FALSE)
  # full pre-polish Nelder-Mead endpoint of every restart; everything else
  # (validity, spectral radius, smallest eigenvalue of Q, near-best spread)
  # can be reconstructed from these vectors
  restart_theta <- matrix(NA_real_, nrow = restarts, ncol = n_par,
                          dimnames = list(NULL, theta_names))
  best_theta <- theta0
  best_val <- Inf

  cores_requested <- cores
  cores_used <- 1L
  if (parallel) {
    if (is.null(cores)) {
      detected <- suppressWarnings(parallel::detectCores())
      # detectCores() may return NA, and on a scheduled node it reflects the
      # hardware rather than the job allocation; set `cores` explicitly there
      cores <- if (is.finite(detected)) max(1L, detected - 1L) else 1L
    }
    cores <- min(as.integer(cores), restarts)
    cores_used <- cores
    cl <- parallel::makeCluster(cores)
    # on.exit is the failure fallback; on success the cluster is closed
    # immediately so that no workers idle through polish and reconstruction
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    # PSOCK workers start with the default library paths, so a package
    # installed in a project, check or user library would not be loadable
    # there; propagate the master's paths before dispatching
    parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())
    res <- parallel::parLapply(cl, starts, run_one)
    parallel::stopCluster(cl)
    for (r in seq_len(restarts)) {
      restart_log$value[r]       <- res[[r]]$value
      restart_log$seconds[r]     <- res[[r]]$seconds
      restart_log$convergence[r] <- res[[r]]$convergence
      restart_log$fn_evals[r]    <- res[[r]]$fn_evals
      restart_log$message[r]     <- res[[r]]$message
      restart_log$status[r]      <- res[[r]]$status
      restart_theta[r, ]         <- res[[r]]$par
      if (res[[r]]$value < best_val) {
        best_val   <- res[[r]]$value
        best_theta <- res[[r]]$par
      }
    }
  } else {
    n_eval <- 0L
    for (r in seq_len(restarts)) {
      res_r <- run_one(starts[[r]])
      restart_log$value[r]       <- res_r$value
      restart_log$seconds[r]     <- res_r$seconds
      restart_log$convergence[r] <- res_r$convergence
      restart_log$fn_evals[r]    <- res_r$fn_evals
      restart_log$message[r]     <- res_r$message
      restart_log$status[r]      <- res_r$status
      restart_theta[r, ]         <- res_r$par
      n_eval <- r
      if (res_r$value < best_val) {
        best_val   <- res_r$value
        best_theta <- res_r$par
        if (verbose) {
          cat(sprintf("New best loss = %.3e (restart %d, %.2fs)\n",
                      best_val, r, res_r$seconds))
        }
        if (best_val < early_tol && best_val < 1e9) break
      } else if (verbose && r %% 10L == 0L) {
        cat(sprintf("Restart %d done (loss = %.3e, %.2fs). Current best = %.3e\n",
                    r, res_r$value, res_r$seconds, best_val))
      }
    }
    restart_log <- restart_log[seq_len(n_eval), , drop = FALSE]
    restart_theta <- restart_theta[seq_len(n_eval), , drop = FALSE]
  }

  # determine the usable restarts explicitly instead of relying only on the
  # incremental update, and fail loudly if none is usable
  usable <- restart_log$status == "completed" &
    is.finite(restart_log$value) &
    apply(restart_theta, 1L, function(z) all(is.finite(z)))
  if (!any(usable)) {
    stop("All optimizer restarts failed; no calibration solution is ",
         "available. See the restart status and messages of the failing ",
         "run.", call. = FALSE)
  }
  i_best <- which(usable)[which.min(restart_log$value[usable])]
  best_val   <- restart_log$value[i_best]
  best_theta <- restart_theta[i_best, ]

  # the pre-polish Nelder-Mead optimum is kept separately: hits and near-best
  # variability are defined on pre-polish solutions
  theta_nm <- best_theta
  loss_nm  <- best_val
  theta_polished <- NULL
  loss_polished <- NA_real_
  polish_convergence <- NA_integer_
  polish_improved <- FALSE
  polish_attempted <- FALSE
  polish_status <- if (!polish) "not requested" else NA_character_
  polish_message <- NA_character_
  polish_fn_evals <- NA_integer_
  polish_seconds <- NA_real_

  if (polish) {
    if (!(is.finite(loss_nm) && loss_nm < 1e9)) {
      polish_status <- "skipped (penalty solution)"
    } else {
      polish_attempted <- TRUE
      tp0 <- proc.time()[["elapsed"]]
      fit_p <- tryCatch(
        do.call(stats::optim,
                c(list(par = theta_nm, fn = loss_varK,
                       method = "BFGS",
                       control = list(maxit = 2000L)),
                  loss_args)),
        error = function(e) e)
      polish_seconds <- proc.time()[["elapsed"]] - tp0
      if (inherits(fit_p, "error")) {
        polish_status  <- "error"
        polish_message <- conditionMessage(fit_p)
      } else if (!(length(fit_p$value) == 1L && is.finite(fit_p$value) &&
                   length(fit_p$par) == n_par && all(is.finite(fit_p$par)))) {
        polish_status  <- "invalid"
        polish_message <- "BFGS returned a non-finite or malformed result"
      } else {
        polish_status      <- "completed"
        theta_polished     <- fit_p$par
        loss_polished      <- fit_p$value
        polish_convergence <- fit_p$convergence
        polish_fn_evals    <- unname(fit_p$counts[["function"]])
        if (loss_polished < loss_nm) {
          polish_improved <- TRUE
          best_val   <- loss_polished
          best_theta <- theta_polished
          if (verbose) {
            cat(sprintf("BFGS polish improved loss: %.3e -> %.3e\n",
                        loss_nm, loss_polished))
          }
        }
      }
    }
  }

  # ---- reconstruct the solution ------------------------------------------
  par_best <- unpack_theta(best_theta, K)
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
      tol_pe <- 1e-10 * max(.Machine$double.eps,
                            abs(Q2[1L, 1L]) + abs(Q2[2L, 2L]))
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
  # same canonical source as the loss: column 1 is Delta_1 = 1 time unit
  implied_b11_1 <- unname(implied["b11", 1L])
  implied_b22_1 <- unname(implied["b22", 1L])
  process_cov <- if (is.null(P)) matrix(NA_real_, 2L, 2L) else P[1:2, 1:2]

  if (!is.finite(best_val) || best_val >= 1e9) {
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

  names(best_theta) <- theta_names
  names(theta_nm) <- theta_names
  if (!is.null(theta_polished)) names(theta_polished) <- theta_names

  out <- list(
    theta = best_theta,
    loss = best_val,
    implied = implied,
    implied_b11_1 = implied_b11_1,
    implied_b22_1 = implied_b22_1,
    error_b11_1 = implied_b11_1 - target_b11_1,
    error_b22_1 = implied_b22_1 - target_b22_1,
    params = list(ax = par_best$ax, ay = par_best$ay,
                  b21 = par_best$b21, b12 = par_best$b12,
                  pe_var1 = Q2[1L, 1L], pe_var2 = Q2[2L, 2L],
                  pe_cov12 = Q2[1L, 2L]),
    mats = mats,
    P = P,
    process_cov = process_cov,
    targets_b21 = targets_b21,
    targets_b12 = targets_b12,
    target_b11_1 = target_b11_1,
    target_b22_1 = target_b22_1,
    settings = list(K = K, delta_max = delta_max,
                    target_b11_1 = target_b11_1,
                    target_b22_1 = target_b22_1,
                    var1 = var1, var2 = var2, cov12 = cov12,
                    var_type = var_type,
                    stab_thresh = stab_thresh, restarts = restarts,
                    maxit = maxit, jitter_sd = jitter_sd,
                    early_tol = early_tol, seed = seed,
                    polish = polish, parallel = parallel,
                    theta0 = theta0,
                    cores_requested = cores_requested,
                    cores_used = cores_used,
                    package_version = as.character(
                      utils::packageVersion("varKgen")),
                    R_version = paste(R.version$major, R.version$minor,
                                      sep = "."),
                    RNGkind = RNGkind()),
    theta_nm = theta_nm,
    loss_nm = loss_nm,
    theta_polished = theta_polished,
    loss_polished = loss_polished,
    polish_convergence = polish_convergence,
    polish_improved = polish_improved,
    polish_attempted = polish_attempted,
    polish_status = polish_status,
    polish_message = polish_message,
    polish_fn_evals = polish_fn_evals,
    polish_seconds = polish_seconds,
    runtime_seconds = sec_total,
    restart_log = restart_log,
    restart_theta = restart_theta
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
  stopifnot(length(digits) == 1L, is.finite(digits), digits >= 0,
            digits == round(digits))
  s <- x$settings
  cat(sprintf("varK calibration (K = %d, delta_max = %d, var_type = \"%s\")\n",
              s$K, s$delta_max, s$var_type))
  cat(sprintf(paste0("AR projection effects at the first target interval ",
                     "Delta_1 = 1 time unit:\n  target  b11 = %g, ",
                     "b22 = %g\n"), s$target_b11_1, s$target_b22_1))
  cat(sprintf("  implied b11 = %.4g, b22 = %.4g\n",
              x$implied_b11_1, x$implied_b22_1))
  cat(sprintf("  error   b11 = %+.2e, b22 = %+.2e\n",
              x$error_b11_1, x$error_b22_1))
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
