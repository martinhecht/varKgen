#' Simulate a stationary bivariate panel from a VAR(K) companion system
#'
#' Simulates \code{N} independent persons over \code{T_obs} time points from
#' the companion-form system \code{s_p = F s_{p-1} + eta_p}, starting each
#' person's state from the stationary distribution \code{N(0, P_stationary)}
#' so that no burn-in is required and all time points are stationary.
#'
#' @details
#' Compared to the original script, the multivariate normal draws use a
#' Cholesky factor computed once (for the initial state and once for the
#' 2 x 2 innovation covariance) instead of calling \code{MASS::mvrnorm} in
#' every time step (which performs an eigendecomposition per call). This
#' removes the MASS dependency and speeds up large simulations. Note that the
#' random number stream therefore differs from the original script for the
#' same seed; distributional properties are identical.
#'
#' @param N Number of persons.
#' @param T_obs Number of time points (named \code{T} in the original script;
#'   renamed to avoid masking \code{TRUE}'s shorthand).
#' @param mats List with companion matrices \code{F} and \code{Q}, e.g. from
#'   \code{\link{make_mats_varK}} or the \code{mats} element of a
#'   \code{\link{calibrate_varK}} result.
#' @param P_stationary Stationary state covariance, e.g. from
#'   \code{\link{stationary_cov}} or the \code{P} element of a calibration
#'   result.
#' @param seed Optional seed set before drawing (default \code{NULL}: the
#'   current RNG state is used).
#'
#' @return A list with \code{N x T_obs} matrices \code{Y1} and \code{Y2}.
#'
#' @examples
#' m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.4, 0.05),
#'                     b21 = c(0.2, 0.05), b12 = c(0.1, 0.02))
#' P <- stationary_cov(m$F, m$Q)
#' sim <- simulate_panel_stationary(N = 100, T_obs = 10, mats = m,
#'                                  P_stationary = P, seed = 1)
#' dim(sim$Y1)
#'
#' @export
simulate_panel_stationary <- function(N, T_obs, mats, P_stationary, seed = NULL) {
  stopifnot(length(N) == 1L, N >= 1L, N == round(N),
            length(T_obs) == 1L, T_obs >= 1L, T_obs == round(T_obs))
  if (is.null(P_stationary)) {
    stop("Need stationary covariance 'P_stationary' for the initial draw.")
  }
  Fmat <- mats$F
  Qmat <- mats$Q
  stopifnot(is.matrix(Fmat), nrow(Fmat) == ncol(Fmat),
            all(dim(Qmat) == dim(Fmat)),
            all(dim(P_stationary) == dim(Fmat)))
  dimS <- nrow(Fmat)

  if (!is.null(seed)) set.seed(seed)

  L_P <- chol_psd(P_stationary)                       # t(L_P) %*% L_P = P
  L_Q <- chol_psd(Qmat[1:2, 1:2, drop = FALSE])       # innovation covariance

  # initial states from the stationary distribution
  S <- matrix(rnorm(N * dimS), N, dimS) %*% L_P

  Y1 <- matrix(NA_real_, N, T_obs)
  Y2 <- matrix(NA_real_, N, T_obs)
  Y1[, 1L] <- S[, 1L]
  Y2[, 1L] <- S[, 2L]

  if (T_obs >= 2L) {
    tF <- t(Fmat)
    for (p in 2L:T_obs) {
      S <- S %*% tF
      S[, 1:2] <- S[, 1:2] + matrix(rnorm(2L * N), N, 2L) %*% L_Q
      Y1[, p] <- S[, 1L]
      Y2[, p] <- S[, 2L]
    }
  }

  list(Y1 = Y1, Y2 = Y2)
}
