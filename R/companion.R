#' Companion-form matrices of a bivariate VAR(K)
#'
#' Builds the companion-form transition matrix \code{F} and innovation
#' covariance \code{Q} of a bivariate VAR(K) process with state vector
#' \code{s_p = (Y1_p, Y2_p, Y1_{p-1}, Y2_{p-1}, ..., Y1_{p-K+1}, Y2_{p-K+1})}.
#'
#' @details
#' The first two rows of \code{F} contain the VAR(K) coefficients: for lag
#' \code{k}, \code{ax[k]} is the autoregressive effect of \code{Y1_{p-k}} on
#' \code{Y1_p}, \code{ay[k]} of \code{Y2_{p-k}} on \code{Y2_p}, \code{b21[k]}
#' is the cross-lagged effect of \code{Y2_{p-k}} on \code{Y1_p} (process
#' 2 -> 1), and \code{b12[k]} of \code{Y1_{p-k}} on \code{Y2_p} (process
#' 1 -> 2). The remaining rows are shift registers. Only the top-left 2 x 2
#' block of \code{Q} is non-zero (innovations enter the current values only).
#'
#' Compared to the original script, the arguments \code{b} and \code{c} were
#' renamed to \code{b21} and \code{b12} (avoiding masking of \code{base::c})
#' and the construction no longer fails for \code{K = 1} (the original loop
#' \code{for (k in 2:K)} evaluated to \code{2:1 = c(2, 1)} in that case and
#' produced an out-of-bounds error).
#'
#' @param ax Numeric vector of length K: autoregressive coefficients of Y1.
#' @param ay Numeric vector of length K: autoregressive coefficients of Y2.
#' @param b21 Numeric vector of length K: cross-lagged effects Y2 -> Y1.
#' @param b12 Numeric vector of length K: cross-lagged effects Y1 -> Y2.
#' @param var_u Innovation variance of Y1 (default 1).
#' @param var_e Innovation variance of Y2 (default 1).
#' @param cov_ue Innovation covariance (default 0); must satisfy
#'   \code{abs(cov_ue) <= sqrt(var_u * var_e)}.
#'
#' @return A list with elements \code{F} (2K x 2K transition matrix) and
#'   \code{Q} (2K x 2K innovation covariance).
#'
#' @examples
#' m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.5, 0.05),
#'                     b21 = c(0.3, 0.05), b12 = c(0.15, 0.03),
#'                     cov_ue = 0.2)
#' dim(m$F)
#'
#' @export
make_mats_varK <- function(ax, ay, b21, b12, var_u = 1, var_e = 1, cov_ue = 0) {
  K <- length(ax)
  stopifnot(K >= 1L,
            length(ay) == K, length(b21) == K, length(b12) == K,
            all(is.finite(c(ax, ay, b21, b12, var_u, var_e, cov_ue))),
            var_u > 0, var_e > 0,
            abs(cov_ue) <= sqrt(var_u * var_e))
  dimS <- 2L * K

  Fmat <- matrix(0, dimS, dimS)

  # top rows: next (Y1, Y2)
  for (k in seq_len(K)) {
    i1 <- 2L * (k - 1L) + 1L  # Y1 at lag k
    i2 <- 2L * (k - 1L) + 2L  # Y2 at lag k
    Fmat[1L, i1] <- ax[k]
    Fmat[1L, i2] <- b21[k]
    Fmat[2L, i1] <- b12[k]
    Fmat[2L, i2] <- ay[k]
  }

  # shift registers (only needed for K >= 2; original code failed for K = 1)
  if (K >= 2L) {
    for (k in 2L:K) {
      Fmat[2L * (k - 1L) + 1L, 2L * (k - 2L) + 1L] <- 1
      Fmat[2L * (k - 1L) + 2L, 2L * (k - 2L) + 2L] <- 1
    }
  }

  Qmat <- matrix(0, dimS, dimS)
  Qmat[1L, 1L] <- var_u
  Qmat[2L, 2L] <- var_e
  Qmat[1L, 2L] <- Qmat[2L, 1L] <- cov_ue

  list(F = Fmat, Q = Qmat)
}
