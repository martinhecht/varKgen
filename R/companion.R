#' Companion transition matrix of a bivariate VAR(K)
#'
#' Internal. Builds only the companion-form transition matrix \code{F}; used
#' by \code{\link{make_mats_varK}} and by the calibration loss, which combines
#' \code{F} with a process error covariance chosen according to
#' \code{var_type}.
#'
#' @param ax,ay Autoregressive coefficient vectors of length K.
#' @param b21,b12 Cross-lagged coefficient vectors of length K
#'   (b21: process 2 -> 1, b12: process 1 -> 2).
#'
#' @return The 2K x 2K companion transition matrix.
#' @keywords internal
companion_F <- function(ax, ay, b21, b12) {
  K <- length(ax)
  stopifnot(K >= 1L,
            length(ay) == K, length(b21) == K, length(b12) == K,
            all(is.finite(c(ax, ay, b21, b12))))
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

  Fmat
}

#' Companion-form matrices of a bivariate VAR(K)
#'
#' Builds the companion-form transition matrix \code{F} and process error
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
#' block of \code{Q} is non-zero (process errors enter the current values
#' only).
#'
#' The arguments \code{pe_var1}, \code{pe_var2} and \code{pe_cov12} specify
#' the process error covariance directly. If instead the stationary process
#' covariance is to be fixed, use \code{\link{calibrate_varK}} with
#' \code{var_type = "process"}, which derives the process error covariance
#' internally (see also the internal solver
#' \code{varKgen:::solve_process_error}).
#'
#' @param ax Numeric vector of length K: autoregressive coefficients of Y1.
#' @param ay Numeric vector of length K: autoregressive coefficients of Y2.
#' @param b21 Numeric vector of length K: cross-lagged effects Y2 -> Y1.
#' @param b12 Numeric vector of length K: cross-lagged effects Y1 -> Y2.
#' @param pe_var1 Process error variance of Y1 (default 1).
#' @param pe_var2 Process error variance of Y2 (default 1).
#' @param pe_cov12 Process error covariance (default 0); must satisfy
#'   \code{abs(pe_cov12) <= sqrt(pe_var1 * pe_var2)}.
#'
#' @return A list with elements \code{F} (2K x 2K transition matrix) and
#'   \code{Q} (2K x 2K process error covariance).
#'
#' @examples
#' m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.5, 0.05),
#'                     b21 = c(0.3, 0.05), b12 = c(0.15, 0.03),
#'                     pe_cov12 = 0.2)
#' dim(m$F)
#'
#' @export
make_mats_varK <- function(ax, ay, b21, b12,
                           pe_var1 = 1, pe_var2 = 1, pe_cov12 = 0) {
  stopifnot(length(pe_var1) == 1L, is.finite(pe_var1), pe_var1 > 0,
            length(pe_var2) == 1L, is.finite(pe_var2), pe_var2 > 0,
            length(pe_cov12) == 1L, is.finite(pe_cov12),
            abs(pe_cov12) <= sqrt(pe_var1 * pe_var2))

  Fmat <- companion_F(ax, ay, b21, b12)
  dimS <- nrow(Fmat)

  Qmat <- matrix(0, dimS, dimS)
  Qmat[1L, 1L] <- pe_var1
  Qmat[2L, 2L] <- pe_var2
  Qmat[1L, 2L] <- Qmat[2L, 1L] <- pe_cov12

  list(F = Fmat, Q = Qmat)
}
