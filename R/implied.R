#' Implied population CLPM coefficients for all intervals 1..delta_max
#'
#' Computes, for each interval \code{delta = 1, ..., delta_max}, the
#' population coefficients of the (typically misspecified) interval-specific
#' cross-lagged panel regressions
#' \code{Y1_{p+delta} ~ Y1_p + Y2_p} and \code{Y2_{p+delta} ~ Y1_p + Y2_p}
#' under the stationary VAR(K) process defined by \code{Fmat} (companion form)
#' with stationary state covariance \code{P}.
#'
#' @details
#' With \code{C(delta) = (F^delta P)[1:2, 1:2]} (the cross-covariance of the
#' current values at distance \code{delta}) and \code{V = P[1:2, 1:2]}, the
#' 2 x 2 coefficient matrix of the population least-squares projection is
#' \code{B(delta) = C(delta) V^{-1}}, so that \code{b11 = B[1,1]},
#' \code{b21 = B[1,2]} (effect of Y2 on Y1), \code{b22 = B[2,2]} and
#' \code{b12 = B[2,1]} (effect of Y1 on Y2).
#'
#' Instead of forming full matrix powers \code{F^delta} for every
#' \code{delta} (O(delta_max^2) full matrix products in the original script),
#' only the first two rows of \code{F^delta} are propagated cumulatively,
#' which reduces the cost to \code{delta_max} products of a 2 x 2K block with
#' \code{F}.
#'
#' @param Fmat Companion-form transition matrix (2K x 2K).
#' @param P Stationary state covariance, e.g. from
#'   \code{\link{stationary_cov}}; may be \code{NULL}, in which case a matrix
#'   of \code{NA} is returned.
#' @param delta_max Largest interval.
#'
#' @return A 4 x \code{delta_max} matrix with rows \code{b11}, \code{b21},
#'   \code{b22}, \code{b12} (same order as in the original script) and one
#'   column per interval.
#'
#' @examples
#' m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.4, 0.05),
#'                     b21 = c(0.2, 0.05), b12 = c(0.1, 0.02))
#' P <- stationary_cov(m$F, m$Q)
#' implied_clpm_profile(m$F, P, delta_max = 5)
#'
#' @export
implied_clpm_profile <- function(Fmat, P, delta_max) {
  stopifnot(is.matrix(Fmat), nrow(Fmat) == ncol(Fmat), nrow(Fmat) >= 2L,
            length(delta_max) == 1L, is.finite(delta_max),
            delta_max >= 1, delta_max == round(delta_max))
  out <- matrix(NA_real_, nrow = 4L, ncol = delta_max,
                dimnames = list(c("b11", "b21", "b22", "b12"),
                                paste0("delta=", seq_len(delta_max))))
  if (is.null(P)) return(out)
  stopifnot(is.matrix(P), all(dim(P) == dim(Fmat)))

  VarNow <- P[1:2, 1:2]
  Pcols  <- P[, 1:2, drop = FALSE]
  Rrow   <- Fmat[1:2, , drop = FALSE]   # (F^1)[1:2, ]

  for (d in seq_len(delta_max)) {
    CovBlock <- Rrow %*% Pcols          # 2 x 2 block of F^d P
    B <- tryCatch(t(solve(VarNow, t(CovBlock))),  # = CovBlock %*% solve(VarNow)
                  error = function(e) NULL)
    if (is.null(B) || !all(is.finite(B))) return(out)
    out["b11", d] <- B[1L, 1L]
    out["b21", d] <- B[1L, 2L]
    out["b22", d] <- B[2L, 2L]
    out["b12", d] <- B[2L, 1L]
    Rrow <- Rrow %*% Fmat               # advance to (F^{d+1})[1:2, ]
  }
  out
}

#' Implied population CLPM coefficients for a single interval
#'
#' Convenience wrapper around \code{\link{implied_clpm_profile}} that returns
#' the population coefficients for one interval \code{delta}. The name and
#' the return format (named vector \code{b11}, \code{b21}, \code{b22},
#' \code{b12}) match the function of the same name in the original script.
#'
#' @param Fmat Companion-form transition matrix.
#' @param P Stationary state covariance (may be \code{NULL}).
#' @param delta Interval (positive integer).
#'
#' @return Named numeric vector with elements \code{b11}, \code{b21},
#'   \code{b22}, \code{b12}.
#'
#' @examples
#' m <- make_mats_varK(ax = 0.6, ay = 0.5, b21 = 0.15, b12 = 0.1)
#' P <- stationary_cov(m$F, m$Q)
#' coeffs_delta(m$F, P, delta = 1)  # equals the VAR(1) coefficients exactly
#'
#' @export
coeffs_delta <- function(Fmat, P, delta) {
  stopifnot(length(delta) == 1L, is.finite(delta),
            delta >= 1, delta == round(delta))
  if (is.null(P)) {
    return(c(b11 = NA_real_, b21 = NA_real_, b22 = NA_real_, b12 = NA_real_))
  }
  prof <- implied_clpm_profile(Fmat, P, delta_max = delta)
  prof[, delta]
}
