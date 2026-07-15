#' Stack all overlapping delta-pairs of a simulated panel
#'
#' For a given interval \code{delta}, stacks all overlapping pairs
#' \code{(p, p + delta)} for \code{p = 1, ..., T_obs - delta} across persons
#' into a long data frame with predictors \code{Y1_p}, \code{Y2_p} and
#' outcomes \code{Y1_q}, \code{Y2_q}.
#'
#' @details
#' Because consecutive pairs share time points, the stacked observations are
#' dependent within persons. Pooled OLS point estimates remain consistent for
#' the population projection coefficients under stationarity, but naive OLS
#' standard errors are not valid; if inference is needed, cluster-robust
#' standard errors (clustering on persons) should be used.
#'
#' The original script silently produced invalid indices for
#' \code{delta >= T_obs}; this version throws an informative error.
#'
#' @param Y1,Y2 \code{N x T_obs} outcome matrices, e.g. from
#'   \code{\link{simulate_panel_stationary}}.
#' @param delta Interval (positive integer smaller than \code{ncol(Y1)}).
#'
#' @return A data frame with columns \code{Y1_p}, \code{Y2_p}, \code{Y1_q},
#'   \code{Y2_q} and \code{N * (T_obs - delta)} rows.
#'
#' @export
make_delta_pairs_overlapping <- function(Y1, Y2, delta) {
  stopifnot(is.matrix(Y1), is.matrix(Y2), all(dim(Y1) == dim(Y2)),
            length(delta) == 1L, is.finite(delta),
            delta >= 1, delta == round(delta))
  TT <- ncol(Y1)
  if (delta >= TT) {
    stop(sprintf("delta (%d) must be smaller than the number of time points (%d).",
                 as.integer(delta), TT))
  }
  p0 <- seq_len(TT - delta)
  data.frame(
    Y1_p = as.vector(Y1[, p0, drop = FALSE]),
    Y2_p = as.vector(Y2[, p0, drop = FALSE]),
    Y1_q = as.vector(Y1[, p0 + delta, drop = FALSE]),
    Y2_q = as.vector(Y2[, p0 + delta, drop = FALSE])
  )
}

#' OLS estimates of the interval-specific CLPM coefficients
#'
#' Fits the two cross-lagged regressions
#' \code{Y1_q ~ Y1_p + Y2_p} and \code{Y2_q ~ Y1_p + Y2_p} by ordinary least
#' squares and returns the four slope coefficients.
#'
#' @details
#' Both regressions share the same design matrix, so a single QR
#' decomposition with a two-column response (\code{stats::lm.fit}) is used
#' instead of two separate \code{lm()} calls with formula overhead. The
#' returned coefficients are numerically identical to those from \code{lm()}.
#'
#' @param df Data frame with columns \code{Y1_p}, \code{Y2_p}, \code{Y1_q},
#'   \code{Y2_q}, e.g. from \code{\link{make_delta_pairs_overlapping}}.
#'
#' @return Named numeric vector with elements \code{b11} (Y1 -> Y1),
#'   \code{b21} (Y2 -> Y1), \code{b22} (Y2 -> Y2), \code{b12} (Y1 -> Y2),
#'   matching the naming of the original script.
#'
#' @export
fit_clpm_ols <- function(df) {
  stopifnot(is.data.frame(df),
            all(c("Y1_p", "Y2_p", "Y1_q", "Y2_q") %in% names(df)),
            nrow(df) >= 3L)
  X <- cbind(1, df$Y1_p, df$Y2_p)
  colnames(X) <- c("(Intercept)", "Y1_p", "Y2_p")
  Ymat <- cbind(df$Y1_q, df$Y2_q)
  fit <- stats::lm.fit(X, Ymat)
  cf <- fit$coefficients
  c(b11 = unname(cf["Y1_p", 1L]),
    b21 = unname(cf["Y2_p", 1L]),
    b22 = unname(cf["Y2_p", 2L]),
    b12 = unname(cf["Y1_p", 2L]))
}
