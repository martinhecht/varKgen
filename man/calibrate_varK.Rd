\name{calibrate_varK}
\alias{calibrate_varK}
\alias{print.varK_calibration}
\title{Calibrate a VAR(K) data-generating process to target cross-lag profiles}
\usage{
calibrate_varK(targets_b21, targets_b12, K, alpha11_1, alpha22_1,
               var_u = 1, var_e = 1, theta0 = NULL, stab_thresh = 0.998,
               restarts = 120L, maxit = 25000L, jitter_sd = 0.10,
               early_tol = 1e-10, seed = 1L, polish = TRUE,
               parallel = FALSE, cores = NULL, verbose = TRUE)

\method{print}{varK_calibration}(x, digits = 4, ...)
}
\arguments{
  \item{targets_b21}{Target profile for the cross effect 2 -> 1, length
    \code{delta_max}.}
  \item{targets_b12}{Target profile for the cross effect 1 -> 2, same
    length.}
  \item{K}{Number of lags of the VAR(K) data-generating process.}
  \item{alpha11_1, alpha22_1}{Fixed lag-1 autoregressive coefficients of
    the DGP (not the implied CLPM autoregressions).}
  \item{var_u, var_e}{Innovation variances (default 1).}
  \item{theta0}{Optional start vector of length \code{4K - 1}; defaults to
    \code{\link{default_theta0}}.}
  \item{stab_thresh}{Stability threshold for the spectral radius.}
  \item{restarts}{Number of jittered Nelder-Mead restarts.}
  \item{maxit}{Maximum Nelder-Mead iterations per restart.}
  \item{jitter_sd}{SD of the normal jitter around \code{theta0}.}
  \item{early_tol}{Early-stopping tolerance on the loss (set to 0 to
    disable early stopping, e.g. for optimizer studies).}
  \item{seed}{Random seed for the restart jitter.}
  \item{polish}{If \code{TRUE}, a BFGS run from the best solution is kept
    if it improves the loss.}
  \item{parallel}{If \code{TRUE}, restarts run on a PSOCK cluster
    (requires the package to be installed; disables early stopping).}
  \item{cores}{Number of workers; defaults to
    \code{parallel::detectCores() - 1}.}
  \item{verbose}{Print progress information.}
  \item{x}{A \code{"varK_calibration"} object.}
  \item{digits}{Digits for the printed implied coefficient table.}
  \item{...}{Ignored.}
}
\value{
An object of class \code{"varK_calibration"}: a list with elements
\code{theta}, \code{loss}, \code{implied} (4 x delta_max matrix, rows
\code{b11}, \code{b21}, \code{b22}, \code{b12}), \code{params},
\code{mats}, \code{P}, \code{targets_b21}, \code{targets_b12},
\code{settings}, \code{runtime_seconds}, \code{restart_log}.
}
\description{
Multi-start Nelder-Mead search (optionally BFGS-polished) for VAR(K)
parameters whose implied population cross-lagged CLPM coefficients match
the supplied target profiles across intervals 1..delta_max.
}
\details{
The loss is the sum of squared deviations of the implied \code{b21(delta)}
and \code{b12(delta)} from the targets. With \code{2 * delta_max} targets
and \code{4K - 1} free parameters, targets can in general only be matched
approximately if \code{2 * delta_max > 4K - 1}. Loss values at stable
points are identical to the original script; results are numerically
equivalent but not bit-identical (exact instead of iterative Lyapunov
solver, graded instead of flat instability penalty).
}
\examples{
\donttest{
targets_b21 <- c(0.30, 0.25, 0.15, 0.12, 0.10)
targets_b12 <- c(0.15, 0.20, 0.10, 0.10, 0.08)
sol <- calibrate_varK(targets_b21, targets_b12, K = 3,
                      alpha11_1 = 0.5, alpha22_1 = 0.5,
                      restarts = 10, maxit = 5000, verbose = FALSE)
print(sol)
}
}
