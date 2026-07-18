\name{calibrate_varK}
\alias{calibrate_varK}
\alias{print.varK_calibration}
\title{Calibrate a VAR(K) data-generating process to target cross-lag profiles}
\usage{
calibrate_varK(targets_b21, targets_b12, K, alpha11_1, alpha22_1,
               var1 = 1, var2 = 1, cov12 = 0,
               var_type = c("process", "process_error"),
               theta0 = NULL, stab_thresh = 0.998,
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
  \item{var1, var2}{Variances (both default 1), interpreted according to
    \code{var_type}; must be positive.}
  \item{cov12}{Covariance (default 0), interpreted according to
    \code{var_type}; must satisfy \code{cov12^2 < var1 * var2}.}
  \item{var_type}{\code{"process"} (default): \code{var1}, \code{var2},
    \code{cov12} fix the stationary process covariance (top-left 2 x 2 block
    of \code{P}); the process error covariance is derived internally.
    \code{"process_error"}: they fix the process error covariance directly.}
  \item{theta0}{Optional start vector of length \code{4K - 2}; defaults to
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
\code{b11}, \code{b21}, \code{b22}, \code{b12}), \code{params} (\code{ax},
\code{ay}, \code{b21}, \code{b12} plus the process error covariance actually
used: \code{pe_var1}, \code{pe_var2}, \code{pe_cov12}), \code{mats},
\code{P}, \code{process_cov} (achieved stationary process covariance, the
top-left 2 x 2 block of \code{P}), \code{targets_b21}, \code{targets_b12},
\code{settings}, \code{runtime_seconds}, \code{restart_log}.
}
\description{
Multi-start Nelder-Mead search (optionally BFGS-polished) for VAR(K)
parameters whose implied population cross-lagged CLPM coefficients match
the supplied target profiles across intervals 1..delta_max. The
contemporaneous 2 x 2 covariance is always fully specified via \code{var1},
\code{var2}, \code{cov12}; \code{var_type} selects whether these values fix
the stationary process covariance or the process error covariance.
}
\details{
With \code{var_type = "process"}, the required process error covariance is
derived exactly in every loss evaluation by an internal linear solver (also
for unequal variances). Setting \code{var1 = var2 = 1} fixes both process
variances at one, so standardized and unstandardized cross-lagged
coefficients coincide and target profiles can be specified in the
standardized metric. Not every target process covariance is admissible for
given dynamics; parameter values whose derived process error covariance is
not positive semi-definite receive a graded penalty. With
\code{var_type = "process_error"}, the stationary process covariance is a
result and is returned as \code{process_cov}.

The process error correlation is no longer a free calibration parameter, so
the number of free parameters is \code{4K - 2} while the number of target
values is \code{2 * delta_max}; if \code{2 * delta_max > 4K - 2}, targets
can in general only be matched approximately. If the best loss is at the
penalty level (>= 1e9), no stable and admissible solution was found and a
warning is issued.
}
\examples{
\donttest{
targets_b21 <- c(0.30, 0.25, 0.15, 0.12, 0.10)
targets_b12 <- c(0.15, 0.20, 0.10, 0.10, 0.08)
sol <- calibrate_varK(targets_b21, targets_b12, K = 3,
                      alpha11_1 = 0.5, alpha22_1 = 0.5,
                      var1 = 1, var2 = 1, cov12 = 0,
                      var_type = "process",
                      restarts = 10, maxit = 5000, verbose = FALSE)
print(sol)
sol$process_cov
}
}
