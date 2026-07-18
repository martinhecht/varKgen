\name{make_mats_varK}
\alias{make_mats_varK}
\title{Companion-form matrices of a bivariate VAR(K)}
\usage{
make_mats_varK(ax, ay, b21, b12, pe_var1 = 1, pe_var2 = 1, pe_cov12 = 0)
}
\arguments{
  \item{ax}{Numeric vector of length K: autoregressive coefficients of Y1.}
  \item{ay}{Numeric vector of length K: autoregressive coefficients of Y2.}
  \item{b21}{Numeric vector of length K: cross-lagged effects Y2 -> Y1.}
  \item{b12}{Numeric vector of length K: cross-lagged effects Y1 -> Y2.}
  \item{pe_var1}{Process error variance of Y1 (default 1).}
  \item{pe_var2}{Process error variance of Y2 (default 1).}
  \item{pe_cov12}{Process error covariance (default 0); must satisfy
    \code{abs(pe_cov12) <= sqrt(pe_var1 * pe_var2)}.}
}
\value{
A list with elements \code{F} (2K x 2K companion transition matrix) and
\code{Q} (2K x 2K process error covariance, non-zero only in the top-left
2 x 2 block).
}
\description{
Builds the companion-form transition matrix \code{F} and process error
covariance \code{Q} of a bivariate VAR(K) with state vector
\code{s_p = (Y1_p, Y2_p, Y1_p-1, Y2_p-1, ..., Y1_p-K+1, Y2_p-K+1)}.
}
\details{
The arguments \code{pe_var1}, \code{pe_var2} and \code{pe_cov12} specify the
process error covariance directly. If instead the stationary process
covariance is to be fixed (e.g., unit process variances so that standardized
and unstandardized coefficients coincide), use \code{\link{calibrate_varK}}
with \code{var_type = "process"}, which derives the process error covariance
internally.

The construction does not fail for \code{K = 1} (the loop
\code{for (k in 2:K)} of the original script produced an out-of-bounds error
in that case).
}
\examples{
m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.5, 0.05),
                    b21 = c(0.3, 0.05), b12 = c(0.15, 0.03),
                    pe_cov12 = 0.2)
dim(m$F)
}
