\name{make_mats_varK}
\alias{make_mats_varK}
\title{Companion-form matrices of a bivariate VAR(K)}
\usage{
make_mats_varK(ax, ay, b21, b12, var_u = 1, var_e = 1, cov_ue = 0)
}
\arguments{
  \item{ax}{Numeric vector of length K: autoregressive coefficients of Y1.}
  \item{ay}{Numeric vector of length K: autoregressive coefficients of Y2.}
  \item{b21}{Numeric vector of length K: cross-lagged effects Y2 -> Y1.}
  \item{b12}{Numeric vector of length K: cross-lagged effects Y1 -> Y2.}
  \item{var_u}{Innovation variance of Y1 (default 1).}
  \item{var_e}{Innovation variance of Y2 (default 1).}
  \item{cov_ue}{Innovation covariance (default 0); must satisfy
    \code{abs(cov_ue) <= sqrt(var_u * var_e)}.}
}
\value{
A list with elements \code{F} (2K x 2K companion transition matrix) and
\code{Q} (2K x 2K innovation covariance, non-zero only in the top-left
2 x 2 block).
}
\description{
Builds the companion-form transition matrix \code{F} and innovation
covariance \code{Q} of a bivariate VAR(K) with state vector
\code{s_p = (Y1_p, Y2_p, Y1_p-1, Y2_p-1, ..., Y1_p-K+1, Y2_p-K+1)}.
}
\details{
Compared to the original script, the arguments \code{b} and \code{c} were
renamed to \code{b21} and \code{b12}, and the construction no longer fails
for \code{K = 1} (the original loop \code{for (k in 2:K)} produced an
out-of-bounds error in that case).
}
\examples{
m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.5, 0.05),
                    b21 = c(0.3, 0.05), b12 = c(0.15, 0.03))
dim(m$F)
}
