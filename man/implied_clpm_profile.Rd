\name{implied_clpm_profile}
\alias{implied_clpm_profile}
\title{Implied population CLPM coefficients for all intervals 1..delta_max}
\usage{
implied_clpm_profile(Fmat, P, delta_max)
}
\arguments{
  \item{Fmat}{Companion-form transition matrix (2K x 2K).}
  \item{P}{Stationary state covariance (may be \code{NULL}, yielding NA).}
  \item{delta_max}{Largest interval.}
}
\value{
A 4 x \code{delta_max} matrix with rows \code{b11}, \code{b21},
\code{b22}, \code{b12} and one column per interval.
}
\description{
For each interval delta, computes the population coefficients of the
interval-specific cross-lagged regressions
\code{Y1_p+delta ~ Y1_p + Y2_p} and \code{Y2_p+delta ~ Y1_p + Y2_p} under
the stationary VAR(K) process, i.e. the least-squares projection
\code{B(delta) = (F^delta P)[1:2, 1:2] solve(P[1:2, 1:2])}.
}
\details{
Only the first two rows of \code{F^delta} are propagated cumulatively
across intervals, avoiding full matrix powers per interval.
}
\examples{
m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.4, 0.05),
                    b21 = c(0.2, 0.05), b12 = c(0.1, 0.02))
P <- stationary_cov(m$F, m$Q)
implied_clpm_profile(m$F, P, delta_max = 5)
}
