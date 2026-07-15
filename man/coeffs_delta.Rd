\name{coeffs_delta}
\alias{coeffs_delta}
\title{Implied population CLPM coefficients for a single interval}
\usage{
coeffs_delta(Fmat, P, delta)
}
\arguments{
  \item{Fmat}{Companion-form transition matrix.}
  \item{P}{Stationary state covariance (may be \code{NULL}).}
  \item{delta}{Interval (positive integer).}
}
\value{
Named numeric vector with elements \code{b11}, \code{b21}, \code{b22},
\code{b12} (same order and names as in the original script).
}
\description{
Convenience wrapper around \code{\link{implied_clpm_profile}} for one
interval.
}
\examples{
m <- make_mats_varK(ax = 0.6, ay = 0.5, b21 = 0.15, b12 = 0.1)
P <- stationary_cov(m$F, m$Q)
coeffs_delta(m$F, P, delta = 1)  # equals the VAR(1) coefficients exactly
}
