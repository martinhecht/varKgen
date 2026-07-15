\name{default_theta0}
\alias{default_theta0}
\title{Default start values for the calibration}
\usage{
default_theta0(K, targets_b21, targets_b12)
}
\arguments{
  \item{K}{Number of lags.}
  \item{targets_b21}{Target profile for b21 (only the first element is used).}
  \item{targets_b12}{Target profile for b12 (only the first element is used).}
}
\value{Numeric vector of length \code{4K - 1}.}
\description{
Reproduces the start-value heuristic of the original script. The parameter
layout is \code{(ax[2..K], ay[2..K], b21[1..K], b12[1..K], rho_raw)}.
}
\examples{
default_theta0(K = 5,
               targets_b21 = c(0.30, 0.25, 0.15, 0.12, 0.10),
               targets_b12 = c(0.15, 0.20, 0.10, 0.10, 0.08))
}
