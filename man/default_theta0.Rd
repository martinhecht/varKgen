\name{default_theta0}
\alias{default_theta0}
\title{Default start values for the calibration}
\usage{
default_theta0(K, targets_b21, targets_b12, target_b11_1, target_b22_1)
}
\arguments{
  \item{K}{Number of lags.}
  \item{targets_b21}{Target profile for b21 (only the first element is used).}
  \item{targets_b12}{Target profile for b12 (only the first element is used).}
  \item{target_b11_1, target_b22_1}{Targets for \code{b11(Delta_1)} and
    \code{b22(Delta_1)} at the first target interval \code{Delta_1 = 1}
    time unit; used as start values for the lag-1 autoregressive
    coefficients.}
}
\value{Numeric vector of length \code{4K}.}
\description{
Adapts the start-value heuristic of the original script to the current
parameter layout \code{(ax[1..K], ay[1..K], b21[1..K], b12[1..K])}. There is
no start value for a process error correlation: that parameter no longer
exists, since the full contemporaneous covariance is set by the user in
\code{\link{calibrate_varK}} (argument \code{var_type}).
}
\examples{
default_theta0(K = 5,
               targets_b21 = c(0.30, 0.25, 0.15, 0.12, 0.10),
               targets_b12 = c(0.15, 0.20, 0.10, 0.10, 0.08),
               target_b11_1 = 0.5, target_b22_1 = 0.5)
}
