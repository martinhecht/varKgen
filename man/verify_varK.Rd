\name{verify_varK}
\alias{verify_varK}
\title{Verify a calibrated VAR(K) generator on simulated data}
\usage{
verify_varK(calibration, N, T_obs, deltas = NULL, seed = NULL,
            return_data = FALSE)
}
\arguments{
  \item{calibration}{Object of class \code{"varK_calibration"}.}
  \item{N}{Number of persons to simulate.}
  \item{T_obs}{Number of time points to simulate.}
  \item{deltas}{Intervals to evaluate (default \code{1:delta_max}); must be
    smaller than \code{T_obs} and at most \code{delta_max}.}
  \item{seed}{Optional simulation seed.}
  \item{return_data}{If \code{TRUE}, return the simulated panel as
    element \code{sim}.}
}
\value{
List with \code{estimates} (4 x length(deltas) matrix) and \code{table}
(data frame with targets, implied values, estimates and differences per
interval).
}
\description{
Simulates one stationary panel, fits the interval-specific CLPM
regressions by pooled OLS, and tabulates estimates against both targets
and implied population values.
}
\details{
The comparison against the implied population coefficients isolates
sampling error (relevant for verifying the pipeline); the comparison
against the targets additionally contains the calibration error.
}
