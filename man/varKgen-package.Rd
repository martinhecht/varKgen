\name{varKgen-package}
\alias{varKgen-package}
\alias{varKgen}
\docType{package}
\title{Calibrated VAR(K) Data-Generating Processes for Interval-Specific Cross-Lagged Effects}
\description{
Tools to construct bivariate VAR(K) data-generating processes whose
model-implied population cross-lagged panel model (CLPM) coefficients
across time intervals delta = 1, ..., delta_max match user-specified
target profiles ("varK generator").
}
\details{
Typical workflow:
\enumerate{
  \item \code{\link{calibrate_varK}}: find VAR(K) parameters whose implied
    population cross-lag profiles match the targets (lag-1 autoregressions
    of the DGP held fixed).
  \item \code{\link{simulate_panel_stationary}}: simulate stationary panels
    from the calibrated process.
  \item \code{\link{make_delta_pairs_overlapping}} and
    \code{\link{fit_clpm_ols}}: estimate interval-specific CLPM
    coefficients from the simulated data.
  \item \code{\link{verify_varK}}: convenience wrapper that compares
    estimates against targets and implied population values.
}
Lower-level building blocks: \code{\link{make_mats_varK}},
\code{\link{stationary_cov}}, \code{\link{implied_clpm_profile}},
\code{\link{coeffs_delta}}, \code{\link{default_theta0}}.
}
\author{Martin Hecht}
\keyword{package}
