\name{simulate_panel_stationary}
\alias{simulate_panel_stationary}
\title{Simulate a stationary bivariate panel from a VAR(K) companion system}
\usage{
simulate_panel_stationary(N, T_obs, mats, P_stationary, seed = NULL)
}
\arguments{
  \item{N}{Number of persons.}
  \item{T_obs}{Number of time points (named \code{T} in the original
    script).}
  \item{mats}{List with companion matrices \code{F} and \code{Q}.}
  \item{P_stationary}{Stationary state covariance for the initial draw.}
  \item{seed}{Optional seed (default \code{NULL}: current RNG state).}
}
\value{A list with \code{N x T_obs} matrices \code{Y1} and \code{Y2}.}
\description{
Simulates \code{N} independent persons over \code{T_obs} time points,
starting from the stationary distribution so that no burn-in is required.
}
\details{
Multivariate normal draws use an exact symmetric positive semi-definite
factor computed once (for the initial state and for the 2 x 2 process error
covariance) instead of \code{MASS::mvrnorm} per time step. The MASS
dependency is removed. The random number stream differs from the original
script, and from package versions up to 0.2.0, for the same seed, while the
distribution of the simulated data is exactly the specified one, including
the singular positive semi-definite case.
}
\examples{
m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.4, 0.05),
                    b21 = c(0.2, 0.05), b12 = c(0.1, 0.02))
P <- stationary_cov(m$F, m$Q)
sim <- simulate_panel_stationary(N = 100, T_obs = 10, mats = m,
                                 P_stationary = P, seed = 1)
dim(sim$Y1)
}
