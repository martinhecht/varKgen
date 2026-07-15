\name{stationary_cov}
\alias{stationary_cov}
\title{Stationary covariance of a stable linear system}
\usage{
stationary_cov(Fmat, Qmat, method = c("kron", "iterate"),
               tol = 1e-12, maxit = 400000L)
}
\arguments{
  \item{Fmat}{Square transition matrix.}
  \item{Qmat}{Innovation covariance (same dimension as \code{Fmat}).}
  \item{method}{\code{"kron"} (exact, default) or \code{"iterate"}
    (fixed-point iteration of the original script).}
  \item{tol}{Convergence tolerance for \code{method = "iterate"}.}
  \item{maxit}{Maximum iterations for \code{method = "iterate"}.}
}
\value{
The stationary covariance \code{P} solving \code{P = F P F' + Q}, or
\code{NULL} if no reliable solution was found (e.g. unstable \code{F}).
}
\description{
Solves the discrete-time Lyapunov equation for the stationary state
covariance. The default method solves the linear system
\code{(I - F \%x\% F) vec(P) = vec(Q)} exactly; its cost does not depend on
the spectral radius, whereas the fixed-point iteration can need thousands
of iterations for highly persistent systems.
}
\examples{
m <- make_mats_varK(ax = c(0.5, 0.1), ay = c(0.4, 0.05),
                    b21 = c(0.2, 0.05), b12 = c(0.1, 0.02))
P <- stationary_cov(m$F, m$Q)
}
