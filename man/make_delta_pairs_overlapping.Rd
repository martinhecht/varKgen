\name{make_delta_pairs_overlapping}
\alias{make_delta_pairs_overlapping}
\title{Stack all overlapping delta-pairs of a simulated panel}
\usage{
make_delta_pairs_overlapping(Y1, Y2, delta)
}
\arguments{
  \item{Y1, Y2}{\code{N x T_obs} outcome matrices.}
  \item{delta}{Interval (positive integer smaller than \code{ncol(Y1)});
    the original script silently produced invalid indices for
    \code{delta >= T_obs}, this version throws an informative error.}
}
\value{
Data frame with columns \code{Y1_p}, \code{Y2_p}, \code{Y1_q},
\code{Y2_q} and \code{N * (T_obs - delta)} rows.
}
\description{
Stacks all overlapping pairs \code{(p, p + delta)} across persons. Pooled
OLS point estimates on these pairs are consistent for the population
projection coefficients under stationarity, but naive OLS standard errors
are invalid (within-person dependence); use cluster-robust standard errors
if inference is required.
}
