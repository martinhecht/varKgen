\name{fit_clpm_ols}
\alias{fit_clpm_ols}
\title{OLS estimates of the interval-specific CLPM coefficients}
\usage{
fit_clpm_ols(df)
}
\arguments{
  \item{df}{Data frame with columns \code{Y1_p}, \code{Y2_p},
    \code{Y1_q}, \code{Y2_q}, e.g. from
    \code{\link{make_delta_pairs_overlapping}}.}
}
\value{
Named numeric vector \code{b11} (Y1 -> Y1), \code{b21} (Y2 -> Y1),
\code{b22} (Y2 -> Y2), \code{b12} (Y1 -> Y2).
}
\description{
Fits \code{Y1_q ~ Y1_p + Y2_p} and \code{Y2_q ~ Y1_p + Y2_p} via a single
QR decomposition with a two-column response (\code{stats::lm.fit});
numerically identical to two separate \code{lm()} calls, but faster.
}
