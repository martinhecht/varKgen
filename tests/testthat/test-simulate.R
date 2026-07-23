test_that("simulation is reproducible with a seed", {
  m <- make_mats_varK(ax = c(0.5, 0.10), ay = c(0.4, 0.05),
                      b21 = c(0.2, 0.05), b12 = c(0.1, 0.02),
                      pe_cov12 = 0.2)
  P <- stationary_cov(m$F, m$Q)
  s1 <- simulate_panel_stationary(500, 4, m, P, seed = 7)
  s2 <- simulate_panel_stationary(500, 4, m, P, seed = 7)
  expect_identical(s1, s2)
  expect_equal(dim(s1$Y1), c(500L, 4L))
})

test_that("simulated data roughly match the stationary covariance", {
  m <- make_mats_varK(ax = c(0.5, 0.10), ay = c(0.4, 0.05),
                      b21 = c(0.2, 0.05), b12 = c(0.1, 0.02),
                      pe_cov12 = 0.2)
  P <- stationary_cov(m$F, m$Q)
  s <- simulate_panel_stationary(20000, 3, m, P, seed = 1)
  # covariance of the last time point should be close to the stationary one
  emp <- stats::cov(cbind(s$Y1[, 3], s$Y2[, 3]))
  expect_lt(max(abs(emp - P[1:2, 1:2])), 0.15 * max(abs(P[1:2, 1:2])))
})

test_that("missing stationary covariance raises an error", {
  m <- make_mats_varK(ax = 0.5, ay = 0.5, b21 = 0.1, b12 = 0.1)
  expect_error(simulate_panel_stationary(10, 5, m, NULL))
})

test_that("cov_factor_psd factors covariance matrices exactly", {
  f <- varKgen:::cov_factor_psd
  # regular
  M <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  expect_equal(crossprod(f(M)), M, tolerance = 1e-12)
  # exactly singular positive semi-definite: chol() would fail here and the
  # old jittered fallback would have simulated from M + epsilon * I
  S <- matrix(c(1, 1, 1, 1), 2, 2)
  expect_error(chol(S))
  expect_equal(crossprod(f(S)), S, tolerance = 1e-12)
  # rank-deficient 3 x 3
  R3 <- outer(c(1, 2, 3), c(1, 2, 3))
  expect_equal(crossprod(f(R3)), R3, tolerance = 1e-10)
  # no correction is applied to genuinely PSD matrices
  expect_equal(attr(f(M), "psd_correction"), 0)
  expect_equal(attr(f(S), "psd_correction"), 0)

  # a negative eigenvalue at noise level is truncated; the factor then
  # belongs to the PSD-adjusted matrix, and the correction is disclosed
  N1 <- S; N1[2, 2] <- 1 - 1e-14
  L1 <- f(N1)
  expect_gt(attr(L1, "psd_correction"), 0)
  expect_lt(attr(L1, "psd_correction"), 1e-10 * max(abs(eigen(N1)$values)))
  ee <- eigen((N1 + t(N1)) / 2, symmetric = TRUE)
  M_used <- ee$vectors %*% diag(pmax(ee$values, 0)) %*% t(ee$vectors)
  expect_equal(crossprod(L1), M_used, tolerance = 1e-12)

  # clearly indefinite is rejected
  expect_error(f(matrix(c(1, 0, 0, -0.5), 2, 2)),
               "not positive semi-definite")
  # the tolerance is purely relative: a matrix whose scale is far below 1
  # must not pass just because the negative eigenvalue is small in absolute
  # terms
  expect_error(f(diag(c(1e-12, -1e-11))), "not positive semi-definite")
  # the zero matrix is handled explicitly
  expect_equal(crossprod(f(matrix(0, 2, 2))), matrix(0, 2, 2))
})

test_that("simulate_panel_stationary validates its matrices", {
  m <- make_mats_varK(c(0.5, 0.1), c(0.4, 0.05), c(0.2, 0.05), c(0.1, 0.03))
  P <- stationary_cov(m$F, m$Q)
  # non-zero entries outside the top-left 2 x 2 block would be ignored
  bad <- m; bad$Q[3, 3] <- 0.5
  expect_error(simulate_panel_stationary(10, 5, bad, P, seed = 1),
               "top-left 2 x 2 block")
  # asymmetric covariance
  bad2 <- m; bad2$Q[1, 2] <- 0.3; bad2$Q[2, 1] <- -0.3
  expect_error(simulate_panel_stationary(10, 5, bad2, P, seed = 1),
               "not symmetric")
  # correction is disclosed and zero for clean input
  sim <- simulate_panel_stationary(20, 4, m, P, seed = 1)
  expect_equal(unname(sim$psd_correction), c(0, 0))
})

test_that("simulation works with a singular process error covariance", {
  # abs(pe_cov12) = sqrt(pe_var1 * pe_var2) is admitted by make_mats_varK
  m <- make_mats_varK(ax = 0.5, ay = 0.4, b21 = 0.2, b12 = 0.1,
                      pe_var1 = 1, pe_var2 = 1, pe_cov12 = 1)
  expect_equal(min(eigen(m$Q[1:2, 1:2], symmetric = TRUE)$values), 0,
               tolerance = 1e-12)
  P <- stationary_cov(m$F, m$Q)
  sim <- simulate_panel_stationary(N = 4000, T_obs = 6, mats = m,
                                   P_stationary = P, seed = 21)
  expect_equal(dim(sim$Y1), c(4000L, 6L))
  expect_true(all(is.finite(sim$Y1)))
  expect_true(all(is.finite(sim$Y2)))
  # the two processes share one innovation direction, so the empirical
  # process covariance must reproduce the (singular) target
  emp <- stats::cov(cbind(as.vector(sim$Y1), as.vector(sim$Y2)))
  expect_equal(emp[1, 1], P[1, 1], tolerance = 0.15)
  expect_equal(emp[1, 2], P[1, 2], tolerance = 0.15)
})
