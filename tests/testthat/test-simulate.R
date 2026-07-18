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
