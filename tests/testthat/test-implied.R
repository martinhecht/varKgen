test_that("for K = 1 and delta = 1 the implied CLPM equals the VAR(1) matrix", {
  m <- make_mats_varK(ax = 0.6, ay = 0.5, b21 = 0.15, b12 = 0.1,
                      pe_cov12 = 0.2)
  P <- stationary_cov(m$F, m$Q)
  cd <- coeffs_delta(m$F, P, 1)
  expect_equal(unname(cd["b11"]), 0.60, tolerance = 1e-10)
  expect_equal(unname(cd["b21"]), 0.15, tolerance = 1e-10)
  expect_equal(unname(cd["b12"]), 0.10, tolerance = 1e-10)
  expect_equal(unname(cd["b22"]), 0.50, tolerance = 1e-10)
})

test_that("profile matches per-delta computation and explicit matrix powers", {
  m <- make_mats_varK(ax = c(0.5, 0.10), ay = c(0.4, 0.05),
                      b21 = c(0.2, 0.05), b12 = c(0.1, 0.02))
  P <- stationary_cov(m$F, m$Q)
  prof <- implied_clpm_profile(m$F, P, delta_max = 4)
  for (d in 1:4) {
    expect_equal(prof[, d], coeffs_delta(m$F, P, d))
  }
  # reference computation at delta = 3 via explicit matrix power
  Fd <- m$F %*% m$F %*% m$F
  C  <- (Fd %*% P)[1:2, 1:2]
  B  <- C %*% solve(P[1:2, 1:2])
  expect_equal(unname(prof["b11", 3]), B[1, 1], tolerance = 1e-10)
  expect_equal(unname(prof["b21", 3]), B[1, 2], tolerance = 1e-10)
  expect_equal(unname(prof["b22", 3]), B[2, 2], tolerance = 1e-10)
  expect_equal(unname(prof["b12", 3]), B[2, 1], tolerance = 1e-10)
})

test_that("NULL stationary covariance yields NA coefficients", {
  m <- make_mats_varK(ax = 0.5, ay = 0.5, b21 = 0.1, b12 = 0.1)
  expect_true(all(is.na(coeffs_delta(m$F, NULL, 2))))
  expect_true(all(is.na(implied_clpm_profile(m$F, NULL, 3))))
})
