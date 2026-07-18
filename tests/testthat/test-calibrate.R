test_that("theta packing/unpacking is consistent (4K - 2 layout)", {
  th <- varKgen:::unpack_theta(c(0.1, 0.05, 0.3, 0.1, 0.15, 0.05),
                               K = 2, alpha11_1 = 0.5, alpha22_1 = 0.4)
  expect_equal(th$ax,  c(0.5, 0.1))
  expect_equal(th$ay,  c(0.4, 0.05))
  expect_equal(th$b21, c(0.3, 0.1))
  expect_equal(th$b12, c(0.15, 0.05))
  expect_equal(length(default_theta0(2, c(0.3, 0.2), c(0.15, 0.1))), 6L)
  expect_equal(length(default_theta0(1, 0.3, 0.15)), 2L)
})

test_that("covariance input validation works", {
  expect_error(calibrate_varK(0.2, 0.1, K = 1,
                              alpha11_1 = 0.4, alpha22_1 = 0.4,
                              var1 = 1, var2 = 1, cov12 = 1))
  expect_error(calibrate_varK(0.2, 0.1, K = 1,
                              alpha11_1 = 0.4, alpha22_1 = 0.4,
                              var1 = -1))
})

test_that("calibration recovers an attainable profile in process_error mode", {
  # ground truth: a stable VAR(2) with a known process error covariance;
  # its implied profile is used as target, so the targets are exactly
  # attainable with the same covariance settings
  ax  <- c(0.5, 0.10)
  ay  <- c(0.5, 0.05)
  b21 <- c(0.25, 0.05)
  b12 <- c(0.12, 0.03)
  m <- make_mats_varK(ax, ay, b21, b12, pe_cov12 = 0.1)
  P <- stationary_cov(m$F, m$Q)
  prof <- implied_clpm_profile(m$F, P, delta_max = 3)
  targets_b21 <- unname(prof["b21", ])
  targets_b12 <- unname(prof["b12", ])

  theta_true <- c(ax[2], ay[2], b21, b12)
  sol <- calibrate_varK(targets_b21, targets_b12, K = 2,
                        alpha11_1 = 0.5, alpha22_1 = 0.5,
                        var1 = 1, var2 = 1, cov12 = 0.1,
                        var_type = "process_error",
                        theta0 = theta_true + 0.02,
                        restarts = 4, maxit = 3000, jitter_sd = 0.03,
                        seed = 3, polish = TRUE, verbose = FALSE)
  expect_s3_class(sol, "varK_calibration")
  expect_lt(sol$loss, 1e-4)
  expect_lt(max(abs(sol$implied["b21", ] - targets_b21)), 5e-3)
  expect_lt(max(abs(sol$implied["b12", ] - targets_b12)), 5e-3)
  expect_equal(sol$params$pe_cov12, 0.1)
  expect_equal(sol$mats$Q[1:2, 1:2], matrix(c(1, 0.1, 0.1, 1), 2, 2))
  expect_equal(dim(sol$implied), c(4L, 3L))
})

test_that("process mode hits the specified stationary process covariance", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        alpha11_1 = 0.5, alpha22_1 = 0.5,
                        var1 = 1, var2 = 1, cov12 = 0,
                        var_type = "process",
                        restarts = 3, maxit = 2000, seed = 2,
                        verbose = FALSE)
  expect_lt(sol$loss, 1e9)
  expect_equal(sol$process_cov[1, 1], 1, tolerance = 1e-5)
  expect_equal(sol$process_cov[2, 2], 1, tolerance = 1e-5)
  expect_lt(abs(sol$process_cov[1, 2]), 1e-5)
  # derived process error covariance is valid and consistent with mats$Q
  Q2 <- sol$mats$Q[1:2, 1:2]
  expect_equal(Q2[1, 1], sol$params$pe_var1)
  expect_equal(Q2[2, 2], sol$params$pe_var2)
  expect_equal(Q2[1, 2], sol$params$pe_cov12)
  expect_gte(min(eigen(Q2, only.values = TRUE)$values), -1e-8)
  # P is consistent with F and the derived Q
  expect_lt(max(abs(sol$P - sol$mats$F %*% sol$P %*% t(sol$mats$F) -
                      sol$mats$Q)), 1e-6)
})

test_that("process mode with unequal targets works end to end (K = 1)", {
  # for K = 1 the implied coefficients at delta = 1 equal the VAR(1)
  # coefficients exactly, so the targets are exactly attainable; the derived
  # process error covariance at the solution was verified by hand:
  # Q = P - F P F' = [[1.164, -0.214], [-0.214, 0.446]]
  sol <- calibrate_varK(0.2, 0.1, K = 1,
                        alpha11_1 = 0.4, alpha22_1 = 0.6,
                        var1 = 1.4, var2 = 0.7, cov12 = -0.1,
                        var_type = "process",
                        restarts = 2, maxit = 2000, seed = 1,
                        verbose = FALSE)
  expect_lt(sol$loss, 1e-6)
  expect_equal(sol$process_cov[1, 1], 1.4, tolerance = 1e-5)
  expect_equal(sol$process_cov[2, 2], 0.7, tolerance = 1e-5)
  expect_equal(sol$process_cov[1, 2], -0.1, tolerance = 1e-4)
  expect_equal(unname(sol$implied["b21", 1]), 0.2, tolerance = 1e-3)
  expect_equal(unname(sol$implied["b12", 1]), 0.1, tolerance = 1e-3)
  expect_equal(sol$params$pe_var1,  1.164, tolerance = 1e-2)
  expect_equal(sol$params$pe_var2,  0.446, tolerance = 1e-2)
  expect_equal(sol$params$pe_cov12, -0.214, tolerance = 1e-2)
})

test_that("verify_varK works with a process mode calibration", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        alpha11_1 = 0.5, alpha22_1 = 0.5,
                        restarts = 3, maxit = 2000, seed = 2,
                        verbose = FALSE)
  v <- verify_varK(sol, N = 3000, T_obs = 12, seed = 99)
  expect_equal(nrow(v$table), 2L)
  expect_true(all(c("diff_b21_vs_implied", "diff_b12_vs_target")
                  %in% names(v$table)))
  expect_lt(max(abs(v$table$diff_b21_vs_implied)), 0.1)
  expect_lt(max(abs(v$table$diff_b12_vs_implied)), 0.1)
})
