test_that("theta packing/unpacking is consistent", {
  th <- varKgen:::unpack_theta(c(0.1, 0.05, 0.3, 0.1, 0.15, 0.05, atanh(0.2)),
                               K = 2, alpha11_1 = 0.5, alpha22_1 = 0.4)
  expect_equal(th$ax,  c(0.5, 0.1))
  expect_equal(th$ay,  c(0.4, 0.05))
  expect_equal(th$b21, c(0.3, 0.1))
  expect_equal(th$b12, c(0.15, 0.05))
  expect_equal(th$cov_ue, 0.2, tolerance = 1e-12)
  expect_equal(length(default_theta0(2, c(0.3, 0.2), c(0.15, 0.1))), 7L)
  expect_equal(length(default_theta0(1, 0.3, 0.15)), 3L)
})

test_that("calibration recovers an attainable target profile (smoke test)", {
  # ground truth: a stable VAR(2); its implied profile is used as target,
  # so the targets are exactly attainable
  ax  <- c(0.5, 0.10)
  ay  <- c(0.5, 0.05)
  b21 <- c(0.25, 0.05)
  b12 <- c(0.12, 0.03)
  m <- make_mats_varK(ax, ay, b21, b12, cov_ue = 0.1)
  P <- stationary_cov(m$F, m$Q)
  prof <- implied_clpm_profile(m$F, P, delta_max = 3)
  targets_b21 <- unname(prof["b21", ])
  targets_b12 <- unname(prof["b12", ])

  theta_true <- c(ax[2], ay[2], b21, b12, atanh(0.1))
  sol <- calibrate_varK(targets_b21, targets_b12, K = 2,
                        alpha11_1 = 0.5, alpha22_1 = 0.5,
                        theta0 = theta_true + 0.02,
                        restarts = 4, maxit = 3000, jitter_sd = 0.03,
                        seed = 3, polish = TRUE, verbose = FALSE)
  expect_s3_class(sol, "varK_calibration")
  expect_lt(sol$loss, 1e-4)
  expect_lt(max(abs(sol$implied["b21", ] - targets_b21)), 5e-3)
  expect_lt(max(abs(sol$implied["b12", ] - targets_b12)), 5e-3)
  expect_equal(dim(sol$implied), c(4L, 3L))
  expect_equal(nrow(sol$restart_log) <= 4, TRUE)
})

test_that("verify_varK returns a coherent table", {
  targets_b21 <- c(0.30, 0.20)
  targets_b12 <- c(0.15, 0.12)
  sol <- calibrate_varK(targets_b21, targets_b12, K = 2,
                        alpha11_1 = 0.5, alpha22_1 = 0.5,
                        restarts = 3, maxit = 2000, seed = 2,
                        verbose = FALSE)
  v <- verify_varK(sol, N = 3000, T_obs = 12, seed = 99)
  expect_equal(nrow(v$table), 2L)
  expect_true(all(c("diff_b21_vs_implied", "diff_b12_vs_target")
                  %in% names(v$table)))
  # with N = 3000 and pooled pairs, estimates should be near the implied values
  expect_lt(max(abs(v$table$diff_b21_vs_implied)), 0.1)
  expect_lt(max(abs(v$table$diff_b12_vs_implied)), 0.1)
})
