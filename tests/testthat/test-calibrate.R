test_that("theta packing/unpacking is consistent (4K layout)", {
  th <- varKgen:::unpack_theta(c(0.5, 0.1, 0.4, 0.05,
                                 0.3, 0.1, 0.15, 0.05), K = 2)
  expect_equal(th$ax,  c(0.5, 0.1))
  expect_equal(th$ay,  c(0.4, 0.05))
  expect_equal(th$b21, c(0.3, 0.1))
  expect_equal(th$b12, c(0.15, 0.05))
  expect_equal(length(default_theta0(2, c(0.3, 0.2), c(0.15, 0.1),
                                     target_b11_1 = 0.5,
                                     target_b22_1 = 0.5)), 8L)
  expect_equal(length(default_theta0(1, 0.3, 0.15,
                                     target_b11_1 = 0.5,
                                     target_b22_1 = 0.5)), 4L)
})

test_that("removed arguments alpha11_1 / alpha22_1 give an informative error", {
  expect_error(
    calibrate_varK(0.2, 0.1, K = 1, target_b11_1 = 0.4, target_b22_1 = 0.6,
                   alpha11_1 = 0.4, verbose = FALSE),
    "no longer supported")
  expect_error(
    calibrate_varK(0.2, 0.1, K = 1, target_b11_1 = 0.4, target_b22_1 = 0.6,
                   alpha22_1 = 0.6, verbose = FALSE),
    "target_b11_1")
})

test_that("covariance input validation works", {
  expect_error(calibrate_varK(0.2, 0.1, K = 1,
                              target_b11_1 = 0.4, target_b22_1 = 0.4,
                              var1 = 1, var2 = 1, cov12 = 1))
  expect_error(calibrate_varK(0.2, 0.1, K = 1,
                              target_b11_1 = 0.4, target_b22_1 = 0.4,
                              var1 = -1))
})

test_that("b11(Delta_1) is a projection effect, not the first VAR matrix", {
  # for K = 1 the projection matrix at Delta_1 equals the VAR(1) matrix
  m1 <- make_mats_varK(0.4, 0.6, 0.2, 0.1)
  P1 <- stationary_cov(m1$F, m1$Q)
  pr1 <- implied_clpm_profile(m1$F, P1, delta_max = 1)
  expect_equal(unname(pr1["b11", 1]), 0.4, tolerance = 1e-10)
  expect_equal(unname(pr1["b22", 1]), 0.6, tolerance = 1e-10)

  # for K = 2 with non-zero higher lags they differ
  m2 <- make_mats_varK(c(0.5, 0.20), c(0.5, 0.15),
                       c(0.25, 0.05), c(0.12, 0.03))
  P2 <- stationary_cov(m2$F, m2$Q)
  pr2 <- implied_clpm_profile(m2$F, P2, delta_max = 3)
  expect_false(isTRUE(all.equal(unname(pr2["b11", 1]), 0.5,
                                tolerance = 1e-4)))
  expect_false(isTRUE(all.equal(unname(pr2["b22", 1]), 0.5,
                                tolerance = 1e-4)))
})

test_that("target intervals are the fixed grid Delta_j = j", {
  # a known VAR(1): the projection matrix at Delta_j is F^j, so column j of
  # the implied profile must equal F^j; this pins Delta_j = j and makes the
  # first column Delta_1 = 1 time unit
  m <- make_mats_varK(0.4, 0.6, 0.2, 0.1)
  P <- stationary_cov(m$F, m$Q)
  prof <- implied_clpm_profile(m$F, P, delta_max = 3)
  Fp <- diag(2)
  for (j in 1:3) {
    Fp <- Fp %*% m$F[1:2, 1:2]
    expect_equal(unname(prof["b11", j]), Fp[1, 1], tolerance = 1e-8)
    expect_equal(unname(prof["b21", j]), Fp[1, 2], tolerance = 1e-8)
    expect_equal(unname(prof["b22", j]), Fp[2, 2], tolerance = 1e-8)
    expect_equal(unname(prof["b12", j]), Fp[2, 1], tolerance = 1e-8)
  }

  # The first element of a target vector belongs to Delta_1, and the
  # autoregressive targets belong to Delta_1 as well. The targets are taken
  # from a known VAR(1) so that they are exactly attainable; note that
  # b21(Delta_1) = 0.20 differs from b21(Delta_2) = 0.16, so a shifted index
  # mapping could not satisfy both.
  Fk <- matrix(c(0.5, 0.1, 0.2, 0.3), 2, 2)   # rows: Y1, Y2
  Fk2 <- Fk %*% Fk
  t21 <- c(Fk[1, 2], Fk2[1, 2])               # 0.20, 0.16
  t12 <- c(Fk[2, 1], Fk2[2, 1])               # 0.10, 0.08
  sol <- calibrate_varK(t21, t12, K = 1,
                        target_b11_1 = Fk[1, 1], target_b22_1 = Fk[2, 2],
                        restarts = 3, maxit = 2000, seed = 4,
                        verbose = FALSE)
  # tolerances are loose enough for optimizer differences across platforms
  # but far tighter than the 0.04 gap between b21(Delta_1) and b21(Delta_2)
  expect_lt(sol$loss, 1e-6)
  expect_equal(unname(sol$implied["b21", ]), t21, tolerance = 1e-3)
  expect_equal(unname(sol$implied["b12", ]), t12, tolerance = 1e-3)
  expect_equal(sol$implied_b11_1, Fk[1, 1], tolerance = 1e-3)
  expect_equal(sol$implied_b22_1, Fk[2, 2], tolerance = 1e-3)
  # the autoregressive targets are read at Delta_1, i.e. the first column
  expect_equal(sol$implied_b11_1, unname(sol$implied["b11", 1]))
  expect_equal(sol$implied_b22_1, unname(sol$implied["b22", 1]))
  # for K = 1 the free VAR coefficients coincide with B(Delta_1)
  expect_equal(sol$params$ax[1], Fk[1, 1], tolerance = 1e-3)
  expect_equal(sol$params$ay[1], Fk[2, 2], tolerance = 1e-3)

  # no user-supplied interval vector is supported
  expect_false(any(c("deltas", "delta", "Delta") %in%
                     names(formals(calibrate_varK))))
})

test_that("the calibration object reports target, implied and error", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        target_b11_1 = 0.5, target_b22_1 = 0.4,
                        restarts = 3, maxit = 2000, seed = 2,
                        verbose = FALSE)
  expect_equal(sol$target_b11_1, 0.5)
  expect_equal(sol$target_b22_1, 0.4)
  expect_equal(sol$error_b11_1, sol$implied_b11_1 - sol$target_b11_1)
  expect_equal(sol$error_b22_1, sol$implied_b22_1 - sol$target_b22_1)
})

test_that("calibration recovers an attainable target in process_error mode", {
  # ground truth: a stable VAR(2) with a known process error covariance;
  # its implied coefficients are used as targets, so the targets are exactly
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
  target_b11_1 <- unname(prof["b11", 1])
  target_b22_1 <- unname(prof["b22", 1])

  theta_true <- c(ax, ay, b21, b12)
  sol <- calibrate_varK(targets_b21, targets_b12, K = 2,
                        target_b11_1 = target_b11_1,
                        target_b22_1 = target_b22_1,
                        var1 = 1, var2 = 1, cov12 = 0.1,
                        var_type = "process_error",
                        theta0 = theta_true + 0.02,
                        restarts = 4, maxit = 3000, jitter_sd = 0.03,
                        seed = 3, polish = TRUE, verbose = FALSE)
  expect_s3_class(sol, "varK_calibration")
  expect_lt(sol$loss, 1e-4)
  expect_lt(max(abs(sol$implied["b21", ] - targets_b21)), 5e-3)
  expect_lt(max(abs(sol$implied["b12", ] - targets_b12)), 5e-3)
  expect_lt(abs(sol$implied_b11_1 - target_b11_1), 5e-3)
  expect_lt(abs(sol$implied_b22_1 - target_b22_1), 5e-3)
  expect_equal(sol$params$pe_cov12, 0.1)
  expect_equal(sol$mats$Q[1:2, 1:2], matrix(c(1, 0.1, 0.1, 1), 2, 2))
  expect_equal(dim(sol$implied), c(4L, 3L))
  expect_equal(length(sol$theta), 8L)
})

test_that("process mode hits the specified stationary process covariance", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        target_b11_1 = 0.5, target_b22_1 = 0.5,
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
  # coefficients exactly, so all four targets are exactly attainable; the
  # derived process error covariance at the solution was verified by hand:
  # Q = P - F P F' = [[1.164, -0.214], [-0.214, 0.446]]
  sol <- calibrate_varK(0.2, 0.1, K = 1,
                        target_b11_1 = 0.4, target_b22_1 = 0.6,
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
  expect_equal(sol$implied_b11_1, 0.4, tolerance = 1e-3)
  expect_equal(sol$implied_b22_1, 0.6, tolerance = 1e-3)
  # for K = 1 the free VAR coefficients equal the projection effects
  expect_equal(sol$params$ax[1], 0.4, tolerance = 1e-3)
  expect_equal(sol$params$ay[1], 0.6, tolerance = 1e-3)
  expect_equal(sol$params$pe_var1,  1.164, tolerance = 1e-2)
  expect_equal(sol$params$pe_var2,  0.446, tolerance = 1e-2)
  expect_equal(sol$params$pe_cov12, -0.214, tolerance = 1e-2)
})

test_that("verify_varK works with a process mode calibration", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        target_b11_1 = 0.5, target_b22_1 = 0.5,
                        restarts = 3, maxit = 2000, seed = 2,
                        verbose = FALSE)
  v <- verify_varK(sol, N = 3000, T_obs = 12, seed = 99)
  expect_equal(nrow(v$table), 2L)
  expect_true(all(c("diff_b21_vs_implied", "diff_b12_vs_target")
                  %in% names(v$table)))
  expect_lt(max(abs(v$table$diff_b21_vs_implied)), 0.1)
  expect_lt(max(abs(v$table$diff_b12_vs_implied)), 0.1)
})

test_that("AR targets are process properties, not fixed VAR coefficients (K > 1)", {
  # K = 3: 12 free VAR coefficients, 2 * 4 + 2 = 10 targets
  t21 <- c(0.25, 0.18, 0.12, 0.09)
  t12 <- c(0.15, 0.11, 0.08, 0.05)
  # deliberately start the lag-1 diagonal far away from the targets, so that
  # the result cannot be an artefact of the default start heuristic
  th0 <- c(0.15, 0.05, 0.05,   # ax[1..3]
           0.15, 0.05, 0.05,   # ay[1..3]
           0.25, 0.03, 0.03,   # b21[1..3]
           0.15, 0.03, 0.03)   # b12[1..3]
  sol <- calibrate_varK(t21, t12, K = 3,
                        target_b11_1 = 0.50, target_b22_1 = 0.40,
                        var1 = 1, var2 = 1, cov12 = 0,
                        var_type = "process",
                        theta0 = th0,
                        restarts = 8, maxit = 4000, seed = 11,
                        verbose = FALSE)

  # (a) the lag-1 diagonal elements are genuine free parameters: they are
  #     elements of the optimized vector theta, not values written in
  expect_equal(length(sol$theta), 12L)
  expect_equal(unname(sol$params$ax[1]), unname(sol$theta[1]))
  expect_equal(unname(sol$params$ay[1]), unname(sol$theta[4]))

  # (b) they are not pinned to the targets: the VAR coefficient ends up
  #     further from the target than the implied projection effect does
  expect_false(isTRUE(all.equal(sol$params$ax[1], 0.50, tolerance = 1e-6)))
  expect_false(isTRUE(all.equal(sol$params$ay[1], 0.40, tolerance = 1e-6)))
  expect_gt(abs(sol$params$ax[1] - 0.50), abs(sol$implied_b11_1 - 0.50))
  expect_gt(abs(sol$params$ay[1] - 0.40), abs(sol$implied_b22_1 - 0.40))

  # (c) the higher autoregressive lags are used, i.e. the full VAR(K) shapes
  #     the projection effects
  expect_gt(max(abs(c(sol$params$ax[-1], sol$params$ay[-1]))), 1e-3)

  # (d) nevertheless the implied projection effects at Delta_1 reproduce
  #     the targets, and the cross-lag targets are still matched
  expect_lt(abs(sol$implied_b11_1 - 0.50), 5e-3)
  expect_lt(abs(sol$implied_b22_1 - 0.40), 5e-3)
  expect_lt(max(abs(sol$implied["b21", ] - t21)), 5e-3)
  expect_lt(max(abs(sol$implied["b12", ] - t12)), 5e-3)

  # (e) the process is stable
  expect_lt(max(Mod(eigen(sol$mats$F, only.values = TRUE)$values)), 0.998)

  # (f) verify_varK reports the interval-1 autoregressive comparison
  v <- verify_varK(sol, N = 4000, T_obs = 12, seed = 5)
  expect_true(all(c("target_b11_1", "implied_b11_1", "estimated_b11_1",
                    "target_b22_1", "implied_b22_1", "estimated_b22_1")
                  %in% names(v$ar1)))
  expect_equal(v$ar1$target_b11_1, 0.50)
  expect_equal(v$ar1$implied_b11_1, sol$implied_b11_1)
  expect_lt(abs(v$ar1$estimated_b11_1 - v$ar1$implied_b11_1), 0.1)
  expect_lt(abs(v$ar1$estimated_b22_1 - v$ar1$implied_b22_1), 0.1)
})

test_that("verify_varK evaluates Delta_1 even if it is not requested", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        target_b11_1 = 0.5, target_b22_1 = 0.5,
                        restarts = 3, maxit = 2000, seed = 2,
                        verbose = FALSE)
  v <- verify_varK(sol, N = 3000, T_obs = 12, deltas = 2, seed = 7)
  expect_equal(nrow(v$ar1), 1L)
  expect_true(is.finite(v$ar1$estimated_b11_1))
  expect_true(is.finite(v$ar1$estimated_b22_1))
})

test_that("the loss equals the independently computed sum of squared deviations", {
  # known VAR(2); in "process_error" mode the covariance path is fully
  # determined, so the implied profile can be reconstructed independently
  ax  <- c(0.45, 0.12)
  ay  <- c(0.55, 0.08)
  b21 <- c(0.22, 0.06)
  b12 <- c(0.14, 0.04)
  K <- 2L; delta_max <- 3L
  m <- make_mats_varK(ax, ay, b21, b12)
  P <- stationary_cov(m$F, m$Q)
  prof <- implied_clpm_profile(m$F, P, delta_max)

  # targets deliberately away from the implied values, so every one of the
  # 2 * delta_max + 2 deviations is non-zero
  t21 <- unname(prof["b21", ]) + c(0.010, -0.020, 0.030)
  t12 <- unname(prof["b12", ]) + c(-0.040, 0.050, -0.060)
  tb11 <- unname(prof["b11", 1]) + 0.070
  tb22 <- unname(prof["b22", 1]) - 0.080

  theta <- c(ax, ay, b21, b12)
  got <- varKgen:::loss_varK(theta, K = K, delta_max = delta_max,
                             targets_b21 = t21, targets_b12 = t12,
                             target_b11_1 = tb11, target_b22_1 = tb22,
                             var1 = 1, var2 = 1, cov12 = 0,
                             var_type = "process_error")

  ref <- sum((unname(prof["b21", ]) - t21)^2) +
    sum((unname(prof["b12", ]) - t12)^2) +
    (unname(prof["b11", 1]) - tb11)^2 +
    (unname(prof["b22", 1]) - tb22)^2
  expect_equal(got, ref, tolerance = 1e-12)

  # every single target enters with weight 1 and there is no block
  # weighting: shifting one target by eps changes the loss by exactly
  # (dev + eps)^2 - dev^2 for that one coefficient
  eps <- 0.01
  chk <- function(t21b, t12b, tb11b, tb22b) {
    varKgen:::loss_varK(theta, K = K, delta_max = delta_max,
                        targets_b21 = t21b, targets_b12 = t12b,
                        target_b11_1 = tb11b, target_b22_1 = tb22b,
                        var1 = 1, var2 = 1, cov12 = 0,
                        var_type = "process_error")
  }
  d_ar <- unname(prof["b11", 1]) - tb11
  expect_equal(chk(t21, t12, tb11 - eps, tb22) - got,
               (d_ar + eps)^2 - d_ar^2, tolerance = 1e-12)
  d_cl <- unname(prof["b21", 2]) - t21[2]
  t21c <- t21; t21c[2] <- t21c[2] - eps
  expect_equal(chk(t21c, t12, tb11, tb22) - got,
               (d_cl + eps)^2 - d_cl^2, tolerance = 1e-12)

  # dropping the autoregressive deviations reproduces the cross-lag part
  # alone, i.e. the two AR targets contribute exactly two squared terms
  got_ar0 <- chk(t21, t12, unname(prof["b11", 1]), unname(prof["b22", 1]))
  expect_equal(got - got_ar0,
               (unname(prof["b11", 1]) - tb11)^2 +
                 (unname(prof["b22", 1]) - tb22)^2,
               tolerance = 1e-12)
})

test_that("numerical penalties are separate from the data-fit component", {
  K <- 2L; delta_max <- 2L
  t21 <- c(0.2, 0.1); t12 <- c(0.1, 0.05)
  # unstable process: spectral radius above stab_thresh
  theta_unstable <- c(1.2, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0)
  Fu <- make_mats_varK(c(1.2, 0.0), c(0.5, 0.0), c(0, 0), c(0, 0))$F
  rho <- max(Mod(eigen(Fu, only.values = TRUE)$values))
  expect_gt(rho, 0.998)
  lu <- varKgen:::loss_varK(theta_unstable, K = K, delta_max = delta_max,
                            targets_b21 = t21, targets_b12 = t12,
                            target_b11_1 = 0.5, target_b22_1 = 0.5,
                            var1 = 1, var2 = 1, cov12 = 0,
                            var_type = "process_error")
  expect_gte(lu, 1e9)
  expect_equal(lu, 1e9 + 1e7 * (rho - 0.998), tolerance = 1e-8)

  # non-finite parameters are caught before any covariance work
  theta_bad <- theta_unstable; theta_bad[1] <- NA_real_
  expect_equal(varKgen:::loss_varK(theta_bad, K = K, delta_max = delta_max,
                                   targets_b21 = t21, targets_b12 = t12,
                                   target_b11_1 = 0.5, target_b22_1 = 0.5,
                                   var_type = "process_error"), 1e12)
})

test_that("restart diagnostics and optimizer stages are fully retained", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        target_b11_1 = 0.5, target_b22_1 = 0.4,
                        restarts = 4, maxit = 2000, seed = 3,
                        early_tol = 0, polish = TRUE, verbose = FALSE)

  # one row per restart, with convergence diagnostics
  expect_equal(nrow(sol$restart_log), 4L)
  expect_true(all(c("value", "seconds", "convergence", "fn_evals", "message")
                  %in% names(sol$restart_log)))
  expect_true(all(is.finite(sol$restart_log$fn_evals)))

  # full pre-polish endpoints of every restart
  expect_equal(dim(sol$restart_theta), c(4L, 8L))
  expect_true(all(is.finite(sol$restart_theta)))
  expect_equal(colnames(sol$restart_theta),
               c("ax1", "ax2", "ay1", "ay2", "b211", "b212", "b121", "b122"))

  # the stored endpoints reproduce the logged losses, so every derived
  # quantity (spectral radius, Q eigenvalues, near-best spread) can be
  # reconstructed from them
  for (r in seq_len(nrow(sol$restart_theta))) {
    lr <- varKgen:::loss_varK(sol$restart_theta[r, ], K = 2, delta_max = 2,
                              targets_b21 = c(0.30, 0.20),
                              targets_b12 = c(0.15, 0.12),
                              target_b11_1 = 0.5, target_b22_1 = 0.4,
                              var1 = 1, var2 = 1, cov12 = 0,
                              var_type = "process")
    expect_equal(lr, sol$restart_log$value[r], tolerance = 1e-10)
  }

  # pre-polish and polished solutions are kept apart
  expect_equal(sol$loss_nm, min(sol$restart_log$value))
  expect_true(is.logical(sol$polish_improved))
  if (isTRUE(sol$polish_improved)) {
    expect_equal(sol$loss, sol$loss_polished)
    expect_lt(sol$loss_polished, sol$loss_nm)
    expect_equal(sol$theta, sol$theta_polished)
  } else {
    expect_equal(sol$loss, sol$loss_nm)
    expect_equal(sol$theta, sol$theta_nm)
  }
})

test_that("input validation rejects malformed control arguments", {
  base <- list(targets_b21 = c(0.3, 0.2), targets_b12 = c(0.15, 0.12),
               K = 2, target_b11_1 = 0.5, target_b22_1 = 0.4,
               restarts = 2, maxit = 500, verbose = FALSE)
  bad <- function(...) do.call(calibrate_varK, utils::modifyList(base, list(...)))
  expect_error(bad(early_tol = NA_real_))
  expect_error(bad(restarts = 0))
  expect_error(bad(restarts = 2.5))
  expect_error(bad(maxit = NA_integer_))
  expect_error(bad(jitter_sd = -1))
  expect_error(bad(seed = NA_real_))
  expect_error(bad(polish = NA))
  expect_error(bad(verbose = "yes"))
  expect_error(bad(cores = 0))
})

test_that("polish = FALSE leaves the polish fields empty", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        target_b11_1 = 0.5, target_b22_1 = 0.4,
                        restarts = 3, maxit = 1500, seed = 6,
                        early_tol = 0, polish = FALSE, verbose = FALSE)
  expect_null(sol$theta_polished)
  expect_true(is.na(sol$loss_polished))
  expect_true(is.na(sol$polish_convergence))
  expect_false(sol$polish_improved)
  expect_equal(sol$loss, sol$loss_nm)
  expect_equal(sol$theta, sol$theta_nm)
  expect_true(all(sol$restart_log$status == "completed"))
})

test_that("check_stability = FALSE matches the checked path for stable systems", {
  m <- make_mats_varK(c(0.5, 0.1), c(0.4, 0.05), c(0.2, 0.05), c(0.1, 0.03))
  P_chk <- stationary_cov(m$F, m$Q, check_stability = TRUE)
  P_raw <- stationary_cov(m$F, m$Q, check_stability = FALSE)
  expect_equal(P_chk, P_raw, tolerance = 1e-12)
  pe_chk <- varKgen:::solve_process_error(m$F, 1, 1, 0, check_stability = TRUE)
  pe_raw <- varKgen:::solve_process_error(m$F, 1, 1, 0, check_stability = FALSE)
  expect_equal(pe_chk$Q2, pe_raw$Q2, tolerance = 1e-12)
  expect_equal(pe_chk$P, pe_raw$P, tolerance = 1e-12)
  # the guard only matters for unstable systems
  Fu <- make_mats_varK(c(1.2, 0), c(0.5, 0), c(0, 0), c(0, 0))$F
  expect_null(stationary_cov(Fu, diag(4), check_stability = TRUE))
  expect_error(stationary_cov(m$F, m$Q, check_stability = NA))
})

test_that("default_theta0 validates its own arguments", {
  expect_error(default_theta0(c(2, 3), 0.3, 0.15, 0.5, 0.5))
  expect_error(default_theta0(NA_real_, 0.3, 0.15, 0.5, 0.5))
  expect_error(default_theta0(0, 0.3, 0.15, 0.5, 0.5))
  expect_error(default_theta0(2.5, 0.3, 0.15, 0.5, 0.5))
  expect_error(default_theta0(2, c(0.3, NA), c(0.15, 0.1), 0.5, 0.5))
  expect_error(default_theta0(2, c(0.3, 0.2), c(0.15, Inf), 0.5, 0.5))
  expect_error(default_theta0(2, 0.3, 0.15, NA_real_, 0.5))
})

test_that("seed and cores must be integral, and theta0 is recorded", {
  base <- list(targets_b21 = c(0.3, 0.2), targets_b12 = c(0.15, 0.12),
               K = 2, target_b11_1 = 0.5, target_b22_1 = 0.4,
               restarts = 2, maxit = 500, verbose = FALSE)
  bad <- function(...) do.call(calibrate_varK, utils::modifyList(base, list(...)))
  expect_error(bad(seed = 1.5))
  expect_error(bad(cores = 2.5))

  th0 <- c(0.30, 0.05, 0.35, 0.05, 0.25, 0.03, 0.12, 0.03)
  sol <- do.call(calibrate_varK, utils::modifyList(base, list(theta0 = th0)))
  expect_equal(sol$settings$theta0, th0)
  expect_equal(sol$settings$cores_used, 1L)
  expect_null(sol$settings$cores_requested)
})

test_that("the parallel path reproduces the sequential one", {
  skip_on_cran()
  skip_if_not_installed("parallel")
  args <- list(targets_b21 = c(0.30, 0.20), targets_b12 = c(0.15, 0.12),
               K = 2, target_b11_1 = 0.5, target_b22_1 = 0.4,
               restarts = 4, maxit = 1500, seed = 8, early_tol = 0,
               verbose = FALSE)
  seq_sol <- do.call(calibrate_varK, args)
  # skip only on a narrow, explicit environment check; from here on any
  # error is a real defect and must fail the test rather than skip it
  cluster_ok <- tryCatch({
    cl <- parallel::makeCluster(2L)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterCall(cl, function(p) .libPaths(p), .libPaths())
    ok <- unlist(parallel::clusterEvalQ(
      cl, requireNamespace("varKgen", quietly = TRUE)))
    parallel::stopCluster(cl)
    isTRUE(all(ok))
  }, error = function(e) FALSE)
  skip_if_not(cluster_ok,
              "workers cannot load varKgen in this environment")
  par_sol <- do.call(calibrate_varK, c(args, list(parallel = TRUE, cores = 2)))
  # identical starts imply identical per-restart results
  expect_equal(par_sol$restart_log$value, seq_sol$restart_log$value,
               tolerance = 1e-10)
  expect_equal(par_sol$restart_theta, seq_sol$restart_theta,
               tolerance = 1e-10)
  expect_equal(par_sol$loss_nm, seq_sol$loss_nm, tolerance = 1e-10)
  expect_equal(par_sol$theta_nm, seq_sol$theta_nm, tolerance = 1e-10)
  expect_equal(par_sol$settings$cores_used, 2L)
})

test_that("run_restart classifies failures instead of aborting", {
  rr <- varKgen:::run_restart
  # an objective that always raises an error
  r1 <- rr(c(0.1, 0.2), function(theta) stop("deliberate failure"), maxit = 10)
  expect_equal(r1$status, "error")
  expect_equal(r1$value, Inf)
  expect_true(is.na(r1$convergence))
  expect_true(is.na(r1$fn_evals))
  expect_match(r1$message, "deliberate failure")
  expect_equal(length(r1$par), 2L)
  expect_true(all(is.na(r1$par)))

  # an objective that is not evaluable at the start
  r2 <- rr(c(0.1, 0.2), function(theta) NaN, maxit = 10)
  expect_true(r2$status %in% c("error", "invalid"))
  expect_equal(r2$value, Inf)
  expect_true(all(is.na(r2$par)))

  # a well-behaved objective completes
  r3 <- rr(c(0.1, 0.2), function(theta) sum((theta - 1)^2), maxit = 500)
  expect_equal(r3$status, "completed")
  expect_true(is.finite(r3$value))
  expect_true(all(is.finite(r3$par)))
  expect_true(is.finite(r3$fn_evals))
})

test_that("a run in which all restarts fail stops with a clear error", {
  skip_if_not(utils::packageVersion("testthat") >= "3.1.8",
              "needs local_mocked_bindings()")
  testthat::local_mocked_bindings(
    run_restart = function(theta_start, fn, maxit, fn_args = list()) {
      list(value = Inf, par = rep(NA_real_, length(theta_start)),
           seconds = 0, convergence = NA_integer_, fn_evals = NA_integer_,
           message = "forced failure", status = "error")
    },
    .package = "varKgen")
  expect_error(
    calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                   target_b11_1 = 0.5, target_b22_1 = 0.4,
                   restarts = 3, maxit = 500, verbose = FALSE),
    "All optimizer restarts failed")
})

test_that("guards on target scale, early_tol and argument types hold", {
  base <- list(targets_b21 = c(0.3, 0.2), targets_b12 = c(0.15, 0.12),
               K = 2, target_b11_1 = 0.5, target_b22_1 = 0.4,
               restarts = 2, maxit = 500, verbose = FALSE)
  bad <- function(...) do.call(calibrate_varK, utils::modifyList(base, list(...)))
  # targets so large that a valid fit could reach the penalty sentinel
  expect_error(bad(targets_b21 = c(1e4, 1e4)), "too large for the penalty")
  # early stopping must not be satisfiable by a penalty solution
  expect_error(bad(early_tol = 1e10))
  # logicals must not pass as numeric scalars
  expect_error(bad(restarts = TRUE))
  expect_error(bad(maxit = TRUE))
  expect_error(bad(K = TRUE))
})

test_that("theta vectors carry canonical names and polish is fully reported", {
  sol <- calibrate_varK(c(0.30, 0.20), c(0.15, 0.12), K = 2,
                        target_b11_1 = 0.5, target_b22_1 = 0.4,
                        restarts = 3, maxit = 1500, seed = 9,
                        early_tol = 0, verbose = FALSE)
  nms <- c("ax1", "ax2", "ay1", "ay2", "b211", "b212", "b121", "b122")
  expect_equal(names(sol$theta), nms)
  expect_equal(names(sol$theta_nm), nms)
  if (!is.null(sol$theta_polished)) expect_equal(names(sol$theta_polished), nms)
  expect_true(sol$polish_attempted)
  expect_equal(sol$polish_status, "completed")
  expect_true(is.finite(sol$polish_seconds))
  expect_true(is.finite(sol$polish_fn_evals))
  # provenance for reproducibility
  expect_true(is.character(sol$settings$package_version))
  expect_true(is.character(sol$settings$R_version))
  expect_length(sol$settings$RNGkind, 3L)
})
