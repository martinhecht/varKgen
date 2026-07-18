test_that("stationary_cov solves the discrete Lyapunov equation", {
  m <- make_mats_varK(ax = c(0.5, 0.10), ay = c(0.4, 0.05),
                      b21 = c(0.2, 0.05), b12 = c(0.1, 0.02),
                      pe_cov12 = 0.3)
  P <- stationary_cov(m$F, m$Q)
  expect_false(is.null(P))
  expect_lt(max(abs(P - m$F %*% P %*% t(m$F) - m$Q)), 1e-9)
  expect_equal(P, t(P))
})

test_that("kron and iterate methods agree", {
  m <- make_mats_varK(ax = c(0.6, 0.05), ay = c(0.5, 0.10),
                      b21 = c(0.15, 0.05), b12 = c(0.10, 0.02))
  P_kron <- stationary_cov(m$F, m$Q, method = "kron")
  P_iter <- stationary_cov(m$F, m$Q, method = "iterate")
  expect_false(is.null(P_kron))
  expect_false(is.null(P_iter))
  expect_lt(max(abs(P_kron - P_iter)), 1e-8)
})

test_that("stationary_cov returns NULL for an unstable system", {
  Fm <- matrix(c(1.05, 0, 0, 0.5), 2, 2)
  Q  <- diag(2)
  expect_null(stationary_cov(Fm, Q, method = "kron"))
})

test_that("solve_process_error reproduces the target process covariance", {
  m <- make_mats_varK(ax = c(0.5, 0.10), ay = c(0.4, 0.05),
                      b21 = c(0.2, 0.05), b12 = c(0.1, 0.02))
  Fm <- m$F
  res <- varKgen:::solve_process_error(Fm, 1.3, 0.9, 0.2)
  expect_false(is.null(res))
  # the derived process error covariance solves the Lyapunov equation ...
  Qf <- matrix(0, nrow(Fm), ncol(Fm))
  Qf[1:2, 1:2] <- res$Q2
  expect_lt(max(abs(res$P - Fm %*% res$P %*% t(Fm) - Qf)), 1e-8)
  # ... hits the target block exactly ...
  expect_equal(res$P[1, 1], 1.3, tolerance = 1e-8)
  expect_equal(res$P[2, 2], 0.9, tolerance = 1e-8)
  expect_equal(res$P[1, 2], 0.2, tolerance = 1e-8)
  # ... agrees with the direct Lyapunov solver ...
  P2 <- stationary_cov(Fm, Qf)
  expect_lt(max(abs(res$P - P2)), 1e-8)
  # ... and is positive semi-definite here
  expect_gte(res$min_eig, 0)
  expect_equal(res$min_eig,
               min(eigen(res$Q2, only.values = TRUE)$values),
               tolerance = 1e-10)
})

test_that("solve_process_error returns NULL for an unstable system", {
  expect_null(varKgen:::solve_process_error(matrix(c(1.05, 0, 0, 0.5), 2, 2),
                                            1, 1, 0))
})
