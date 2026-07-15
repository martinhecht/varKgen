test_that("stationary_cov solves the discrete Lyapunov equation", {
  m <- make_mats_varK(ax = c(0.5, 0.10), ay = c(0.4, 0.05),
                      b21 = c(0.2, 0.05), b12 = c(0.1, 0.02),
                      cov_ue = 0.3)
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
