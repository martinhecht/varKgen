test_that("companion matrix has correct structure for K = 1 (bug fix)", {
  m <- make_mats_varK(ax = 0.5, ay = 0.4, b21 = 0.2, b12 = 0.1)
  expect_equal(dim(m$F), c(2L, 2L))
  expect_equal(m$F, matrix(c(0.5, 0.1, 0.2, 0.4), 2, 2))
  expect_equal(m$Q, diag(2))
})

test_that("companion matrix has correct structure for K = 3", {
  ax  <- c(0.5, 0.10, 0.05)
  ay  <- c(0.4, 0.08, 0.02)
  b21 <- c(0.3, 0.20, 0.10)
  b12 <- c(0.15, 0.10, 0.05)
  m <- make_mats_varK(ax, ay, b21, b12,
                      pe_var1 = 2, pe_var2 = 3, pe_cov12 = 0.5)
  Fm <- m$F
  expect_equal(dim(Fm), c(6L, 6L))
  for (k in 1:3) {
    expect_equal(Fm[1, 2 * (k - 1) + 1], ax[k])
    expect_equal(Fm[1, 2 * (k - 1) + 2], b21[k])
    expect_equal(Fm[2, 2 * (k - 1) + 1], b12[k])
    expect_equal(Fm[2, 2 * (k - 1) + 2], ay[k])
  }
  # shift registers
  expect_equal(Fm[3, 1], 1)
  expect_equal(Fm[4, 2], 1)
  expect_equal(Fm[5, 3], 1)
  expect_equal(Fm[6, 4], 1)
  expect_equal(sum(Fm[3:6, ] != 0), 4)
  # process error covariance only in the top-left block
  expect_equal(m$Q[1:2, 1:2], matrix(c(2, 0.5, 0.5, 3), 2, 2))
  expect_true(all(m$Q[3:6, ] == 0))
  expect_true(all(m$Q[, 3:6] == 0))
})

test_that("input validation works", {
  expect_error(make_mats_varK(ax = c(0.5, 0.1), ay = 0.4,
                              b21 = c(0.2, 0.1), b12 = c(0.1, 0.05)))
  expect_error(make_mats_varK(ax = 0.5, ay = 0.4, b21 = 0.2, b12 = 0.1,
                              pe_cov12 = 2))
})
