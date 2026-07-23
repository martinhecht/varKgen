test_that("fit_clpm_ols matches lm()", {
  set.seed(1)
  n <- 200
  df <- data.frame(Y1_p = rnorm(n), Y2_p = rnorm(n))
  df$Y1_q <- 0.5 * df$Y1_p + 0.2 * df$Y2_p + rnorm(n)
  df$Y2_q <- 0.1 * df$Y1_p + 0.4 * df$Y2_p + rnorm(n)
  co <- fit_clpm_ols(df)
  f1 <- lm(Y1_q ~ Y1_p + Y2_p, data = df)
  f2 <- lm(Y2_q ~ Y2_p + Y1_p, data = df)
  expect_equal(unname(co["b11"]), unname(coef(f1)["Y1_p"]), tolerance = 1e-10)
  expect_equal(unname(co["b21"]), unname(coef(f1)["Y2_p"]), tolerance = 1e-10)
  expect_equal(unname(co["b22"]), unname(coef(f2)["Y2_p"]), tolerance = 1e-10)
  expect_equal(unname(co["b12"]), unname(coef(f2)["Y1_p"]), tolerance = 1e-10)
})

test_that("make_delta_pairs_overlapping stacks correctly and validates delta", {
  Y1 <- matrix(1:6, nrow = 2)   # 2 persons, 3 time points
  Y2 <- Y1 + 100
  d1 <- make_delta_pairs_overlapping(Y1, Y2, 1)
  expect_equal(nrow(d1), 4L)
  expect_equal(d1$Y1_p, c(1, 2, 3, 4))
  expect_equal(d1$Y1_q, c(3, 4, 5, 6))
  expect_equal(d1$Y2_p, c(101, 102, 103, 104))
  expect_equal(d1$id, c(1, 2, 1, 2))
  expect_equal(d1$p, c(1, 1, 2, 2))
  expect_equal(d1$q, c(2, 2, 3, 3))
  d2 <- make_delta_pairs_overlapping(Y1, Y2, 2)
  expect_equal(nrow(d2), 2L)
  expect_error(make_delta_pairs_overlapping(Y1, Y2, 3))
})
