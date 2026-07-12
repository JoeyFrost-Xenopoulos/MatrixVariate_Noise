test_that("mv_validate_x_list accepts valid input", {
  x1 <- matrix(1:6, nrow = 2, ncol = 3)
  x2 <- matrix(7:12, nrow = 2, ncol = 3)
  result <- mv_validate_x_list(list(x1, x2))
  expect_identical(result, list(x1, x2))
})

test_that("mv_validate_x_list rejects empty list", {
  expect_error(mv_validate_x_list(list()), "non-empty list")
})

test_that("mv_validate_x_list rejects non-list", {
  expect_error(mv_validate_x_list(matrix(1:4, 2, 2)), "non-empty list")
})

test_that("mv_validate_x_list rejects dimension mismatch", {
  x1 <- matrix(1:6, nrow = 2, ncol = 3)
  x2 <- matrix(1:4, nrow = 2, ncol = 2)
  expect_error(mv_validate_x_list(list(x1, x2)), "dimensions")
})

test_that("mv_validate_x_list rejects non-matrix element", {
  x1 <- matrix(1:6, nrow = 2, ncol = 3)
  expect_error(mv_validate_x_list(list(x1, c(1, 2, 3))), "not a matrix")
})

test_that("mv_log_sum_exp basic computation", {
  # log(exp(1) + exp(2)) = log(e + e^2) = 2 + log(1 + exp(-1))
  result <- mv_log_sum_exp(c(1, 2))
  expected <- log(exp(1) + exp(2))
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("mv_log_sum_exp handles single value", {
  expect_equal(mv_log_sum_exp(5), 5, tolerance = 1e-10)
})

test_that("mv_log_sum_exp handles large values without overflow", {
  # Large values that would overflow with naive exp()
  vals <- c(1000, 1001, 1002)
  result <- mv_log_sum_exp(vals)
  expected <- 1002 + log(exp(-2) + exp(-1) + 1)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("mv_log_sum_exp handles all -Inf", {
  expect_equal(mv_log_sum_exp(c(-Inf, -Inf)), -Inf)
})

test_that("mv_log_sum_exp handles mix of finite and -Inf", {
  result <- mv_log_sum_exp(c(-Inf, 3, -Inf))
  expect_equal(result, 3, tolerance = 1e-10)
})

test_that("mv_log_sum_exp handles empty finite after filtering", {
  expect_equal(mv_log_sum_exp(c(-Inf)), -Inf)
})
