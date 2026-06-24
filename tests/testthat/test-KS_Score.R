## ---- matrix_noise_ks_score ----

test_that("ks_score returns correct structure", {
  set.seed(42)
  x_list <- c(
    lapply(1:15, function(i) matrix(rnorm(6, mean = 3), 2, 3)),
    lapply(1:15, function(i) matrix(rnorm(6, mean = -3), 2, 3))
  )
  fit <- matrix_variate_mixture_fit(x_list, g = 2, max_iter = 30, verbose = FALSE)
  ks <- matrix_noise_ks_score(fit, x_list)

  expect_type(ks, "list")
  expect_named(ks, c("statistic", "p.value", "n_used"))
  expect_true(is.numeric(ks$statistic))
  expect_true(is.numeric(ks$p.value))
  expect_true(is.numeric(ks$n_used))
})

test_that("ks_score statistic is between 0 and 1", {
  set.seed(42)
  x_list <- lapply(1:20, function(i) matrix(rnorm(6), 2, 3))
  fit <- matrix_variate_mixture_fit(x_list, g = 2, max_iter = 20, verbose = FALSE)
  ks <- matrix_noise_ks_score(fit, x_list)
  expect_gte(ks$statistic, 0)
  expect_lte(ks$statistic, 1)
})

test_that("ks_score handles fit with all noise (cluster = 0)", {
  set.seed(42)
  x_list <- lapply(1:10, function(i) matrix(rnorm(6), 2, 3))
  # Fake a fit where all points are noise
  fake_fit <- list(
    cluster = rep(0L, 10),
    M = list(matrix(0, 2, 3)),
    U = list(diag(2)),
    V = list(diag(3))
  )
  ks <- matrix_noise_ks_score(fake_fit, x_list)
  expect_equal(ks$statistic, Inf)
  expect_true(is.na(ks$p.value))
})

test_that("ks_score n_used matches non-noise count", {
  set.seed(42)
  x_list <- lapply(1:20, function(i) matrix(rnorm(6), 2, 3))
  fit <- matrix_variate_mixture_fit(x_list, g = 2, max_iter = 20, verbose = FALSE)
  ks <- matrix_noise_ks_score(fit, x_list)
  expect_equal(ks$n_used, sum(fit$cluster > 0))
})

## ---- matrix_noise_hc_heuristic_grid ----

test_that("heuristic grid returns sorted positive values", {
  set.seed(42)
  x_list <- lapply(1:10, function(i) matrix(rnorm(6), 2, 3))
  grid <- matrix_noise_hc_heuristic_grid(x_list)
  expect_true(all(grid > 0))
  expect_true(all(is.finite(grid)))
  expect_equal(grid, sort(grid))
})

test_that("heuristic grid respects n_points", {
  set.seed(42)
  x_list <- lapply(1:10, function(i) matrix(rnorm(6), 2, 3))
  grid <- matrix_noise_hc_heuristic_grid(x_list, n_points = 15)
  expect_lte(length(grid), 15)
  expect_gte(length(grid), 2)
})

test_that("heuristic grid adapts to matrix dimension", {
  set.seed(42)
  x_small <- lapply(1:10, function(i) matrix(rnorm(4), 2, 2))
  x_large <- lapply(1:10, function(i) matrix(rnorm(20), 4, 5))
  grid_small <- matrix_noise_hc_heuristic_grid(x_small)
  grid_large <- matrix_noise_hc_heuristic_grid(x_large)
  # Larger dimension should push grid toward smaller k values
  expect_lt(median(grid_large), median(grid_small))
})
