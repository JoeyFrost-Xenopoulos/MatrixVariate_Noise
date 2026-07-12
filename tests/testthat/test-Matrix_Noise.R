## ---- mv_noise_fit (HC) ----

test_that("noise_fit HC returns correct structure", {
  set.seed(42)
  x_list <- c(
    lapply(1:10, function(i) matrix(rnorm(6, mean = 3), 2, 3)),
    lapply(1:10, function(i) matrix(rnorm(6, mean = -3), 2, 3))
  )
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                                   max_iter = 30, verbose = FALSE)

  expect_type(fit, "list")
  expect_true("noise" %in% names(fit))
  expect_equal(fit$noise$type, "hc")
  expect_true(is.numeric(fit$noise$k))
  expect_length(fit$pi, 3)  # 2 components + noise
  expect_true(all(fit$cluster %in% 0:2))
  expect_length(fit$cluster, 20)
})

test_that("noise_fit HC mixing proportions sum to 1", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                                   max_iter = 20, verbose = FALSE)
  expect_equal(sum(fit$pi), 1, tolerance = 1e-10)
})

test_that("noise_fit HC cluster 0 represents noise", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                                   max_iter = 20, verbose = FALSE)
  # cluster=0 is valid for noise

  expect_true(all(fit$cluster %in% 0:2))
})

test_that("noise_fit HC responsibilities sum to 1 per row", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                                   max_iter = 20, verbose = FALSE)
  row_sums <- rowSums(fit$z)
  expect_equal(row_sums, rep(1, 15), tolerance = 1e-8)
})

test_that("noise_fit HC log-likelihood trace is non-decreasing", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                                   max_iter = 30, verbose = FALSE)
  ll <- fit$logLik
  if (length(ll) > 1) {
    diffs <- diff(ll)
    expect_true(all(diffs >= -1e-4))
  }
})

## ---- mv_noise_fit (BR) ----

test_that("noise_fit BR returns correct structure", {
  set.seed(42)
  # BR needs enough data for convex hull in 2*2=4 dimensions
  x_list <- lapply(1:20, function(i) matrix(rnorm(4, sd = 2), 2, 2))
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "br",
                                   max_iter = 20, verbose = FALSE)

  expect_type(fit, "list")
  expect_equal(fit$noise$type, "br")
  expect_true(all(fit$cluster %in% 0:2))
})

## ---- noise_fit with init options ----

test_that("noise_fit works with kmeans init", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                                   max_iter = 10, init = "kmeans", verbose = FALSE)
  expect_length(fit$cluster, 15)
})

test_that("noise_fit works with emrefine init", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                                   max_iter = 10, init = "emrefine", verbose = FALSE)
  expect_length(fit$cluster, 15)
})

test_that("noise_fit works with dbscan init", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                                   max_iter = 10, init = "dbscan", verbose = FALSE)
  expect_length(fit$cluster, 15)
})
