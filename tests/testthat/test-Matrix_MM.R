## ---- matrix_mm_fit (EM mode) ----

test_that("mm_fit EM mode returns correct structure without noise", {
  set.seed(42)
  x_list <- c(
    lapply(1:10, function(i) matrix(rnorm(6, mean = 3), 2, 3)),
    lapply(1:10, function(i) matrix(rnorm(6, mean = -3), 2, 3))
  )
  fit <- matrix_mm_fit(x_list, g = 2, method = "em", noise_type = NULL,
                        max_iter = 20, verbose = FALSE)

  expect_type(fit, "list")
  expect_equal(fit$method, "em")
  expect_length(fit$pi, 2)
  expect_length(fit$cluster, 20)
  expect_true(all(fit$cluster %in% 1:2))
  expect_false("noise" %in% names(fit))
})

test_that("mm_fit MM mode returns correct structure without noise", {
  set.seed(42)
  x_list <- c(
    lapply(1:10, function(i) matrix(rnorm(6, mean = 3), 2, 3)),
    lapply(1:10, function(i) matrix(rnorm(6, mean = -3), 2, 3))
  )
  fit <- matrix_mm_fit(x_list, g = 2, method = "mm", noise_type = NULL,
                        max_iter = 20, verbose = FALSE)

  expect_type(fit, "list")
  expect_equal(fit$method, "mm")
  expect_length(fit$pi, 2)
  expect_length(fit$cluster, 20)
})

test_that("mm_fit EM with HC noise returns noise info", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- matrix_mm_fit(x_list, g = 2, method = "em", noise_type = "hc",
                        max_iter = 20, verbose = FALSE)

  expect_true("noise" %in% names(fit))
  expect_equal(fit$noise$type, "hc")
  expect_length(fit$pi, 3)  # 2 components + noise
  expect_true(all(fit$cluster %in% 0:2))
})

test_that("mm_fit MM with HC noise returns noise info", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- matrix_mm_fit(x_list, g = 2, method = "mm", noise_type = "hc",
                        max_iter = 20, verbose = FALSE)

  expect_true("noise" %in% names(fit))
  expect_equal(fit$noise$type, "hc")
})

test_that("mm_fit mixing proportions sum to 1", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- matrix_mm_fit(x_list, g = 2, method = "mm", noise_type = NULL,
                        max_iter = 20, verbose = FALSE)
  expect_equal(sum(fit$pi), 1, tolerance = 1e-10)
})

test_that("mm_fit responsibilities sum to 1 per row", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- matrix_mm_fit(x_list, g = 2, method = "em", noise_type = NULL,
                        max_iter = 20, verbose = FALSE)
  row_sums <- rowSums(fit$z)
  expect_equal(row_sums, rep(1, 15), tolerance = 1e-8)
})

test_that("mm_fit recovers well-separated clusters (EM)", {
  set.seed(123)
  x_list <- c(
    lapply(1:15, function(i) matrix(rnorm(6, mean = 10, sd = 0.1), 2, 3)),
    lapply(1:15, function(i) matrix(rnorm(6, mean = -10, sd = 0.1), 2, 3))
  )
  fit <- matrix_mm_fit(x_list, g = 2, method = "em", noise_type = NULL,
                        max_iter = 50, verbose = FALSE)
  group1 <- fit$cluster[1:15]
  group2 <- fit$cluster[16:30]
  expect_true(length(unique(group1)) == 1)
  expect_true(length(unique(group2)) == 1)
  expect_true(group1[1] != group2[1])
})

test_that("mm_fit recovers well-separated clusters (MM)", {
  set.seed(123)
  x_list <- c(
    lapply(1:15, function(i) matrix(rnorm(6, mean = 10, sd = 0.1), 2, 3)),
    lapply(1:15, function(i) matrix(rnorm(6, mean = -10, sd = 0.1), 2, 3))
  )
  fit <- matrix_mm_fit(x_list, g = 2, method = "mm", noise_type = NULL,
                        max_iter = 50, verbose = FALSE)
  group1 <- fit$cluster[1:15]
  group2 <- fit$cluster[16:30]
  expect_true(length(unique(group1)) == 1)
  expect_true(length(unique(group2)) == 1)
  expect_true(group1[1] != group2[1])
})

test_that("mm_fit with BR noise works", {
  set.seed(42)
  # BR needs enough data for convex hull
  x_list <- lapply(1:20, function(i) matrix(rnorm(4, sd = 2), 2, 2))
  fit <- matrix_mm_fit(x_list, g = 2, method = "em", noise_type = "br",
                        max_iter = 20, verbose = FALSE)
  expect_true("noise" %in% names(fit))
  expect_equal(fit$noise$type, "br")
})

test_that("mm_fit works with different init methods", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))

  fit_kmeans <- matrix_mm_fit(x_list, g = 2, method = "em", noise_type = NULL,
                               init = "kmeans", max_iter = 10, verbose = FALSE)
  fit_ecme <- matrix_mm_fit(x_list, g = 2, method = "em", noise_type = NULL,
                             init = "ecme", max_iter = 10, verbose = FALSE)

  expect_length(fit_kmeans$cluster, 15)
  expect_length(fit_ecme$cluster, 15)
})
