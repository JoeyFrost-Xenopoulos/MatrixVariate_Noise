## ---- diagnose_kmeans_wcss ----

test_that("diagnose_kmeans_wcss works with matrix data and g", {
  set.seed(42)
  data <- matrix(rnorm(100), nrow = 20, ncol = 5)
  result <- diagnose_kmeans_wcss(data, g = 2, verbose = FALSE)

  expect_type(result, "list")
  expect_named(result, c("A", "B", "recommendation", "message"))
  expect_true(is.numeric(result$A))
  expect_true(is.na(result$B))
  expect_equal(result$recommendation, "unknown")
})

test_that("diagnose_kmeans_wcss with true_labels computes B", {
  set.seed(42)
  data <- rbind(
    matrix(rnorm(50, mean = 5), nrow = 10, ncol = 5),
    matrix(rnorm(50, mean = -5), nrow = 10, ncol = 5)
  )
  true_labels <- rep(1:2, each = 10)
  result <- diagnose_kmeans_wcss(data, true_labels = true_labels, g = 2, verbose = FALSE)
  expect_true(is.numeric(result$B))
  expect_true(result$recommendation %in% c("try_kmeanspp", "euclidean_insufficient"))
})

test_that("diagnose_kmeans_wcss works with list of matrices", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  result <- diagnose_kmeans_wcss(x_list, g = 2, verbose = FALSE)
  expect_true(is.numeric(result$A))
})

test_that("diagnose_kmeans_wcss with kmeans_labels computes WCSS manually", {
  set.seed(42)
  data <- matrix(rnorm(60), nrow = 12, ncol = 5)
  labels <- rep(1:2, each = 6)
  result <- diagnose_kmeans_wcss(data, kmeans_labels = labels, verbose = FALSE)
  expect_true(is.numeric(result$A))
  expect_gt(result$A, 0)
})

test_that("diagnose_kmeans_wcss errors on length mismatch", {
  data <- matrix(rnorm(40), 10, 4)
  expect_error(
    diagnose_kmeans_wcss(data, kmeans_labels = 1:5, verbose = FALSE),
    "length must match"
  )
  expect_error(
    diagnose_kmeans_wcss(data, true_labels = 1:5, g = 2, verbose = FALSE),
    "length must match"
  )
})

test_that("diagnose_kmeans_wcss errors on list with dimension mismatch", {
  x1 <- matrix(1, 2, 3)
  x2 <- matrix(1, 3, 3)
  expect_error(
    diagnose_kmeans_wcss(list(x1, x2), g = 2, verbose = FALSE),
    "same dimensions"
  )
})

test_that("diagnose_kmeans_wcss errors on list with non-matrices", {
  expect_error(
    diagnose_kmeans_wcss(list(c(1, 2), c(3, 4)), g = 2, verbose = FALSE),
    "must be matrices"
  )
})

test_that("diagnose_kmeans_wcss errors when no g and no kmeans_labels", {
  data <- matrix(rnorm(20), 5, 4)
  expect_error(diagnose_kmeans_wcss(data, verbose = FALSE), "either.*kmeans_labels.*or.*g")
})

test_that("diagnose_kmeans_wcss recommendation is try_kmeanspp when B < A", {
  set.seed(1)
  # Fabricate a scenario: true labels should give lower WCSS
  data <- rbind(
    matrix(rnorm(100, mean = 10, sd = 0.1), nrow = 20, ncol = 5),
    matrix(rnorm(100, mean = -10, sd = 0.1), nrow = 20, ncol = 5)
  )
  true_labels <- rep(1:2, each = 20)
  # Intentionally bad kmeans_labels
  bad_labels <- rep(c(1, 2), 20)
  result <- diagnose_kmeans_wcss(data, true_labels = true_labels, kmeans_labels = bad_labels, verbose = FALSE)
  expect_equal(result$recommendation, "try_kmeanspp")
})

## ---- matrix_noise_ecdf_vs_cdf_plot ----

test_that("ecdf_vs_cdf_plot errors with too few non-noise points", {
  set.seed(42)
  x_list <- lapply(1:10, function(i) matrix(rnorm(6), 2, 3))
  fake_fit <- list(
    cluster = rep(0L, 10),
    M = list(matrix(0, 2, 3)),
    U = list(diag(2)),
    V = list(diag(3))
  )
  expect_error(
    matrix_noise_ecdf_vs_cdf_plot(fake_fit, x_list),
    "Not enough non-noise"
  )
})

test_that("ecdf_vs_cdf_plot runs without error on valid fit", {
  set.seed(42)
  x_list <- c(
    lapply(1:15, function(i) matrix(rnorm(6, mean = 3), 2, 3)),
    lapply(1:15, function(i) matrix(rnorm(6, mean = -3), 2, 3))
  )
  fit <- matrix_variate_mixture_fit(x_list, g = 2, max_iter = 30, verbose = FALSE)
  # Should produce a plot without error (warnings from plot devices are OK)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(
    suppressWarnings(matrix_noise_ecdf_vs_cdf_plot(fit, x_list))
  )
})
