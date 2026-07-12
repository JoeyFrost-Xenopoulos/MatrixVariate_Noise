## ---- make_spd ----

test_that("make_spd returns positive definite matrix for SPD input", {
  mat <- diag(3)
  result <- make_spd(mat)
  expect_true(is.matrix(result))
  # Cholesky should succeed
  expect_silent(chol(result))
})

test_that("make_spd symmetrizes asymmetric input", {
  mat <- matrix(c(4, 1, 2, 5), 2, 2)
  result <- make_spd(mat)
  expect_equal(result, t(result), tolerance = 1e-10)
})

test_that("make_spd rescues nearly singular matrix", {
  # Create a near-singular matrix
  mat <- matrix(c(1, 0.9999999, 0.9999999, 1), 2, 2)
  result <- make_spd(mat)
  expect_silent(chol(result))
})

test_that("make_spd errors on hopeless matrix", {
  # Negative definite matrix with tiny jitter budget
  mat <- -diag(3) * 1e12
  expect_error(make_spd(mat, jitter = 1e-8, max_tries = 2), "positive definite")
})

## ---- mv_mahalanobis ----

test_that("mv_mahalanobis is zero at the mean", {
  mean_mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  row_cov <- diag(2)
  col_cov <- diag(3)
  d <- mv_mahalanobis(mean_mat, mean_mat, row_cov, col_cov)
  expect_equal(d, 0, tolerance = 1e-8)
})

test_that("mv_mahalanobis is positive away from mean", {
  mean_mat <- matrix(0, 2, 3)
  x <- matrix(1, 2, 3)
  row_cov <- diag(2)
  col_cov <- diag(3)
  d <- mv_mahalanobis(x, mean_mat, row_cov, col_cov)
  expect_gt(d, 0)
})

test_that("mv_mahalanobis with identity cov equals Frobenius norm squared", {
  set.seed(42)
  mean_mat <- matrix(0, 2, 3)
  x <- matrix(rnorm(6), 2, 3)
  row_cov <- diag(2)
  col_cov <- diag(3)
  d <- mv_mahalanobis(x, mean_mat, row_cov, col_cov)
  frobenius_sq <- sum(x^2)
  expect_equal(d, frobenius_sq, tolerance = 1e-6)
})

test_that("mv_mahalanobis scales with covariance", {
  set.seed(123)
  x <- matrix(rnorm(6), 2, 3)
  mean_mat <- matrix(0, 2, 3)

  d1 <- mv_mahalanobis(x, mean_mat, diag(2), diag(3))
  d2 <- mv_mahalanobis(x, mean_mat, 2 * diag(2), diag(3))
  # Larger row covariance => smaller distance
  expect_lt(d2, d1)
})

## ---- mv_log_density ----

test_that("mv_log_density returns finite scalar", {
  set.seed(1)
  x <- matrix(rnorm(6), 2, 3)
  mean_mat <- matrix(0, 2, 3)
  row_cov <- diag(2)
  col_cov <- diag(3)
  ld <- mv_log_density(x, mean_mat, row_cov, col_cov)
  expect_true(is.finite(ld))
  expect_length(ld, 1)
})

test_that("mv_log_density is maximized at the mean", {
  set.seed(10)
  mean_mat <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
  row_cov <- diag(2)
  col_cov <- diag(3)
  ld_at_mean <- mv_log_density(mean_mat, mean_mat, row_cov, col_cov)
  ld_away <- mv_log_density(mean_mat + 5, mean_mat, row_cov, col_cov)
  expect_gt(ld_at_mean, ld_away)
})

test_that("mv_log_density is always negative", {
  set.seed(5)
  x <- matrix(rnorm(6), 2, 3)
  mean_mat <- matrix(0, 2, 3)
  ld <- mv_log_density(x, mean_mat, diag(2), diag(3))
  expect_lt(ld, 0)
})

## ---- mv_mixture_fit ----

test_that("mv_mixture_fit returns correct structure", {
  set.seed(42)
  # Create simple 2-component data
  x_list <- c(
    lapply(1:10, function(i) matrix(rnorm(6, mean = 3), 2, 3)),
    lapply(1:10, function(i) matrix(rnorm(6, mean = -3), 2, 3))
  )
  fit <- mv_mixture_fit(x_list, g = 2, max_iter = 20, verbose = FALSE)

  expect_type(fit, "list")
  expect_named(fit, c("pi", "M", "U", "V", "z", "cluster", "logLik", "iterations", "converged"))
  expect_length(fit$pi, 2)
  expect_length(fit$M, 2)
  expect_length(fit$U, 2)
  expect_length(fit$V, 2)
  expect_equal(nrow(fit$z), 20)
  expect_equal(ncol(fit$z), 2)
  expect_length(fit$cluster, 20)
  expect_true(all(fit$cluster %in% 1:2))
  expect_true(length(fit$logLik) > 0)
})

test_that("mv_mixture_fit mixing proportions sum to 1", {
  set.seed(42)
  x_list <- c(
    lapply(1:10, function(i) matrix(rnorm(6, mean = 3), 2, 3)),
    lapply(1:10, function(i) matrix(rnorm(6, mean = -3), 2, 3))
  )
  fit <- mv_mixture_fit(x_list, g = 2, max_iter = 20, verbose = FALSE)
  expect_equal(sum(fit$pi), 1, tolerance = 1e-10)
})

test_that("mv_mixture_fit recovers well-separated clusters", {
  set.seed(123)
  # Very well-separated groups
  x_list <- c(
    lapply(1:15, function(i) matrix(rnorm(6, mean = 10, sd = 0.1), 2, 3)),
    lapply(1:15, function(i) matrix(rnorm(6, mean = -10, sd = 0.1), 2, 3))
  )
  fit <- mv_mixture_fit(x_list, g = 2, max_iter = 50, verbose = FALSE)
  # Check cluster purity: each group should be homogeneous
  group1 <- fit$cluster[1:15]
  group2 <- fit$cluster[16:30]
  expect_true(length(unique(group1)) == 1)
  expect_true(length(unique(group2)) == 1)
  expect_true(group1[1] != group2[1])
})

test_that("mv_mixture_fit log-likelihood is non-decreasing", {
  set.seed(7)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_mixture_fit(x_list, g = 2, max_iter = 30, verbose = FALSE)
  ll <- fit$logLik
  if (length(ll) > 1) {
    diffs <- diff(ll)
    # Allow small numerical noise
    expect_true(all(diffs >= -1e-4))
  }
})

test_that("mv_mixture_fit works with kmeans init", {
  set.seed(99)
  x_list <- lapply(1:12, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_mixture_fit(x_list, g = 2, max_iter = 10, init = "kmeans", verbose = FALSE)
  expect_length(fit$cluster, 12)
})

test_that("mv_mixture_fit works with emrefine init", {
  set.seed(99)
  x_list <- lapply(1:12, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_mixture_fit(x_list, g = 2, max_iter = 10, init = "emrefine", verbose = FALSE)
  expect_length(fit$cluster, 12)
})

test_that("mv_mixture_fit responsibilities sum to 1 per row", {
  set.seed(55)
  x_list <- lapply(1:12, function(i) matrix(rnorm(6), 2, 3))
  fit <- mv_mixture_fit(x_list, g = 2, max_iter = 20, verbose = FALSE)
  row_sums <- rowSums(fit$z)
  expect_equal(row_sums, rep(1, 12), tolerance = 1e-8)
})
