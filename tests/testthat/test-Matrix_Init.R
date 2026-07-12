## ---- mv_mixture_kmeans_init ----

test_that("kmeans init returns correct structure", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_kmeans_init(x_list, g = 2)

  expect_type(params, "list")
  expect_named(params, c("pi", "M", "U", "V", "cluster"))
  expect_length(params$pi, 2)
  expect_length(params$M, 2)
  expect_length(params$U, 2)
  expect_length(params$V, 2)
  expect_length(params$cluster, 15)
})

test_that("kmeans init mixing proportions sum to 1", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_kmeans_init(x_list, g = 2)
  expect_equal(sum(params$pi), 1, tolerance = 1e-10)
})

test_that("kmeans init produces SPD covariance matrices", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_kmeans_init(x_list, g = 2)
  for (j in 1:2) {
    expect_silent(chol(params$U[[j]]))
    expect_silent(chol(params$V[[j]]))
  }
})

test_that("kmeans init mean matrices have correct dimensions", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_kmeans_init(x_list, g = 3)
  for (j in 1:3) {
    expect_equal(dim(params$M[[j]]), c(2, 3))
    expect_equal(dim(params$U[[j]]), c(2, 2))
    expect_equal(dim(params$V[[j]]), c(3, 3))
  }
})

## ---- mv_mixture_kmeans_init_kmeanspp_behavior ----

test_that("kmeans init returns correct structure", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_kmeans_init(x_list, g = 2)

  expect_type(params, "list")
  expect_named(params, c("pi", "M", "U", "V", "cluster"))
  expect_length(params$pi, 2)
  expect_length(params$M, 2)
  expect_length(params$U, 2)
  expect_length(params$V, 2)
  expect_length(params$cluster, 15)
})

test_that("kmeans init mixing proportions sum to 1", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_kmeans_init(x_list, g = 2)
  expect_equal(sum(params$pi), 1, tolerance = 1e-10)
})

test_that("kmeans init produces SPD covariance matrices", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_kmeans_init(x_list, g = 2)
  for (j in 1:2) {
    expect_silent(chol(params$U[[j]]))
    expect_silent(chol(params$V[[j]]))
  }
})

test_that("kmeans init mean matrices have correct dimensions", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_kmeans_init(x_list, g = 3)
  for (j in 1:3) {
    expect_equal(dim(params$M[[j]]), c(2, 3))
    expect_equal(dim(params$U[[j]]), c(2, 2))
    expect_equal(dim(params$V[[j]]), c(3, 3))
  }
})

test_that("kmeans init works via mixture_fit with init='kmeans'", {
  set.seed(123)
  x_list <- c(
    lapply(1:10, function(i) matrix(rnorm(6, mean = 5), 2, 3)),
    lapply(1:10, function(i) matrix(rnorm(6, mean = -5), 2, 3))
  )
  fit <- mv_mixture_fit(x_list, g = 2, max_iter = 20,
                                     init = "kmeans", verbose = FALSE)
  expect_length(fit$cluster, 20)
  expect_true(all(fit$cluster %in% 1:2))
  expect_equal(sum(fit$pi), 1, tolerance = 1e-10)
})

test_that("kmeans init selects spread-out centers", {
  set.seed(42)
  # Well-separated data — centers should not be from same cluster
  x_list <- c(
    lapply(1:15, function(i) matrix(rnorm(6, mean = 10, sd = 0.1), 2, 3)),
    lapply(1:15, function(i) matrix(rnorm(6, mean = -10, sd = 0.1), 2, 3))
  )
  params <- mv_mixture_kmeans_init(x_list, g = 2)
  # Both clusters should be represented
  group1 <- params$cluster[1:15]
  group2 <- params$cluster[16:30]
  expect_true(length(unique(group1)) == 1)
  expect_true(length(unique(group2)) == 1)
  expect_true(group1[1] != group2[1])
})

## ---- mv_mixture_emrefine_init ----

test_that("emrefine init returns correct structure", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_emrefine_init(x_list, g = 2)

  expect_named(params, c("pi", "M", "U", "V", "cluster"))
  expect_length(params$pi, 2)
  expect_equal(sum(params$pi), 1, tolerance = 1e-10)
  expect_length(params$cluster, 15)
})

test_that("emrefine init produces valid covariance matrices", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  params <- mv_mixture_emrefine_init(x_list, g = 2, max_iter = 3)
  for (j in 1:2) {
    expect_silent(chol(params$U[[j]]))
    expect_silent(chol(params$V[[j]]))
  }
})

test_that("emrefine init respects max_iter", {
  set.seed(42)
  x_list <- lapply(1:15, function(i) matrix(rnorm(6), 2, 3))
  # Should not error with just 1 iteration
  params <- mv_mixture_emrefine_init(x_list, g = 2, max_iter = 1)
  expect_type(params, "list")
})
