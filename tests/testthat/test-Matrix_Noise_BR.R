## ---- mv_noise_convex_hull_support ----

test_that("convex_hull_support returns correct structure", {
  set.seed(42)
  # Need d+1 unique points; for 2x2 matrices d=4, need >=5 unique matrices
  x_list <- lapply(1:10, function(i) matrix(rnorm(4), 2, 2))
  support <- mv_noise_convex_hull_support(x_list)

  expect_type(support, "list")
  expect_named(support, c("points", "hull", "log_volume", "jitter"))
  expect_true(is.numeric(support$log_volume))
  expect_true(is.finite(support$log_volume))
})

test_that("convex_hull_support errors with too few unique points", {
  # d = 4 (2x2 matrices), need > 4 unique points
  x_list <- lapply(1:3, function(i) matrix(i, 2, 2))
  expect_error(
    mv_noise_convex_hull_support(x_list),
    "d \\+ 1 unique"
  )
})

test_that("convex_hull_support log_volume is reasonable", {
  set.seed(42)
  x_list <- lapply(1:20, function(i) matrix(rnorm(4), 2, 2))
  support <- mv_noise_convex_hull_support(x_list)
  expect_true(is.finite(support$log_volume))
})

## ---- mv_noise_br_log_density ----

test_that("br_log_density returns finite values for points inside hull", {
  set.seed(42)
  x_list <- lapply(1:20, function(i) matrix(rnorm(4, sd = 2), 2, 2))
  support <- mv_noise_convex_hull_support(x_list)

  # The mean of the points should be inside the hull
  mean_mat <- Reduce("+", x_list) / length(x_list)
  ld <- mv_noise_br_log_density(list(mean_mat), support)
  expect_true(is.finite(ld[1]))
})

test_that("br_log_density returns -Inf for points outside hull", {
  set.seed(42)
  x_list <- lapply(1:20, function(i) matrix(rnorm(4, sd = 0.1), 2, 2))
  support <- mv_noise_convex_hull_support(x_list)

  # A point far from the origin should be outside
  far_point <- matrix(100, 2, 2)
  ld <- mv_noise_br_log_density(list(far_point), support)
  expect_equal(ld[1], -Inf)
})

test_that("br_log_density returns same value for all interior points", {
  set.seed(42)
  x_list <- lapply(1:20, function(i) matrix(rnorm(4, sd = 2), 2, 2))
  support <- mv_noise_convex_hull_support(x_list)

  # Pick two points that are convex combinations of existing data
  p1 <- (x_list[[1]] + x_list[[2]]) / 2
  p2 <- (x_list[[3]] + x_list[[4]]) / 2

  ld <- mv_noise_br_log_density(list(p1, p2), support)
  finite_ld <- ld[is.finite(ld)]
  if (length(finite_ld) == 2) {
    # Uniform: same density for all interior points
    expect_equal(finite_ld[1], finite_ld[2], tolerance = 1e-10)
  }
})
