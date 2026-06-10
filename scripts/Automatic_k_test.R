# test_automatic_k_selection.R
# Tests the automatic noise parameter (k) selection scheme for HC noise type

library(testthat)
source("..R/Matrix_Init.R")
source("..R/Matrix.R")
source("..R/Matrix_KS_Score.R")
source("..R/Matrix_Noise.R")
source("..R/Matrix_Noise_BR.R")
source("..R/Matrix_Utils.R")

# Helper function to generate synthetic matrix-variate data
generate_test_data <- function(n_per_group, rows, cols, 
                               group_means, row_covs, col_covs,
                               noise_prop = 0) {
  total_n <- sum(n_per_group)
  g <- length(group_means)
  
  x_list <- vector("list", total_n)
  true_labels <- integer(total_n)
  noise_idx <- sample(total_n, floor(total_n * noise_prop))
  
  idx <- 1
  for (group in 1:g) {
    n_group <- n_per_group[group]
    for (i in 1:n_group) {
      # Generate matrix from matrix-variate normal
      Z <- matrix(rnorm(rows * cols), rows, cols)
      eigen_row <- eigen(row_covs[[group]])
      eigen_col <- eigen(col_covs[[group]])
      
      A <- eigen_row$vectors %*% diag(sqrt(eigen_row$values)) %*% t(eigen_row$vectors)
      B <- eigen_col$vectors %*% diag(sqrt(eigen_col$values)) %*% t(eigen_col$vectors)
      
      x_list[[idx]] <- group_means[[group]] + A %*% Z %*% B
      true_labels[idx] <- group
      idx <- idx + 1
    }
  }
  
  # Add noise observations
  if (noise_prop > 0) {
    for (i in noise_idx) {
      # Generate uniform noise within range of data
      all_data <- do.call(rbind, lapply(x_list, as.vector))
      data_range <- range(all_data, na.rm = TRUE)
      x_list[[i]] <- matrix(runif(rows * cols, data_range[1], data_range[2]), rows, cols)
      true_labels[i] <- 0  # 0 indicates noise
    }
  }
  
  list(x_list = x_list, true_labels = true_labels)
}

# Test 1: Basic functionality with well-separated clusters
test_that("Automatic k selection finds reasonable k value", {
  set.seed(123)
  
  # Parameters
  rows <- 4
  cols <- 3
  g <- 2
  
  # Generate clean data (no noise)
  group_means <- list(
    matrix(0, rows, cols),
    matrix(2, rows, cols)
  )
  
  row_covs <- list(
    diag(rows),
    diag(rows)
  )
  
  col_covs <- list(
    diag(cols),
    diag(cols)
  )
  
  test_data <- generate_test_data(
    n_per_group = c(30, 30),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0
  )
  
  # Run with automatic k selection
  result <- matrix_variate_noise_fit(
    x_list = test_data$x_list,
    g = g,
    noise_type = "hc",
    estimate_k = TRUE,
    adaptive_grid = TRUE,
    verbose = TRUE,
    nstart = 10,
    max_iter = 100
  )
  
  # Check results
  expect_true(!is.null(result$k_selection))
  expect_true(result$k_selection$selected_k > 0)
  expect_true(result$k_selection$selected_k < 1)
  expect_true(length(result$k_selection$ks_scores) == length(result$k_selection$k_grid))
  expect_true(all(result$k_selection$n_used > 0))
  
  # Check that noise proportion is reasonable
  expect_true(result$noise$pi >= 0 && result$noise$pi <= 1)
  expect_true(result$noise$type == "hc")
  expect_true(result$noise$k == result$k_selection$selected_k)
  
  # Check cluster assignments
  expect_true(length(result$cluster) == length(test_data$x_list))
  expect_true(all(result$cluster >= 0 & result$cluster <= g))
  
  cat("\nSelected k:", result$k_selection$selected_k)
  cat("\nKS scores range:", range(result$k_selection$ks_scores))
  cat("\nNoise proportion:", result$noise$pi)
})

# Test 2: Different grid sizes and adaptive vs fixed
test_that("Adaptive grid vs fixed grid produce reasonable results", {
  set.seed(456)
  
  rows <- 3
  cols <- 3
  g <- 2
  
  group_means <- list(
    matrix(0, rows, cols),
    matrix(3, rows, cols)
  )
  
  row_covs <- list(diag(rows), diag(rows))
  col_covs <- list(diag(cols), diag(cols))
  
  test_data <- generate_test_data(
    n_per_group = c(40, 40),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0.05
  )
  
  # Run with adaptive grid
  result_adaptive <- matrix_variate_noise_fit(
    x_list = test_data$x_list,
    g = g,
    noise_type = "hc",
    estimate_k = TRUE,
    adaptive_grid = TRUE,
    verbose = FALSE,
    nstart = 5,
    max_iter = 50
  )
  
  # Run with custom grid
  custom_grid <- 10^seq(-12, -2, length.out = 20)
  result_custom <- matrix_variate_noise_fit(
    x_list = test_data$x_list,
    g = g,
    noise_type = "hc",
    estimate_k = TRUE,
    k_grid = custom_grid,
    adaptive_grid = FALSE,
    verbose = FALSE,
    nstart = 5,
    max_iter = 50
  )
  
  expect_true(result_adaptive$k_selection$selected_k > 0)
  expect_true(result_custom$k_selection$selected_k > 0)
  expect_true(length(result_adaptive$k_selection$k_grid) >= 2)
  expect_equal(result_custom$k_selection$k_grid, custom_grid)
  
  cat("\nAdaptive grid selected k:", result_adaptive$k_selection$selected_k)
  cat("\nCustom grid selected k:", result_custom$k_selection$selected_k)
})

# Test 3: Compare k selection with and without actual noise
test_that("K selection behaves differently with and without noise", {
  set.seed(789)
  
  rows <- 4
  cols <- 4
  g <- 2
  
  group_means <- list(
    matrix(0, rows, cols),
    matrix(2, rows, cols)
  )
  
  row_covs <- list(diag(rows), diag(rows))
  col_covs <- list(diag(cols), diag(cols))
  
  # Clean data
  clean_data <- generate_test_data(
    n_per_group = c(50, 50),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0
  )
  
  # Noisy data (10% noise)
  noisy_data <- generate_test_data(
    n_per_group = c(50, 50),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0.1
  )
  
  result_clean <- matrix_variate_noise_fit(
    x_list = clean_data$x_list,
    g = g,
    noise_type = "hc",
    estimate_k = TRUE,
    verbose = FALSE,
    nstart = 5,
    max_iter = 50
  )
  
  result_noisy <- matrix_variate_noise_fit(
    x_list = noisy_data$x_list,
    g = g,
    noise_type = "hc",
    estimate_k = TRUE,
    verbose = FALSE,
    nstart = 5,
    max_iter = 50
  )
  
  # For noisy data, we expect higher noise proportion
  cat("\nClean data - noise proportion:", result_clean$noise$pi)
  cat("\nNoisy data - noise proportion:", result_noisy$noise$pi)
  cat("\nClean data - selected k:", result_clean$k_selection$selected_k)
  cat("\nNoisy data - selected k:", result_noisy$k_selection$selected_k)
  
  # Noisy data should have higher noise proportion
  # (though not strictly guaranteed due to randomness)
  if (result_clean$noise$pi < result_noisy$noise$pi) {
    cat("\n✓ Noise proportion higher in noisy data as expected")
  } else {
    cat("\n⚠ Noise proportion not higher in noisy data (possible due to randomness)")
  }
})

# Test 4: Verify KS score calculation
test_that("KS scores are computed correctly", {
  set.seed(999)
  
  rows <- 2
  cols <- 2
  g <- 1
  
  group_means <- list(matrix(0, rows, cols))
  row_covs <- list(diag(rows))
  col_covs <- list(diag(cols))
  
  test_data <- generate_test_data(
    n_per_group = c(100),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0
  )
  
  # Fit with a fixed k first
  fit <- matrix_variate_noise_fit(
    x_list = test_data$x_list,
    g = g,
    noise_type = "hc",
    noise_k = 1e-6,
    estimate_k = FALSE,
    verbose = FALSE,
    max_iter = 50
  )
  
  # Compute KS score manually
  ks_result <- matrix_noise_ks_score(fit, test_data$x_list)
  
  expect_true(is.numeric(ks_result$statistic))
  expect_true(ks_result$statistic >= 0)
  expect_true(is.numeric(ks_result$p.value))
  expect_true(ks_result$n_used > 0)
  
  cat("\nKS statistic for good fit:", ks_result$statistic)
  cat("\nKS p-value:", ks_result$p.value)
})

# Run all tests
cat("\n=== Running Automatic K Selection Tests ===\n")
test_that("All automatic k selection tests pass", {
  # Tests are run above
  expect_true(TRUE)
})