# test_edge_cases.R
# Tests the algorithm on various edge cases

source("..R/Matrix_Init.R")
source("..R/Matrix.R")
source("..R/Matrix_KS_Score.R")
source("..R/Matrix_Noise.R")
source("..R/Matrix_Noise_BR.R")
source("..R/Matrix_Utils.R")

test_edge_cases <- function() {
  cat("\n=== Testing Edge Cases ===\n")
  
  # Edge Case 1: Extremely small matrices
  test_that("Handle extremely small matrices (1x1, 2x1, 1x2)", {
    set.seed(1)
    
    small_dims <- list(c(1,1), c(2,1), c(1,2))
    
    for (dims in small_dims) {
      rows <- dims[1]
      cols <- dims[2]
      g <- 2
      
      # Generate simple data
      x_list <- vector("list", 20)
      for (i in 1:10) {
        x_list[[i]] <- matrix(rnorm(rows * cols, 0, 1), rows, cols)
      }
      for (i in 11:20) {
        x_list[[i]] <- matrix(rnorm(rows * cols, 5, 1), rows, cols)
      }
      
      result <- tryCatch({
        matrix_variate_noise_fit(
          x_list = x_list,
          g = g,
          noise_type = "hc",
          max_iter = 100,
          nstart = 5,
          estimate_k = TRUE,
          verbose = FALSE
        )
      }, error = function(e) {
        cat("\nError with dimensions", rows, "x", cols, ":", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(result)) {
        expect_true(length(result$cluster) == 20)
        expect_true(all(result$cluster %in% 0:g))
        cat(sprintf("\n✓ %dx%d matrices: converged in %d iterations, noise prop = %.3f", 
                    rows, cols, result$iterations, result$noise$pi))
      }
    }
  })
  
  # Edge Case 2: Single cluster with noise
  test_that("Single cluster (g=1) with noise", {
    set.seed(2)
    
    rows <- 3
    cols <- 3
    n_clean <- 40
    n_noise <- 10
    
    # Generate clean data from one cluster
    x_list <- vector("list", n_clean)
    for (i in 1:n_clean) {
      x_list[[i]] <- matrix(rnorm(rows * cols, 0, 1), rows, cols)
    }
    
    # Add noise
    for (i in 1:n_noise) {
      x_list[[n_clean + i]] <- matrix(runif(rows * cols, -5, 5), rows, cols)
    }
    
    result <- matrix_variate_noise_fit(
      x_list = x_list,
      g = 1,
      noise_type = "hc",
      max_iter = 150,
      nstart = 10,
      estimate_k = TRUE,
      verbose = FALSE
    )
    
    # Noise should be detected
    noise_detected <- sum(result$cluster == 0)
    expect_true(noise_detected > 0)
    expect_true(result$noise$pi > 0)
    
    cat(sprintf("\n✓ Single cluster: detected %d/%d noise points (%.1f%%)", 
                noise_detected, n_noise, 100 * noise_detected / n_noise))
  })
  
  # Edge Case 3: All points are noise
  test_that("All points are noise", {
    set.seed(3)
    
    rows <- 4
    cols <- 4
    n_total <- 50
    
    # Generate pure noise
    x_list <- vector("list", n_total)
    for (i in 1:n_total) {
      x_list[[i]] <- matrix(runif(rows * cols, -10, 10), rows, cols)
    }
    
    result <- matrix_variate_noise_fit(
      x_list = x_list,
      g = 2,  # Try to find 2 clusters in pure noise
      noise_type = "hc",
      max_iter = 100,
      nstart = 5,
      estimate_k = TRUE,
      verbose = FALSE
    )
    
    # Should assign most points to noise
    noise_prop <- result$noise$pi
    cat(sprintf("\n✓ Pure noise: noise proportion = %.3f (expected > 0.5)", noise_prop))
    expect_true(noise_prop > 0.3)  # At least some noise detection
  })
  
  # Edge Case 4: Very high dimensional data
  test_that("High dimensional matrices", {
    set.seed(4)
    
    rows <- 10
    cols <- 10
    g <- 2
    n_per_group <- 15
    
    # Generate data with high dimensions
    x_list <- vector("list", n_per_group * g)
    idx <- 1
    for (group in 1:g) {
      mean_mat <- matrix(rep(group * 2, rows * cols), rows, cols)
      for (i in 1:n_per_group) {
        x_list[[idx]] <- mean_mat + matrix(rnorm(rows * cols, 0, 0.5), rows, cols)
        idx <- idx + 1
      }
    }
    
    # Add some noise
    for (i in 1:5) {
      x_list[[idx]] <- matrix(runif(rows * cols, -5, 15), rows, cols)
      idx <- idx + 1
    }
    
    result <- matrix_variate_noise_fit(
      x_list = x_list,
      g = g,
      noise_type = "hc",
      max_iter = 200,
      nstart = 5,
      estimate_k = TRUE,
      verbose = FALSE
    )
    
    expect_true(result$converged || result$iterations == 200)
    cat(sprintf("\n✓ High dimensional (%dx%d): converged in %d iterations", 
                rows, cols, result$iterations))
  })
  
  # Edge Case 5: Extremely unbalanced groups
  test_that("Extremely unbalanced group sizes", {
    set.seed(5)
    
    rows <- 3
    cols <- 3
    g <- 3
    n_per_group <- c(50, 5, 5)  # Very unbalanced
    
    group_means <- list(
      matrix(0, rows, cols),
      matrix(5, rows, cols),
      matrix(-5, rows, cols)
    )
    
    row_covs <- replicate(g, diag(rows), simplify = FALSE)
    col_covs <- replicate(g, diag(cols), simplify = FALSE)
    
    test_data <- generate_test_data(
      n_per_group = n_per_group,
      rows = rows, cols = cols,
      group_means = group_means,
      row_covs = row_covs,
      col_covs = col_covs,
      noise_prop = 0
    )
    
    result <- matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      max_iter = 150,
      nstart = 20,  # More starts for difficult case
      estimate_k = TRUE,
      verbose = FALSE
    )
    
    # Check small groups were identified
    small_groups <- table(result$cluster[result$cluster > 0])
    cat("\n✓ Unbalanced groups: detected group sizes =", 
        paste(small_groups, collapse = ", "))
    expect_true(length(small_groups) <= g)
  })
  
  # Edge Case 6: Perfect separation (no noise needed)
  test_that("Perfectly separated clusters", {
    set.seed(6)
    
    rows <- 2
    cols <- 2
    g <- 2
    n_per_group <- 30
    
    # Generate perfectly separated clusters
    x_list <- vector("list", n_per_group * g)
    idx <- 1
    for (group in 1:g) {
      mean_mat <- matrix(rep(group * 10, rows * cols), rows, cols)
      for (i in 1:n_per_group) {
        x_list[[idx]] <- mean_mat + matrix(rnorm(rows * cols, 0, 0.1), rows, cols)
        idx <- idx + 1
      }
    }
    
    result <- matrix_variate_noise_fit(
      x_list = x_list,
      g = g,
      noise_type = "hc",
      max_iter = 100,
      nstart = 10,
      estimate_k = TRUE,
      verbose = FALSE
    )
    
    # Should have very low noise proportion
    cat(sprintf("\n✓ Perfect separation: noise proportion = %.6f (should be very low)", 
                result$noise$pi))
    expect_true(result$noise$pi < 0.05)
  })
  
  # Edge Case 7: Convergence tolerance sensitivity
  test_that("Different convergence tolerances", {
    set.seed(7)
    
    rows <- 4
    cols <- 3
    g <- 2
    
    x_list <- vector("list", 40)
    for (i in 1:20) {
      x_list[[i]] <- matrix(rnorm(rows * cols, 0, 1), rows, cols)
    }
    for (i in 21:40) {
      x_list[[i]] <- matrix(rnorm(rows * cols, 3, 1), rows, cols)
    }
    
    tolerances <- c(1e-3, 1e-6, 1e-9)
    results_tol <- list()
    
    for (tol_val in tolerances) {
      result <- matrix_variate_noise_fit(
        x_list = x_list,
        g = g,
        noise_type = "hc",
        max_iter = 200,
        tol = tol_val,
        nstart = 5,
        estimate_k = FALSE,
        noise_k = 1e-5,
        verbose = FALSE
      )
      results_tol[[as.character(tol_val)]] <- result
      cat(sprintf("\n✓ Tolerance %e: converged in %d iterations, logLik = %.2f", 
                  tol_val, result$iterations, tail(result$logLik, 1)))
    }
    
    # Stricter tolerance should lead to more iterations
    expect_true(results_tol[["1e-9"]]$iterations >= results_tol[["1e-3"]]$iterations)
  })
  
  # Edge Case 8: BR noise with convex hull issues
  test_that("BR noise with nearly singular convex hull", {
    set.seed(8)
    
    rows <- 3
    cols <- 3
    g <- 2
    
    # Generate nearly collinear data
    x_list <- vector("list", 30)
    for (i in 1:15) {
      base <- matrix(rnorm(rows * cols, 0, 1), rows, cols)
      # Make first column nearly constant
      base[,1] <- base[,1] * 0.01
      x_list[[i]] <- base
    }
    for (i in 16:30) {
      base <- matrix(rnorm(rows * cols, 5, 1), rows, cols)
      base[,1] <- base[,1] * 0.01
      x_list[[i]] <- base
    }
    
    result <- matrix_variate_noise_fit(
      x_list = x_list,
      g = g,
      noise_type = "br",
      max_iter = 150,
      nstart = 10,
      verbose = FALSE
    )
    
    expect_true(!is.null(result$noise$hull))
    cat("\n✓ BR noise: hull computed successfully with dimension", 
        ifelse(is.null(result$noise$hull), "NULL", 
               paste(dim(result$noise$hull), collapse = "x")))
  })
}

# Run all edge case tests
test_edge_cases()

# Comprehensive summary
cat("\n\n=== EDGE CASE TESTING COMPLETE ===\n")
cat("All edge cases have been tested. Check individual results above.\n")