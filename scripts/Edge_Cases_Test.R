# test_edge_cases.R
# Standalone test for edge cases with CSV export

source("../R/Matrix_Init.R")
source("../R/Matrix.R")
source("../R/Matrix_KS_Score.R")
source("../R/Matrix_Noise.R")
source("../R/Matrix_Noise_BR.R")
source("../R/Matrix_Utils.R")

# Helper function for generating test data
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
  
  if (noise_prop > 0) {
    for (i in noise_idx) {
      all_data <- do.call(rbind, lapply(x_list, as.vector))
      data_range <- range(all_data, na.rm = TRUE)
      x_list[[i]] <- matrix(runif(rows * cols, data_range[1], data_range[2]), rows, cols)
      true_labels[i] <- 0
    }
  }
  
  list(x_list = x_list, true_labels = true_labels)
}

# Save edge case results to CSV - FIXED VERSION
save_edge_results <- function(results_list, filename) {
  if (length(results_list) == 0) {
    cat("No results to save\n")
    return(NULL)
  }
  
  # Filter out NULL or invalid entries
  valid_results <- list()
  for (i in seq_along(results_list)) {
    if (!is.null(results_list[[i]]) && length(results_list[[i]]) > 0) {
      # Check if required fields exist
      if (!is.null(results_list[[i]]$test_case)) {
        valid_results[[length(valid_results) + 1]] <- results_list[[i]]
      }
    }
  }
  
  if (length(valid_results) == 0) {
    cat("No valid results to save\n")
    return(NULL)
  }
  
  # Build data frame row by row to avoid errors
  df_rows <- list()
  for (i in seq_along(valid_results)) {
    x <- valid_results[[i]]
    row <- data.frame(
      test_case = ifelse(is.null(x$test_case), "unknown", x$test_case),
      success = ifelse(is.null(x$success), FALSE, x$success),
      noise_proportion = ifelse(is.null(x$noise_proportion), 0, as.numeric(x$noise_proportion)),
      iterations = ifelse(is.null(x$iterations), 0, as.integer(x$iterations)),
      converged = ifelse(is.null(x$converged), FALSE, x$converged),
      selected_k = ifelse(is.null(x$selected_k), 0, as.numeric(x$selected_k)),
      message = ifelse(is.null(x$message), "", as.character(x$message)),
      stringsAsFactors = FALSE
    )
    df_rows[[i]] <- row
  }
  
  # Combine all rows
  if (length(df_rows) > 0) {
    final_df <- do.call(rbind, df_rows)
    write.csv(final_df, filename, row.names = FALSE)
    cat("Edge case results saved to:", filename, "\n")
    cat("  Saved", nrow(final_df), "results\n")
    return(final_df)
  } else {
    cat("No data frames created\n")
    return(NULL)
  }
}

cat("\n========================================\n")
cat("TEST 3: Edge Cases\n")
cat("========================================\n")

# Store results
edge_results <- list()

# Edge Case 1: Extremely small matrices
cat("\n--- Test 3.1: Extremely small matrices ---\n")
set.seed(1)

small_dims <- list(c(1,1), c(2,1), c(1,2))

for (dims in small_dims) {
  rows <- dims[1]
  cols <- dims[2]
  g <- 2
  
  cat("\nTesting", rows, "x", cols, "matrices:\n")
  
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
    cat("  ✗ Error:", e$message, "\n")
    return(NULL)
  })
  
  test_name <- paste0("small_", rows, "x", cols)
  if (!is.null(result) && !is.null(result$noise)) {
    cat("  ✓ Converged:", result$converged, "\n")
    cat("  ✓ Iterations:", result$iterations, "\n")
    cat("  ✓ Noise proportion:", result$noise$pi, "\n")
    
    edge_results[[test_name]] <- list(
      test_case = test_name,
      success = TRUE,
      noise_proportion = result$noise$pi,
      iterations = result$iterations,
      converged = result$converged,
      selected_k = if (!is.null(result$k_selection)) result$k_selection$selected_k else NA,
      message = "Success"
    )
  } else {
    cat("  ✗ Test failed\n")
    edge_results[[test_name]] <- list(
      test_case = test_name,
      success = FALSE,
      noise_proportion = NA,
      iterations = NA,
      converged = FALSE,
      selected_k = NA,
      message = "Failed to converge or NULL result"
    )
  }
}

# Edge Case 2: Single cluster with noise
cat("\n--- Test 3.2: Single cluster (g=1) with noise ---\n")
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

result <- tryCatch({
  matrix_variate_noise_fit(
    x_list = x_list,
    g = 1,
    noise_type = "hc",
    max_iter = 150,
    nstart = 10,
    estimate_k = TRUE,
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(result) && !is.null(result$noise)) {
  noise_detected <- sum(result$cluster == 0, na.rm = TRUE)
  cat("✓ Noise points detected:", noise_detected, "out of", n_noise, "\n")
  cat("✓ Noise proportion:", result$noise$pi, "\n")
  
  edge_results[["single_cluster"]] <- list(
    test_case = "single_cluster",
    success = TRUE,
    noise_proportion = result$noise$pi,
    iterations = result$iterations,
    converged = result$converged,
    selected_k = if (!is.null(result$k_selection)) result$k_selection$selected_k else NA,
    message = paste("Detected", noise_detected, "/", n_noise, "noise points")
  )
} else {
  cat("✗ Test failed\n")
  edge_results[["single_cluster"]] <- list(
    test_case = "single_cluster",
    success = FALSE,
    noise_proportion = NA,
    iterations = NA,
    converged = FALSE,
    selected_k = NA,
    message = "Test failed"
  )
}

# Edge Case 3: All points are noise
cat("\n--- Test 3.3: All points are noise ---\n")
set.seed(3)

rows <- 4
cols <- 4
n_total <- 50

# Generate pure noise
x_list <- vector("list", n_total)
for (i in 1:n_total) {
  x_list[[i]] <- matrix(runif(rows * cols, -10, 10), rows, cols)
}

result <- tryCatch({
  matrix_variate_noise_fit(
    x_list = x_list,
    g = 2,
    noise_type = "hc",
    max_iter = 100,
    nstart = 5,
    estimate_k = TRUE,
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(result) && !is.null(result$noise)) {
  cat("✓ Noise proportion:", result$noise$pi, "\n")
  cat("✓ Points assigned to noise:", sum(result$cluster == 0, na.rm = TRUE), "out of", n_total, "\n")
  
  edge_results[["pure_noise"]] <- list(
    test_case = "pure_noise",
    success = TRUE,
    noise_proportion = result$noise$pi,
    iterations = result$iterations,
    converged = result$converged,
    selected_k = if (!is.null(result$k_selection)) result$k_selection$selected_k else NA,
    message = paste("Assigned", sum(result$cluster == 0, na.rm = TRUE), "/", n_total, "to noise")
  )
} else {
  cat("✗ Test failed\n")
  edge_results[["pure_noise"]] <- list(
    test_case = "pure_noise",
    success = FALSE,
    noise_proportion = NA,
    iterations = NA,
    converged = FALSE,
    selected_k = NA,
    message = "Test failed"
  )
}

# Edge Case 4: High dimensional data
cat("\n--- Test 3.4: High dimensional matrices (10x10) ---\n")
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

result <- tryCatch({
  matrix_variate_noise_fit(
    x_list = x_list,
    g = g,
    noise_type = "hc",
    max_iter = 200,
    nstart = 5,
    estimate_k = TRUE,
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(result) && !is.null(result$noise)) {
  cat("✓ Converged:", result$converged, "\n")
  cat("✓ Iterations:", result$iterations, "\n")
  cat("✓ Noise proportion:", result$noise$pi, "\n")
  
  edge_results[["high_dim"]] <- list(
    test_case = "high_dim",
    success = result$converged,
    noise_proportion = result$noise$pi,
    iterations = result$iterations,
    converged = result$converged,
    selected_k = if (!is.null(result$k_selection)) result$k_selection$selected_k else NA,
    message = if(result$converged) "Converged" else "Did not converge"
  )
} else {
  cat("✗ Test failed\n")
  edge_results[["high_dim"]] <- list(
    test_case = "high_dim",
    success = FALSE,
    noise_proportion = NA,
    iterations = NA,
    converged = FALSE,
    selected_k = NA,
    message = "Test failed"
  )
}

# Edge Case 5: Extremely unbalanced groups
cat("\n--- Test 3.5: Extremely unbalanced group sizes ---\n")
set.seed(5)

rows <- 3
cols <- 3
g <- 3
n_per_group <- c(50, 5, 5)

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

result <- tryCatch({
  matrix_variate_noise_fit(
    x_list = test_data$x_list,
    g = g,
    noise_type = "hc",
    max_iter = 150,
    nstart = 20,
    estimate_k = TRUE,
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(result) && !is.null(result$noise)) {
  cluster_sizes <- table(result$cluster[result$cluster > 0])
  cat("✓ Detected cluster sizes:", paste(cluster_sizes, collapse = ", "), "\n")
  cat("✓ Number of clusters found:", length(cluster_sizes), "\n")
  
  edge_results[["unbalanced"]] <- list(
    test_case = "unbalanced",
    success = TRUE,
    noise_proportion = result$noise$pi,
    iterations = result$iterations,
    converged = result$converged,
    selected_k = if (!is.null(result$k_selection)) result$k_selection$selected_k else NA,
    message = paste("Found", length(cluster_sizes), "clusters")
  )
} else {
  cat("✗ Test failed\n")
  edge_results[["unbalanced"]] <- list(
    test_case = "unbalanced",
    success = FALSE,
    noise_proportion = NA,
    iterations = NA,
    converged = FALSE,
    selected_k = NA,
    message = "Test failed"
  )
}

# Edge Case 6: BR noise type
cat("\n--- Test 3.6: BR noise type ---\n")
set.seed(7)

rows <- 3
cols <- 3
g <- 2

# Generate normal data
x_list <- vector("list", 40)
for (i in 1:20) {
  x_list[[i]] <- matrix(rnorm(rows * cols, 0, 1), rows, cols)
}
for (i in 21:40) {
  x_list[[i]] <- matrix(rnorm(rows * cols, 5, 1), rows, cols)
}

result_br <- tryCatch({
  matrix_variate_noise_fit(
    x_list = x_list,
    g = g,
    noise_type = "br",
    max_iter = 150,
    nstart = 10,
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(result_br) && !is.null(result_br$noise)) {
  cat("✓ BR noise - converged:", result_br$converged, "\n")
  cat("✓ BR noise - iterations:", result_br$iterations, "\n")
  cat("✓ BR noise - noise proportion:", result_br$noise$pi, "\n")
  
  edge_results[["br_noise"]] <- list(
    test_case = "br_noise",
    success = result_br$converged,
    noise_proportion = result_br$noise$pi,
    iterations = result_br$iterations,
    converged = result_br$converged,
    selected_k = NA,  # BR doesn't use k
    message = "BR noise type test"
  )
} else {
  cat("✗ Test failed\n")
  edge_results[["br_noise"]] <- list(
    test_case = "br_noise",
    success = FALSE,
    noise_proportion = NA,
    iterations = NA,
    converged = FALSE,
    selected_k = NA,
    message = "Test failed"
  )
}

# Save all edge case results
cat("\n--- Saving Edge Case Results ---\n")
if (length(edge_results) > 0) {
  edge_df <- save_edge_results(edge_results, "edge_case_results.csv")
  
  if (!is.null(edge_df) && nrow(edge_df) > 0) {
    # Print summary
    cat("\n--- Edge Case Summary ---\n")
    print(edge_df)
  } else {
    cat("No valid results to display\n")
  }
} else {
  cat("No results to save\n")
}

cat("\n=== Edge Case Testing Complete ===\n")