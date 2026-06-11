# extensive_tests_complete.R
# Complete extensive testing suite for matrix-variate noise mixture clustering

source("../R/Matrix_Init.R")
source("../R/Matrix.R")
source("../R/Matrix_KS_Score.R")
source("../R/Matrix_Noise.R")
source("../R/Matrix_Noise_BR.R")
source("../R/Matrix_Utils.R")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Enhanced data generator with more control
generate_complex_data <- function(
    n_per_group, rows, cols, 
    group_means, row_covs, col_covs,
    noise_prop = 0,
    noise_type = "uniform",
    overlap = 0,
    outliers = 0
) {
  g <- length(group_means)
  total_clean <- sum(n_per_group)
  n_noise <- floor(total_clean * noise_prop)
  n_outliers <- floor(total_clean * outliers)
  total_n <- total_clean + n_noise + n_outliers
  
  x_list <- vector("list", total_n)
  true_labels <- integer(total_n)
  
  # Generate clean clusters
  idx <- 1
  for (group in 1:g) {
    n_group <- n_per_group[group]
    for (i in 1:n_group) {
      if (overlap > 0 && group < g) {
        mix_prob <- runif(1)
        if (mix_prob < overlap) {
          actual_mean <- (group_means[[group]] + group_means[[group + 1]]) / 2
        } else {
          actual_mean <- group_means[[group]]
        }
      } else {
        actual_mean <- group_means[[group]]
      }
      
      Z <- matrix(rnorm(rows * cols), rows, cols)
      eigen_row <- eigen(row_covs[[group]])
      eigen_col <- eigen(col_covs[[group]])
      
      A <- eigen_row$vectors %*% diag(sqrt(eigen_row$values)) %*% t(eigen_row$vectors)
      B <- eigen_col$vectors %*% diag(sqrt(eigen_col$values)) %*% t(eigen_col$vectors)
      
      x_list[[idx]] <- actual_mean + A %*% Z %*% B
      true_labels[idx] <- group
      idx <- idx + 1
    }
  }
  
  # Add standard noise
  if (n_noise > 0) {
    all_data <- do.call(rbind, lapply(x_list[1:(total_clean - n_outliers)], as.vector))
    data_range <- range(all_data, na.rm = TRUE)
    data_sd <- sd(all_data, na.rm = TRUE)
    
    for (i in 1:n_noise) {
      if (noise_type == "uniform") {
        noise_mat <- matrix(runif(rows * cols, data_range[1], data_range[2]), rows, cols)
      } else if (noise_type == "gaussian") {
        noise_mat <- matrix(rnorm(rows * cols, mean(all_data), data_sd), rows, cols)
      } else if (noise_type == "cauchy") {
        noise_mat <- matrix(rcauchy(rows * cols, 0, data_sd), rows, cols)
      }
      x_list[[idx]] <- noise_mat
      true_labels[idx] <- 0
      idx <- idx + 1
    }
  }
  
  # Add extreme outliers
  if (n_outliers > 0) {
    all_data <- do.call(rbind, lapply(x_list[1:idx-1], as.vector))
    data_range <- range(all_data, na.rm = TRUE)
    extreme_range <- data_range + c(-10, 10) * sd(all_data)
    
    for (i in 1:n_outliers) {
      x_list[[idx]] <- matrix(runif(rows * cols, extreme_range[1], extreme_range[2]), rows, cols)
      true_labels[idx] <- -1
      idx <- idx + 1
    }
  }
  
  list(x_list = x_list, true_labels = true_labels)
}

# Generate scalable data
generate_scalable_data <- function(n_points, rows, cols, g, noise_prop = 0.1) {
  set.seed(42)
  
  group_means <- list()
  row_covs <- list()
  col_covs <- list()
  
  for (i in 1:g) {
    group_means[[i]] <- matrix(rnorm(rows * cols, i*2, 1), rows, cols)
    row_covs[[i]] <- diag(rows) * (1 + runif(1, 0, 1))
    col_covs[[i]] <- diag(cols) * (1 + runif(1, 0, 1))
  }
  
  n_per_group <- rep(floor(n_points / g), g)
  n_per_group[g] <- n_points - sum(n_per_group[1:(g-1)])
  
  x_list <- vector("list", n_points)
  true_labels <- integer(n_points)
  
  idx <- 1
  for (group in 1:g) {
    for (i in 1:n_per_group[group]) {
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
  
  n_noise <- floor(n_points * noise_prop)
  if (n_noise > 0) {
    all_data <- do.call(rbind, lapply(x_list, as.vector))
    data_range <- range(all_data)
    for (i in 1:n_noise) {
      noise_idx <- sample(n_points, 1)
      x_list[[noise_idx]] <- matrix(runif(rows * cols, data_range[1], data_range[2]), rows, cols)
      true_labels[noise_idx] <- 0
    }
  }
  
  list(x_list = x_list, true_labels = true_labels)
}

# Comprehensive evaluation function
evaluate_comprehensive <- function(result, true_labels, x_list, test_name) {
  noise_true <- true_labels == 0
  noise_pred <- result$cluster == 0
  
  noise_accuracy <- mean(noise_true == noise_pred)
  noise_precision <- if(sum(noise_pred) > 0) sum(noise_true & noise_pred) / sum(noise_pred) else NA
  noise_recall <- if(sum(noise_true) > 0) sum(noise_true & noise_pred) / sum(noise_true) else NA
  noise_f1 <- if(!is.na(noise_precision) && !is.na(noise_recall) && (noise_precision + noise_recall) > 0) {
    2 * noise_precision * noise_recall / (noise_precision + noise_recall)
  } else { NA }
  
  non_noise_idx <- which(true_labels > 0 & result$cluster > 0)
  if (length(non_noise_idx) > 0) {
    true_non_noise <- true_labels[non_noise_idx]
    pred_non_noise <- result$cluster[non_noise_idx]
    
    n <- length(non_noise_idx)
    agreements <- 0
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        same_true <- (true_non_noise[i] == true_non_noise[j])
        same_pred <- (pred_non_noise[i] == pred_non_noise[j])
        if (same_true == same_pred) agreements <- agreements + 1
      }
    }
    total_pairs <- n * (n - 1) / 2
    ari <- if(total_pairs > 0) agreements / total_pairs else 1
  } else {
    ari <- NA
  }
  
  selected_k <- if(!is.null(result$k_selection)) result$k_selection$selected_k else NA
  ks_stat <- if(!is.null(result$k_selection)) min(unlist(result$k_selection$ks_scores), na.rm = TRUE) else NA
  ks_pval <- if(!is.null(result$k_selection)) max(unlist(result$k_selection$ks_pvalues), na.rm = TRUE) else NA
  
  return(data.frame(
    test_name = test_name,
    selected_k = selected_k,
    noise_proportion_estimated = result$noise$pi,
    noise_proportion_true = mean(noise_true),
    noise_accuracy = noise_accuracy,
    noise_f1 = noise_f1,
    ari = ari,
    iterations = result$iterations,
    converged = result$converged,
    ks_pvalue = ks_pval,
    stringsAsFactors = FALSE
  ))
}

# ============================================================================
# TEST SUITE 1: Varying Noise Proportions
# ============================================================================

cat("\n")
cat(rep("=", 80), collapse = "")
cat("\nEXTENSIVE TEST SUITE: Complete Analysis\n")
cat(rep("=", 80), collapse = "")
cat("\n")

all_results <- list()

# Test 1: Varying noise proportions
cat("\n--- Test 1: Varying Noise Proportions (0% to 30%) ---\n")

rows <- 4
cols <- 3
g <- 2
group_means <- list(matrix(0, rows, cols), matrix(3, rows, cols))
row_covs <- list(diag(rows), diag(rows))
col_covs <- list(diag(cols), diag(cols))

noise_levels <- c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

for (noise_level in noise_levels) {
  cat(sprintf("  Noise: %.0f%%", noise_level * 100))
  
  test_data <- generate_complex_data(
    n_per_group = c(50, 50),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = noise_level,
    noise_type = "uniform"
  )
  
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 10,
      max_iter = 100
    )
  }, error = function(e) {
    cat(" ✗ Error\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("noise_", noise_level*100))
    all_results[[length(all_results)+1]] <- eval
    cat(sprintf(" ✓ Detected: %.1f%%, ARI: %.3f, k: %.2e\n", 
                eval$noise_proportion_estimated * 100, eval$ari, eval$selected_k))
  }
}

# Test 2: Different noise distributions
cat("\n--- Test 2: Different Noise Distributions ---\n")

noise_distributions <- c("uniform", "gaussian", "cauchy")

for (dist in noise_distributions) {
  cat(sprintf("  Distribution: %s", dist))
  
  test_data <- generate_complex_data(
    n_per_group = c(50, 50),
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0.15,
    noise_type = dist
  )
  
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 10,
      max_iter = 100
    )
  }, error = function(e) {
    cat(" ✗ Error\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("dist_", dist))
    all_results[[length(all_results)+1]] <- eval
    cat(sprintf(" ✓ F1: %.3f, ARI: %.3f, k: %.2e\n", 
                eval$noise_f1, eval$ari, eval$selected_k))
  }
}

# Test 3: Varying cluster separation
cat("\n--- Test 3: Varying Cluster Separation ---\n")

separations <- c(0.5, 1, 2, 3, 4, 5)

for (sep in separations) {
  cat(sprintf("  Separation: %.1f", sep))
  
  group_means_sep <- list(matrix(0, rows, cols), matrix(sep, rows, cols))
  
  test_data <- generate_complex_data(
    n_per_group = c(50, 50),
    rows = rows, cols = cols,
    group_means = group_means_sep,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0.10
  )
  
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 10,
      max_iter = 100
    )
  }, error = function(e) {
    cat(" ✗ Error\n")
    return(NULL)
  })
  
  if (!is.null(result)) {
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("sep_", sep))
    all_results[[length(all_results)+1]] <- eval
    cat(sprintf(" ✓ ARI: %.3f, Noise: %.1f%%, k: %.2e\n", 
                eval$ari, eval$noise_proportion_estimated * 100, eval$selected_k))
  }
}

# Test 4: Scalability with points
cat("\n--- Test 4: Scaling with Number of Points ---\n")

point_counts <- c(50, 100, 200, 500)
rows <- 3
cols <- 3
g <- 2

for (n_pts in point_counts) {
  cat(sprintf("  Points: %d", n_pts))
  
  test_data <- generate_scalable_data(n_pts, rows, cols, g, noise_prop = 0.1)
  
  start_time <- Sys.time()
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 5,
      max_iter = 100
    )
  }, error = function(e) {
    cat(" ✗ Error\n")
    return(NULL)
  })
  end_time <- Sys.time()
  
  if (!is.null(result)) {
    exec_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("n", n_pts))
    all_results[[length(all_results)+1]] <- eval
    cat(sprintf(" ✓ Time: %.2f sec, Iter: %d, ARI: %.3f\n", 
                exec_time, result$iterations, eval$ari))
  }
}

# Test 5: Scalability with dimensions
cat("\n--- Test 5: Scaling with Matrix Dimensions ---\n")

dimensions <- list(c(2,2), c(3,3), c(5,5), c(8,8), c(10,10))
n_points <- 100

for (dims in dimensions) {
  rows <- dims[1]
  cols <- dims[2]
  cat(sprintf("  Size: %dx%d", rows, cols))
  
  test_data <- generate_scalable_data(n_points, rows, cols, g = 2, noise_prop = 0.1)
  
  start_time <- Sys.time()
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = 2,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 5,
      max_iter = 100
    )
  }, error = function(e) {
    cat(" ✗ Error\n")
    return(NULL)
  })
  end_time <- Sys.time()
  
  if (!is.null(result)) {
    exec_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("dim", rows, "x", cols))
    all_results[[length(all_results)+1]] <- eval
    cat(sprintf(" ✓ Time: %.2f sec, k: %.2e, Conv: %s\n", 
                exec_time, eval$selected_k, result$converged))
  }
}

# Test 6: Varying number of clusters
cat("\n--- Test 6: Scaling with Number of Clusters ---\n")

cluster_counts <- c(2, 3, 4, 5)
rows <- 4
cols <- 4
n_points <- 200

for (g in cluster_counts) {
  cat(sprintf("  Clusters: %d", g))
  
  test_data <- generate_scalable_data(n_points, rows, cols, g, noise_prop = 0.1)
  
  start_time <- Sys.time()
  result <- tryCatch({
    matrix_variate_noise_fit(
      x_list = test_data$x_list,
      g = g,
      noise_type = "hc",
      estimate_k = TRUE,
      verbose = FALSE,
      nstart = 10,
      max_iter = 150
    )
  }, error = function(e) {
    cat(" ✗ Error\n")
    return(NULL)
  })
  end_time <- Sys.time()
  
  if (!is.null(result)) {
    exec_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    eval <- evaluate_comprehensive(result, test_data$true_labels, test_data$x_list, 
                                   paste0("g", g))
    all_results[[length(all_results)+1]] <- eval
    cat(sprintf(" ✓ Time: %.2f sec, Iter: %d, Noise: %.2f\n", 
                exec_time, result$iterations, eval$noise_proportion_estimated))
  }
}

# Test 7: Stress tests - pathological data
cat("\n--- Test 7: Stress Tests with Pathological Data ---\n")

# 7.1 Collinear data
cat("  Collinear data...")
rows <- 5
cols <- 4
x_list <- vector("list", 100)
for (i in 1:100) {
  base <- matrix(rnorm(rows * cols, 0, 1), rows, cols)
  for (j in 2:cols) {
    base[, j] <- base[, 1] + rnorm(rows, 0, 0.01)
  }
  x_list[[i]] <- base
}

result <- tryCatch({
  matrix_variate_noise_fit(
    x_list = x_list,
    g = 2,
    noise_type = "hc",
    estimate_k = TRUE,
    verbose = FALSE,
    nstart = 10,
    max_iter = 100
  )
}, error = function(e) {
  cat(" ✗ Error\n")
  return(NULL)
})

if (!is.null(result)) {
  all_results[[length(all_results)+1]] <- data.frame(
    test_name = "stress_collinear",
    selected_k = if(!is.null(result$k_selection)) result$k_selection$selected_k else NA,
    noise_proportion_estimated = result$noise$pi,
    noise_proportion_true = NA,
    noise_accuracy = NA,
    noise_f1 = NA,
    ari = NA,
    iterations = result$iterations,
    converged = result$converged,
    ks_pvalue = NA,
    stringsAsFactors = FALSE
  )
  cat(sprintf(" ✓ Converged: %s, Iter: %d, Noise: %.3f\n", 
              result$converged, result$iterations, result$noise$pi))
}

# 7.2 Sparse data
cat("  Sparse data (80% zeros)...")
x_list <- vector("list", 100)
for (i in 1:100) {
  mat <- matrix(rnorm(rows * cols, 0, 0.5), rows, cols)
  mat[sample(1:(rows*cols), 0.8*rows*cols)] <- 0
  x_list[[i]] <- mat
}

result <- tryCatch({
  matrix_variate_noise_fit(
    x_list = x_list,
    g = 2,
    noise_type = "hc",
    estimate_k = TRUE,
    verbose = FALSE,
    nstart = 10,
    max_iter = 100
  )
}, error = function(e) {
  cat(" ✗ Error\n")
  return(NULL)
})

if (!is.null(result)) {
  all_results[[length(all_results)+1]] <- data.frame(
    test_name = "stress_sparse",
    selected_k = if(!is.null(result$k_selection)) result$k_selection$selected_k else NA,
    noise_proportion_estimated = result$noise$pi,
    noise_proportion_true = NA,
    noise_accuracy = NA,
    noise_f1 = NA,
    ari = NA,
    iterations = result$iterations,
    converged = result$converged,
    ks_pvalue = NA,
    stringsAsFactors = FALSE
  )
  cat(sprintf(" ✓ Converged: %s, Iter: %d, Noise: %.3f\n", 
              result$converged, result$iterations, result$noise$pi))
}

# 7.3 Heavy-tailed data
cat("  Heavy-tailed (t-distribution)...")
x_list <- vector("list", 100)
for (i in 1:100) {
  x_list[[i]] <- matrix(rt(rows * cols, df = 2), rows, cols)
}

result <- tryCatch({
  matrix_variate_noise_fit(
    x_list = x_list,
    g = 2,
    noise_type = "hc",
    estimate_k = TRUE,
    verbose = FALSE,
    nstart = 10,
    max_iter = 100
  )
}, error = function(e) {
  cat(" ✗ Error\n")
  return(NULL)
})

if (!is.null(result)) {
  all_results[[length(all_results)+1]] <- data.frame(
    test_name = "stress_heavytail",
    selected_k = if(!is.null(result$k_selection)) result$k_selection$selected_k else NA,
    noise_proportion_estimated = result$noise$pi,
    noise_proportion_true = NA,
    noise_accuracy = NA,
    noise_f1 = NA,
    ari = NA,
    iterations = result$iterations,
    converged = result$converged,
    ks_pvalue = NA,
    stringsAsFactors = FALSE
  )
  cat(sprintf(" ✓ Converged: %s, Iter: %d, Noise: %.3f\n", 
              result$converged, result$iterations, result$noise$pi))
}

# ============================================================================
# SAVE AND SUMMARIZE RESULTS
# ============================================================================

cat("\n")
cat(rep("=", 80), collapse = "")
cat("\nSAVING RESULTS\n")
cat(rep("=", 80), collapse = "")
cat("\n")

if (length(all_results) > 0) {
  # Combine all results
  final_df <- do.call(rbind, all_results)
  
  # Save to CSV
  write.csv(final_df, "extensive_test_results.csv", row.names = FALSE)
  cat("\n✓ Results saved to: extensive_test_results.csv\n")
  
  # Print summary statistics
  cat("\n")
  cat(rep("=", 80), collapse = "")
  cat("\nSUMMARY STATISTICS\n")
  cat(rep("=", 80), collapse = "")
  cat("\n")
  
  # Filter numeric columns for summary
  numeric_cols <- sapply(final_df, is.numeric)
  summary_stats <- data.frame(
    Metric = names(final_df)[numeric_cols],
    Mean = sapply(final_df[numeric_cols], mean, na.rm = TRUE),
    SD = sapply(final_df[numeric_cols], sd, na.rm = TRUE),
    Min = sapply(final_df[numeric_cols], min, na.rm = TRUE),
    Max = sapply(final_df[numeric_cols], max, na.rm = TRUE)
  )
  
  print(summary_stats)
  
  # Key performance indicators
  cat("\n")
  cat(rep("=", 80), collapse = "")
  cat("\nKEY PERFORMANCE INDICATORS\n")
  cat(rep("=", 80), collapse = "")
  cat("\n")
  
  # Noise detection accuracy (where ground truth available)
  noise_tests <- final_df[!is.na(final_df$noise_proportion_true) & final_df$noise_proportion_true > 0, ]
  if (nrow(noise_tests) > 0) {
    noise_error <- mean(abs(noise_tests$noise_proportion_estimated - noise_tests$noise_proportion_true))
    cat(sprintf("  Average Noise Detection Error: %.3f (%.1f%%)\n", 
                noise_error, noise_error * 100))
  }
  
  # Convergence rate
  conv_rate <- mean(final_df$converged, na.rm = TRUE) * 100
  cat(sprintf("  Convergence Rate: %.1f%%\n", conv_rate))
  
  # Average iterations
  avg_iter <- mean(final_df$iterations, na.rm = TRUE)
  cat(sprintf("  Average Iterations to Converge: %.1f\n", avg_iter))
  
  # Average ARI where available
  ari_tests <- final_df[!is.na(final_df$ari), ]
  if (nrow(ari_tests) > 0) {
    avg_ari <- mean(ari_tests$ari, na.rm = TRUE)
    cat(sprintf("  Average ARI (Cluster Recovery): %.3f\n", avg_ari))
  }
  
  # K selection range
  k_tests <- final_df[!is.na(final_df$selected_k) & final_df$selected_k > 0, ]
  if (nrow(k_tests) > 0) {
    cat(sprintf("  Selected k range: %.2e to %.2e\n", 
                min(k_tests$selected_k), max(k_tests$selected_k)))
  }
}

cat("\n")
cat(rep("=", 80), collapse = "")
cat("\nEXTENSIVE TESTS COMPLETE\n")
cat(rep("=", 80), collapse = "")
cat("\n\n")