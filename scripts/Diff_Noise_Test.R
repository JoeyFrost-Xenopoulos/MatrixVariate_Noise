# test_noise_types.R
# Standalone test for different noise patterns with CSV export

source("../R/Matrix_Init.R")
source("../R/Matrix.R")
source("../R/Matrix_KS_Score.R")
source("../R/Matrix_Noise.R")
source("../R/Matrix_Noise_BR.R")
source("../R/Matrix_Utils.R")

# Helper to generate clean data
generate_clean_data <- function(n_per_group, rows, cols, group_means, row_covs, col_covs) {
  g <- length(group_means)
  total_n <- sum(n_per_group)
  x_list <- vector("list", total_n)
  true_labels <- integer(total_n)
  
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
  
  list(x_list = x_list, true_labels = true_labels)
}

# Generate data with different noise patterns
generate_noisy_data <- function(n_clean_per_group, n_noise, rows, cols,
                                group_means, row_covs, col_covs,
                                noise_pattern = "uniform") {
  
  # Generate clean data
  clean_result <- generate_clean_data(n_clean_per_group, rows, cols, 
                                      group_means, row_covs, col_covs)
  x_list <- clean_result$x_list
  true_labels <- clean_result$true_labels
  
  # Get data statistics for noise generation
  all_data <- do.call(rbind, lapply(x_list, as.vector))
  data_mean <- mean(all_data, na.rm = TRUE)
  data_sd <- sd(all_data, na.rm = TRUE)
  data_range <- range(all_data, na.rm = TRUE)
  
  # Add noise with specific pattern
  for (i in 1:n_noise) {
    if (noise_pattern == "uniform") {
      noise_mat <- matrix(runif(rows * cols, data_range[1], data_range[2]), rows, cols)
      
    } else if (noise_pattern == "gaussian") {
      noise_mat <- matrix(rnorm(rows * cols, data_mean, data_sd), rows, cols)
      
    } else if (noise_pattern == "outlier") {
      extreme_range <- data_range + c(-10, 10) * data_sd
      noise_mat <- matrix(runif(rows * cols, extreme_range[1], extreme_range[2]), rows, cols)
      
    } else if (noise_pattern == "row_structured") {
      noise_mat <- matrix(0, rows, cols)
      for (row in 1:rows) {
        if (row %% 2 == 0) {
          noise_mat[row, ] <- runif(cols, data_range[1] - data_sd, data_range[2] - data_sd)
        } else {
          noise_mat[row, ] <- runif(cols, data_range[1] + data_sd, data_range[2] + data_sd)
        }
      }
      
    } else if (noise_pattern == "col_structured") {
      noise_mat <- matrix(0, rows, cols)
      for (col in 1:cols) {
        if (col %% 2 == 0) {
          noise_mat[, col] <- runif(rows, data_range[1] - data_sd, data_range[2] - data_sd)
        } else {
          noise_mat[, col] <- runif(rows, data_range[1] + data_sd, data_range[2] + data_sd)
        }
      }
      
    } else if (noise_pattern == "block") {
      noise_mat <- matrix(rnorm(rows * cols, data_mean, data_sd), rows, cols)
      block_rows <- 1:floor(rows/2)
      block_cols <- 1:floor(cols/2)
      noise_mat[block_rows, block_cols] <- noise_mat[block_rows, block_cols] * 5
      
    } else {
      stop("Unknown noise pattern")
    }
    
    x_list[[length(x_list) + 1]] <- noise_mat
    true_labels <- c(true_labels, 0)
  }
  
  list(x_list = x_list, true_labels = true_labels)
}

# Helper function to evaluate fit quality
evaluate_fit <- function(result, true_labels, x_list) {
  # Calculate simple accuracy for non-noise points
  non_noise_idx <- which(true_labels > 0 & result$cluster > 0)
  
  ari <- if (length(non_noise_idx) > 1) {
    true_non_noise <- true_labels[non_noise_idx]
    pred_non_noise <- result$cluster[non_noise_idx]
    
    # Simple Rand index approximation
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
    if (total_pairs > 0) agreements / total_pairs else 1
  } else {
    1
  }
  
  # Noise detection accuracy
  true_noise <- true_labels == 0
  pred_noise <- result$cluster == 0
  
  noise_accuracy <- mean(true_noise == pred_noise)
  noise_precision <- if (sum(pred_noise) > 0) sum(true_noise & pred_noise) / sum(pred_noise) else NA
  noise_recall <- if (sum(true_noise) > 0) sum(true_noise & pred_noise) / sum(true_noise) else NA
  
  # Get selected k if available
  selected_k <- if (!is.null(result$k_selection)) result$k_selection$selected_k else NA
  
  list(
    pattern = NA,  # Will be filled in
    ari = ari,
    noise_accuracy = noise_accuracy,
    noise_precision = noise_precision,
    noise_recall = noise_recall,
    noise_proportion_estimated = result$noise$pi,
    noise_proportion_true = mean(true_noise),
    selected_k = selected_k,
    iterations = result$iterations,
    converged = result$converged,
    logLik = tail(result$logLik, 1)
  )
}

# Save results to CSV
save_noise_test_results <- function(results_list, filename) {
  if (length(results_list) == 0) {
    cat("No results to save\n")
    return(NULL)
  }
  
  df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(
      pattern = x$pattern,
      ari = x$ari,
      noise_accuracy = x$noise_accuracy,
      noise_precision = ifelse(is.na(x$noise_precision), 0, x$noise_precision),
      noise_recall = ifelse(is.na(x$noise_recall), 0, x$noise_recall),
      noise_proportion_estimated = x$noise_proportion_estimated,
      noise_proportion_true = x$noise_proportion_true,
      selected_k = ifelse(is.na(x$selected_k), 0, x$selected_k),
      iterations = x$iterations,
      converged = x$converged,
      logLik = x$logLik,
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(df, filename, row.names = FALSE)
  cat("Results saved to:", filename, "\n")
  return(df)
}

# Main test function
test_noise_patterns <- function() {
  set.seed(42)
  
  rows <- 5
  cols <- 4
  g <- 2
  n_clean_per_group <- c(50, 50)
  n_noise <- 20
  
  group_means <- list(
    matrix(0, rows, cols),
    matrix(3, rows, cols)
  )
  
  row_covs <- list(diag(rows), diag(rows))
  col_covs <- list(diag(cols), diag(cols))
  
  noise_patterns <- c("uniform", "gaussian", "outlier", "row_structured", "col_structured", "block")
  results <- list()
  
  cat("\n========================================\n")
  cat("TEST 2: Different Noise Patterns\n")
  cat("========================================\n")
  cat(sprintf("\n%-20s %10s %10s %10s %10s %10s %12s\n", 
              "Pattern", "ARI", "NoiseAcc", "Precision", "Recall", "NoiseProp", "SelectedK"))
  cat(rep("-", 95), "\n", sep="")
  
  for (pattern in noise_patterns) {
    cat(sprintf("%-20s", pattern))
    
    # Generate data with this noise pattern
    test_data <- generate_noisy_data(
      n_clean_per_group = n_clean_per_group,
      n_noise = n_noise,
      rows = rows, cols = cols,
      group_means = group_means,
      row_covs = row_covs,
      col_covs = col_covs,
      noise_pattern = pattern
    )
    
    # Test HC noise model
    result <- tryCatch({
      matrix_variate_noise_fit(
        x_list = test_data$x_list,
        g = g,
        noise_type = "hc",
        max_iter = 200,
        tol = 1e-5,
        nstart = 10,
        estimate_k = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      cat("\n  Error in", pattern, ":", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      eval_result <- evaluate_fit(result, test_data$true_labels, test_data$x_list)
      eval_result$pattern <- pattern
      results[[pattern]] <- eval_result
      
      # Format selected k for display
      selected_k_str <- if (!is.na(eval_result$selected_k) && eval_result$selected_k > 0) 
        sprintf("%.2e", eval_result$selected_k) else "N/A"
      
      cat(sprintf(" %10.3f %10.3f %10.3f %10.3f %10.3f %12s\n",
                  eval_result$ari, 
                  eval_result$noise_accuracy, 
                  ifelse(is.na(eval_result$noise_precision), 0, eval_result$noise_precision),
                  ifelse(is.na(eval_result$noise_recall), 0, eval_result$noise_recall),
                  eval_result$noise_proportion_estimated,
                  selected_k_str))
    } else {
      cat(" %10s %10s %10s %10s %10s %12s\n", "FAILED", "FAILED", "FAILED", "FAILED", "FAILED", "FAILED")
    }
  }
  
  # Save results to CSV
  if (length(results) > 0) {
    save_noise_test_results(results, "noise_pattern_test_results.csv")
  }
  
  # Summary
  cat("\n--- Summary of Results ---\n")
  for (pattern in names(results)) {
    cat(sprintf("%s: ARI=%.3f, NoiseAcc=%.3f, NoiseProp=%.3f, SelectedK=%.2e\n",
                pattern, 
                results[[pattern]]$ari,
                results[[pattern]]$noise_accuracy,
                results[[pattern]]$noise_proportion_estimated,
                results[[pattern]]$selected_k))
  }
  
  return(results)
}

# Run the tests
test_results <- test_noise_patterns()
cat("\n=== Noise Pattern Testing Complete ===\n")