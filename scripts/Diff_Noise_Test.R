# test_noise_types.R
# Tests the algorithm against different types of noise structures

source("..R/Matrix_Init.R")
source("..R/Matrix.R")
source("..R/Matrix_KS_Score.R")
source("..R/Matrix_Noise.R")
source("..R/Matrix_Noise_BR.R")
source("..R/Matrix_Utils.R")

# Generate data with different noise patterns
generate_noisy_data <- function(n_clean_per_group, n_noise, rows, cols,
                                group_means, row_covs, col_covs,
                                noise_pattern = "uniform") {
  g <- length(group_means)
  total_n <- sum(n_clean_per_group) + n_noise
  
  # Generate clean data
  clean_data <- generate_test_data(
    n_per_group = n_clean_per_group,
    rows = rows, cols = cols,
    group_means = group_means,
    row_covs = row_covs,
    col_covs = col_covs,
    noise_prop = 0
  )
  
  x_list <- clean_data$x_list
  true_labels <- clean_data$true_labels
  
  # Add noise with specific pattern
  all_data <- do.call(rbind, lapply(x_list, as.vector))
  data_mean <- mean(all_data, na.rm = TRUE)
  data_sd <- sd(all_data, na.rm = TRUE)
  data_range <- range(all_data, na.rm = TRUE)
  
  for (i in 1:n_noise) {
    if (noise_pattern == "uniform") {
      # Uniform noise across entire range
      noise_mat <- matrix(runif(rows * cols, data_range[1], data_range[2]), rows, cols)
      
    } else if (noise_pattern == "gaussian") {
      # Gaussian noise with same mean/variance as data
      noise_mat <- matrix(rnorm(rows * cols, data_mean, data_sd), rows, cols)
      
    } else if (noise_pattern == "outlier") {
      # Extreme outliers
      extreme_range <- data_range + c(-10, 10) * data_sd
      noise_mat <- matrix(runif(rows * cols, extreme_range[1], extreme_range[2]), rows, cols)
      
    } else if (noise_pattern == "row_structured") {
      # Row-specific noise (some rows have different distributions)
      noise_mat <- matrix(0, rows, cols)
      for (row in 1:rows) {
        if (row %% 2 == 0) {
          noise_mat[row, ] <- runif(cols, data_range[1] - data_sd, data_range[2] - data_sd)
        } else {
          noise_mat[row, ] <- runif(cols, data_range[1] + data_sd, data_range[2] + data_sd)
        }
      }
      
    } else if (noise_pattern == "col_structured") {
      # Column-specific noise
      noise_mat <- matrix(0, rows, cols)
      for (col in 1:cols) {
        if (col %% 2 == 0) {
          noise_mat[, col] <- runif(rows, data_range[1] - data_sd, data_range[2] - data_sd)
        } else {
          noise_mat[, col] <- runif(rows, data_range[1] + data_sd, data_range[2] + data_sd)
        }
      }
      
    } else if (noise_pattern == "block") {
      # Block-structured noise
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
  # Calculate adjusted Rand index for non-noise points
  non_noise_idx <- which(true_labels > 0 & result$cluster > 0)
  
  ari <- if (length(non_noise_idx) > 1) {
    # Simple Rand index approximation
    true_non_noise <- true_labels[non_noise_idx]
    pred_non_noise <- result$cluster[non_noise_idx]
    
    # Count agreements
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
  
  list(
    ari = ari,
    noise_accuracy = noise_accuracy,
    noise_precision = noise_precision,
    noise_recall = noise_recall,
    noise_proportion_estimated = result$noise$pi,
    noise_proportion_true = mean(true_noise),
    logLik = tail(result$logLik, 1)
  )
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
  
  cat("\n=== Testing Different Noise Patterns ===\n")
  cat(sprintf("%-20s %10s %10s %10s %10s %10s\n", 
              "Pattern", "ARI", "NoiseAcc", "Precision", "Recall", "NoiseProp"))
  cat(rep("-", 70), "\n", sep="")
  
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
    
    # Test both HC and BR noise models
    for (noise_type in c("hc", "br")) {
      result <- tryCatch({
        matrix_variate_noise_fit(
          x_list = test_data$x_list,
          g = g,
          noise_type = noise_type,
          max_iter = 200,
          tol = 1e-5,
          nstart = 10,
          noise_k = if(noise_type == "hc") 1e-5 else NULL,
          estimate_k = (noise_type == "hc"),
          verbose = FALSE
        )
      }, error = function(e) {
        cat("\nError in", pattern, "with", noise_type, ":", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(result)) {
        eval <- evaluate_fit(result, test_data$true_labels, test_data$x_list)
        
        if (noise_type == "hc") {
          cat(sprintf(" %10.3f %10.3f %10.3f %10.3f %10.3f\n",
                      eval$ari, eval$noise_accuracy, 
                      ifelse(is.na(eval$noise_precision), 0, eval$noise_precision),
                      ifelse(is.na(eval$noise_recall), 0, eval$noise_recall),
                      eval$noise_proportion_estimated))
          results[[paste(pattern, noise_type)]] <- eval
        }
      }
    }
  }
  
  return(results)
}

# Run tests
test_results <- test_noise_patterns()

# Compare HC vs BR performance
cat("\n=== HC vs BR Noise Model Comparison ===\n")
# Note: Add specific comparison logic here