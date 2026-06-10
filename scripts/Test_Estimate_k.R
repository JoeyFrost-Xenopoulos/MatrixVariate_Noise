# test_adaptive_k_selection.R
# Comprehensive testing framework for adaptive k selection in matrix-variate noise clustering

library(ggplot2)
library(reshape2)
library(dplyr)

# Try to load mclust for ARI, but don't fail if not available
if (!require(mclust, quietly = TRUE)) {
  warning("mclust package not available. ARI calculations will return NA.")
  # Provide a fallback
  adjustedRandIndex <- function(x, y) {
    return(NA)
  }
}

#' Source all necessary R scripts for the matrix-variate noise clustering
#'
#' @param script_dir Directory containing the R scripts. If NULL, assumes scripts
#'   are in the current working directory.
#' @return Logical indicating success
source_all_scripts <- function(script_dir = NULL) {

  # Prefer the project's R/ directory and known filename patterns
  script_files <- c(
    "Matrix_Utils.R",
    "Matrix.R",
    "Matrix_Noise.R",
    "Matrix_Noise_BR.R",
    "Matrix_KS_Score.R",
    "Matrix_Init.R"
  )

  # Auto-detect a sensible R/ directory when script_dir is not provided
  if (is.null(script_dir)) {
    candidate_dirs <- c("R", "./R", "../R")
    found <- NULL
    for (d in candidate_dirs) {
      if (dir.exists(d)) { found <- d; break }
    }
    if (!is.null(found)) {
      script_files <- file.path(found, script_files)
    } else {
      # Fall back to current working directory files if no R/ directory
      script_files <- script_files
    }
  } else {
    script_files <- file.path(script_dir, script_files)
  }

  # Source each script if present; otherwise try to source any nearby R files
  sourced_any <- FALSE
  for (script in script_files) {
    if (file.exists(script)) {
      cat(sprintf("Sourcing: %s\n", basename(script)))
      tryCatch({
        source(script)
        sourced_any <- TRUE
      }, error = function(e) {
        stop(sprintf("Error sourcing %s: %s\n", script, e$message))
      })
    }
  }

  if (!sourced_any) {
    cat("No expected script files found; attempting to source any R files in the project root or script directory.\n")
    r_files <- list.files(path = ifelse(is.null(script_dir), ".", script_dir), pattern = "\\.R$", recursive = FALSE)
    if (length(r_files) > 0) {
      for (rf in r_files) {
        cat(sprintf("  Sourcing: %s\n", rf))
        tryCatch({
          source(rf)
        }, error = function(e) {
          cat(sprintf("    Warning sourcing %s: %s\n", rf, e$message))
        })
      }
    } else {
      stop("No R script files found. Please ensure the algorithm scripts are in the project's R/ directory or working directory.")
    }
  }

  # Verify key functions are available
  required_functions <- c(
    "matrix_variate_noise_fit",
    "matrix_validate_x_list",
    "matrix_log_sum_exp",
    "matrix_variate_log_density"
  )

  missing_functions <- required_functions[!sapply(required_functions, exists)]

  if (length(missing_functions) > 0) {
    stop(sprintf("Missing required functions after sourcing: %s\n", paste(missing_functions, collapse = ", ")))
  }

  cat("All required functions loaded successfully!\n")
  return(TRUE)
}

# If the original functions are in a single file, you can source that instead
# source("matrix_variate_noise_clustering.R")  # Uncomment and adjust path as needed

# Otherwise, source all individual scripts
# Try to source from current directory or specify the correct path
if (!exists("matrix_variate_noise_fit")) {
  cat("Loading matrix-variate noise clustering functions...\n")
  
  # Option 1: If all functions are in one file
  if (file.exists("matrix_variate_noise_clustering.R")) {
    source("matrix_variate_noise_clustering.R")
  } 
  # Option 2: If functions are split across multiple files
  else {
    # Auto-detect the project's R/ directory and source files from there
    source_all_scripts(NULL)
  }
}

#' Generate synthetic matrix-variate data with known structure
#'
#' @param n_per_group Number of observations per Gaussian component
#' @param noise_prop Proportion of noise observations (0-1)
#' @param dims Vector of length 2: c(rows, cols)
#' @param group_means List of mean matrices for each Gaussian component
#' @param row_covs List of row covariance matrices
#' @param col_covs List of column covariance matrices
#' @param noise_type Type of noise to generate ("uniform" or "outlier")
#' @return List containing data, true labels, and parameters
generate_test_data <- function(n_per_group = 50,
                               noise_prop = 0.1,
                               dims = c(5, 5),
                               n_groups = 3,
                               group_means = NULL,
                               row_covs = NULL,
                               col_covs = NULL,
                               noise_type = "uniform") {
  
  r <- dims[1]
  p <- dims[2]
  n_clean <- n_per_group * n_groups
  n_noise <- round(n_clean * noise_prop / (1 - noise_prop))
  n_total <- n_clean + n_noise
  
  # Default parameters if not provided
  if (is.null(group_means)) {
    group_means <- list()
    for (k in 1:n_groups) {
      # Create distinct mean patterns
      M <- matrix(rnorm(r * p, mean = k * 2, sd = 1), r, p)
      group_means[[k]] <- M
    }
  }
  
  if (is.null(row_covs)) {
    row_covs <- list()
    for (k in 1:n_groups) {
      # AR(1) structure for row covariance
      rho <- 0.7
      R <- matrix(NA, r, r)
      for (i in 1:r) {
        for (j in 1:r) {
          R[i, j] <- rho^abs(i - j)
        }
      }
      row_covs[[k]] <- R
    }
  }
  
  if (is.null(col_covs)) {
    col_covs <- list()
    for (k in 1:n_groups) {
      # Compound symmetry for column covariance
      sigma <- 1
      rho <- 0.5
      C <- matrix(rho, p, p)
      diag(C) <- sigma
      col_covs[[k]] <- C
    }
  }
  
  # Generate clean data
  data_list <- list()
  true_labels <- integer(n_total)
  component_assignments <- integer(n_total)
  
  idx <- 1
  for (group in 1:n_groups) {
    for (i in 1:n_per_group) {
      # Generate matrix-variate normal
      Z <- matrix(rnorm(r * p), r, p)
      # Transform using row and column covariances
      R_chol <- chol(row_covs[[group]])
      C_chol <- chol(col_covs[[group]])
      X <- group_means[[group]] + t(R_chol) %*% Z %*% C_chol
      data_list[[idx]] <- X
      true_labels[idx] <- group
      component_assignments[idx] <- group
      idx <- idx + 1
    }
  }
  
  # Generate noise data
  for (i in 1:n_noise) {
    if (noise_type == "uniform") {
      # Uniform noise within a bounding box
      all_clean <- do.call(rbind, lapply(data_list[1:n_clean], as.vector))
      bounds <- apply(all_clean, 2, range)
      noise_vec <- runif(r * p, min = bounds[1, ], max = bounds[2, ])
      X <- matrix(noise_vec, r, p)
    } else if (noise_type == "outlier") {
      # Outlier noise far from means
      X <- matrix(rnorm(r * p, mean = 100, sd = 10), r, p)
    } else {
      stop("Unknown noise_type")
    }
    data_list[[idx]] <- X
    true_labels[idx] <- 0  # 0 indicates noise
    component_assignments[idx] <- 0
    idx <- idx + 1
  }
  
  # Shuffle data
  shuffle_idx <- sample(n_total)
  data_list <- data_list[shuffle_idx]
  true_labels <- true_labels[shuffle_idx]
  component_assignments <- component_assignments[shuffle_idx]
  
  list(
    data = data_list,
    true_labels = true_labels,
    true_components = component_assignments,
    parameters = list(
      n_groups = n_groups,
      n_per_group = n_per_group,
      n_noise = n_noise,
      dims = dims,
      group_means = group_means,
      row_covs = row_covs,
      col_covs = col_covs
    )
  )
}

#' Evaluate clustering performance
#'
#' @param fit Fitted model from matrix_variate_noise_fit
#' @param truth True labels (0 for noise, 1..K for groups)
#' @return List of evaluation metrics
evaluate_clustering <- function(fit, truth) {
  n <- length(truth)
  predicted <- fit$cluster
  
  # 1. Noise identification metrics
  true_noise <- truth == 0
  pred_noise <- predicted == 0
  
  noise_accuracy <- sum(true_noise == pred_noise) / n
  noise_precision <- ifelse(sum(pred_noise) > 0, 
                            sum(pred_noise & true_noise) / sum(pred_noise), 
                            0)
  noise_recall <- ifelse(sum(true_noise) > 0,
                         sum(pred_noise & true_noise) / sum(true_noise),
                         0)
  noise_f1 <- ifelse(noise_precision + noise_recall > 0,
                     2 * noise_precision * noise_recall / (noise_precision + noise_recall),
                     0)
  
  # 2. Group parsing metrics (only for non-noise points)
  clean_idx <- truth > 0
  if (sum(clean_idx) > 0) {
    truth_clean <- truth[clean_idx]
    pred_clean <- predicted[clean_idx]
    
    # Adjusted Rand Index for groups
    ari <- tryCatch(
      mclust::adjustedRandIndex(truth_clean, pred_clean),
      error = function(e) NA
    )
    
    # Number of groups found
    n_groups_found <- length(unique(pred_clean[pred_clean > 0]))
    
    # Group purity
    if (n_groups_found > 0) {
      group_purity <- 0
      for (g in unique(pred_clean)) {
        if (g > 0) {
          mask <- pred_clean == g
          if (sum(mask) > 0) {
            true_in_group <- truth_clean[mask]
            if (length(unique(true_in_group)) > 0) {
              majority <- max(table(true_in_group))
              group_purity <- group_purity + majority / sum(mask)
            }
          }
        }
      }
      group_purity <- group_purity / n_groups_found
    } else {
      group_purity <- NA
    }
  } else {
    ari <- NA
    n_groups_found <- 0
    group_purity <- NA
  }
  
  # 3. Overall assignment accuracy
  # Need to match predicted group labels to true groups
  if (sum(clean_idx) > 0 && n_groups_found == length(unique(truth_clean))) {
    # Create confusion matrix and find best label mapping
    conf_matrix <- table(truth_clean, pred_clean)
    assignment_accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
  } else {
    assignment_accuracy <- NA
  }
  
  # 4. Parameter estimation accuracy (if true parameters known)
  param_accuracy <- list()
  # This would require storing true parameters in fit object
  
  list(
    noise_identification = list(
      accuracy = noise_accuracy,
      precision = noise_precision,
      recall = noise_recall,
      f1 = noise_f1
    ),
    group_parsing = list(
      ari = ari,
      n_groups_found = n_groups_found,
      group_purity = group_purity
    ),
    assignment = list(
      overall_accuracy = assignment_accuracy
    ),
    convergence = list(
      converged = fit$converged,
      iterations = fit$iterations,
      final_loglik = tail(fit$logLik, 1)
    ),
    k_selection = if (!is.null(fit$k_selection)) {
      list(
        selected_k = fit$k_selection$selected_k,
        best_ks = min(fit$k_selection$ks_scores),
        grid_size = length(fit$k_selection$k_grid)
      )
    } else {
      NULL
    }
  )
}

#' Run comprehensive tests
#'
#' @param test_configs List of test configurations
#' @param n_replicates Number of replicates per configuration
#' @return Data frame with all results
run_adaptive_k_tests <- function(test_configs, n_replicates = 10) {
  
  results <- list()
  
  for (config_idx in seq_along(test_configs)) {
    config <- test_configs[[config_idx]]
    
    cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
    cat(sprintf("Testing configuration %d/%d\n", config_idx, length(test_configs)))
    cat(sprintf("  Groups: %d, N per group: %d, Noise proportion: %.2f\n",
                config$n_groups, config$n_per_group, config$noise_prop))
    cat(sprintf("  Dimensions: %d x %d\n", config$dims[1], config$dims[2]))
    cat(sprintf("  Noise type: %s\n", config$noise_type))
    cat(sprintf("  Adaptive grid: %s\n", config$adaptive_grid))
    cat(sprintf("  Replicates: %d\n", n_replicates))
    
    config_results <- list()
    
    for (rep in 1:n_replicates) {
      if (rep %% 5 == 0) cat(sprintf("    Replicate %d/%d\n", rep, n_replicates))
      
      # Generate data
      test_data <- generate_test_data(
        n_per_group = config$n_per_group,
        noise_prop = config$noise_prop,
        dims = config$dims,
        n_groups = config$n_groups,
        noise_type = config$noise_type
      )
      
      # Fit model with adaptive k selection
      tryCatch({
        fit <- matrix_variate_noise_fit(
          x_list = test_data$data,
          g = config$n_groups,  # True number of groups
          noise_type = "hc",
          max_iter = 500,
          tol = 1e-6,
          nstart = 20,
          estimate_k = TRUE,
          adaptive_grid = config$adaptive_grid,
          k_grid = config$k_grid,
          verbose = FALSE
        )
        
        # Evaluate
        eval <- evaluate_clustering(fit, test_data$true_labels)
        
        # Store results
        config_results[[rep]] <- list(
          config_id = config_idx,
          replicate = rep,
          n_groups_true = config$n_groups,
          n_per_group = config$n_per_group,
          noise_prop = config$noise_prop,
          dim_r = config$dims[1],
          dim_p = config$dims[2],
          noise_type = config$noise_type,
          adaptive_grid = config$adaptive_grid,
          noise_accuracy = eval$noise_identification$accuracy,
          noise_precision = eval$noise_identification$precision,
          noise_recall = eval$noise_identification$recall,
          noise_f1 = eval$noise_identification$f1,
          group_ari = eval$group_parsing$ari,
          n_groups_found = eval$group_parsing$n_groups_found,
          group_purity = eval$group_parsing$group_purity,
          assignment_accuracy = eval$assignment$overall_accuracy,
          converged = eval$convergence$converged,
          iterations = eval$convergence$iterations,
          final_loglik = eval$convergence$final_loglik,
          selected_k = if (!is.null(eval$k_selection)) eval$k_selection$selected_k else NA,
          ks_best = if (!is.null(eval$k_selection)) eval$k_selection$best_ks else NA
        )
        
        cat(sprintf("    Replicate %d: F1=%.3f, ARI=%.3f, Groups=%d, k=%.2e\n",
                    rep, eval$noise_identification$f1, 
                    ifelse(is.na(eval$group_parsing$ari), 0, eval$group_parsing$ari),
                    eval$group_parsing$n_groups_found,
                    if (!is.null(eval$k_selection)) eval$k_selection$selected_k else config$noise_k))
        
      }, error = function(e) {
        cat(sprintf("      Error in replicate %d: %s\n", rep, e$message))
        config_results[[rep]] <- NULL
      })
    }
    
    # Combine results for this configuration
    config_results <- config_results[!sapply(config_results, is.null)]
    if (length(config_results) > 0) {
      results <- c(results, config_results)
    }
  }
  
  # Convert to data frame
  if (length(results) > 0) {
    results_df <- do.call(rbind, lapply(results, as.data.frame))
    return(results_df)
  } else {
    warning("No successful results to return")
    return(data.frame())
  }
}

#' Visualize test results
#'
#' @param results_df Data frame from run_adaptive_k_tests
visualize_results <- function(results_df) {
  
  if (nrow(results_df) == 0) {
    cat("No results to visualize\n")
    return()
  }
  
  # 1. Noise identification performance
  p1 <- ggplot(results_df, aes(x = factor(noise_prop), y = noise_f1, fill = factor(adaptive_grid))) +
    geom_boxplot() +
    facet_wrap(~ n_groups_true + dim_r + dim_p, labeller = label_both) +
    labs(title = "Noise Identification F1 Score",
         x = "Noise Proportion", y = "F1 Score",
         fill = "Adaptive Grid") +
    theme_minimal()
  
  # 2. Group parsing performance
  p2 <- ggplot(results_df, aes(x = factor(noise_prop), y = group_ari, fill = factor(adaptive_grid))) +
    geom_boxplot() +
    facet_wrap(~ n_groups_true + dim_r + dim_p, labeller = label_both) +
    labs(title = "Group Parsing (Adjusted Rand Index)",
         x = "Noise Proportion", y = "ARI",
         fill = "Adaptive Grid") +
    theme_minimal()
  
  # 3. Number of groups found
  p3 <- ggplot(results_df, aes(x = factor(noise_prop), y = n_groups_found, fill = factor(adaptive_grid))) +
    geom_boxplot() +
    geom_hline(aes(yintercept = unique(n_groups_true), color = "True"), linetype = "dashed") +
    facet_wrap(~ n_groups_true + dim_r + dim_p, labeller = label_both) +
    labs(title = "Number of Groups Found",
         x = "Noise Proportion", y = "Number of Groups",
         fill = "Adaptive Grid", color = "") +
    theme_minimal()
  
  # 4. Convergence and iterations
  p4 <- ggplot(results_df, aes(x = factor(noise_prop), y = iterations, fill = factor(adaptive_grid))) +
    geom_boxplot() +
    facet_wrap(~ n_groups_true + dim_r + dim_p, labeller = label_both) +
    labs(title = "EM Iterations to Convergence",
         x = "Noise Proportion", y = "Iterations",
         fill = "Adaptive Grid") +
    theme_minimal()
  
  # 5. Selected k values (log scale)
  if (any(!is.na(results_df$selected_k))) {
    p5 <- ggplot(results_df, aes(x = factor(noise_prop), y = log10(selected_k), fill = factor(adaptive_grid))) +
      geom_boxplot() +
      facet_wrap(~ n_groups_true + dim_r + dim_p, labeller = label_both) +
      labs(title = "Selected k Values (log10 scale)",
           x = "Noise Proportion", y = "log10(k)",
           fill = "Adaptive Grid") +
      theme_minimal()
    print(p5)
  }
  
  # Print plots
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
  # Summary statistics
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("SUMMARY STATISTICS\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  summary_stats <- results_df %>%
    group_by(n_groups_true, n_per_group, noise_prop, dim_r, dim_p, adaptive_grid) %>%
    summarise(
      mean_noise_f1 = mean(noise_f1, na.rm = TRUE),
      sd_noise_f1 = sd(noise_f1, na.rm = TRUE),
      mean_group_ari = mean(group_ari, na.rm = TRUE),
      sd_group_ari = sd(group_ari, na.rm = TRUE),
      prop_correct_groups = mean(n_groups_found == n_groups_true, na.rm = TRUE),
      mean_iterations = mean(iterations, na.rm = TRUE),
      convergence_rate = mean(converged, na.rm = TRUE),
      .groups = 'drop'
    )
  
  print(summary_stats)
  
  # Return summary stats invisibly
  invisible(summary_stats)
}

#' Main execution function
main <- function() {
  
  # Check if required packages are installed
  required_packages <- c("ggplot2", "reshape2", "dplyr", "mclust")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages)
  }
  
  # Load required packages
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(mclust)
  
  # Define test configurations
  # Varying: dimensions, number of groups, sample sizes, noise proportions
  test_configs <- list(
    # Small dimensions, few groups
    list(
      n_groups = 2,
      n_per_group = 30,
      noise_prop = 0.1,
      dims = c(3, 3),
      noise_type = "uniform",
      adaptive_grid = TRUE,
      k_grid = NULL
    ),
    list(
      n_groups = 2,
      n_per_group = 30,
      noise_prop = 0.2,
      dims = c(3, 3),
      noise_type = "uniform",
      adaptive_grid = TRUE,
      k_grid = NULL
    ),
    # Medium dimensions
    list(
      n_groups = 3,
      n_per_group = 50,
      noise_prop = 0.15,
      dims = c(5, 5),
      noise_type = "uniform",
      adaptive_grid = TRUE,
      k_grid = NULL
    ),
    list(
      n_groups = 3,
      n_per_group = 50,
      noise_prop = 0.25,
      dims = c(5, 5),
      noise_type = "uniform",
      adaptive_grid = TRUE,
      k_grid = NULL
    ),
    # Compare with fixed grid (non-adaptive)
    list(
      n_groups = 3,
      n_per_group = 50,
      noise_prop = 0.15,
      dims = c(5, 5),
      noise_type = "uniform",
      adaptive_grid = FALSE,
      k_grid = 10^seq(-16, -1, length.out = 30)
    ),
    # Larger dimensions (computational)
    list(
      n_groups = 3,
      n_per_group = 40,
      noise_prop = 0.1,
      dims = c(8, 8),
      noise_type = "uniform",
      adaptive_grid = TRUE,
      k_grid = NULL
    ),
    # Higher noise proportion
    list(
      n_groups = 3,
      n_per_group = 60,
      noise_prop = 0.3,
      dims = c(5, 5),
      noise_type = "uniform",
      adaptive_grid = TRUE,
      k_grid = NULL
    )
  )
  
  # Run tests
  cat("Starting adaptive k selection tests...\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  results <- run_adaptive_k_tests(test_configs, n_replicates = 5)  # Reduced replicates for demo
  
  # Save results
  if (nrow(results) > 0) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    results_file <- paste0("adaptive_k_test_results_", timestamp, ".csv")
    write.csv(results, results_file, row.names = FALSE)
    cat(sprintf("\nResults saved to: %s\n", results_file))
    
    # Visualize
    cat("\nGenerating visualizations...\n")
    pdf(paste0("adaptive_k_visualization_", timestamp, ".pdf"), width = 12, height = 8)
    visualize_results(results)
    dev.off()
    cat("Visualizations saved to PDF\n")
    
    # Return results for further analysis
    invisible(results)
  } else {
    cat("\nNo results were generated. Please check the errors above.\n")
  }
}

# Run the tests
if (interactive()) {
  main()
} else {
  # For batch execution
  cat("Running in batch mode\n")
  results <- main()
}