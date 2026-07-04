
# Benchmark: Automatic HC Noise Constant Estimation
#
# Compares automatic estimate_k performance using
#   • K-means initialization
#   • ECME initialization
#
# Outputs
# --------
# estimate_k_summary.csv
# estimate_k_candidates.csv


rm(list = ls())

# library(Ampharos)

library(mclust)

library(clusterGeneration)

library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)
library(stringr)

set.seed(12345)

# Output directory

output_dir <- "benchmark_results"

if (!dir.exists(output_dir))
    dir.create(output_dir)

# Benchmark Parameters

benchmark_settings <- list(
    
    # Monte Carlo
    
    n_replications = 100,

    # Matrix dimensions
    
    dimensions = list(

        c(2,2),
        c(2,4),
        c(3,3),
        c(3,5),
        c(4,4),
        c(5,5),
        c(6,6)

    ),

    # Number of Gaussian components

    n_groups = c(2,3),
    
    # Sample size  

    sample_size = c(
        150,
        300
    ),
    
    # Initializations
    
    initialization = c(

        "kmeans",
        "ecme"

    ),
    
    # Noise mechanisms
    
    noise_mechanism = c(

        "matrix",
        "column",
        "element",
        "permutation"

    ),

    # Noise intensity

    contamination = c(

        0.00,
        0.05,
        0.10,
        0.20,
        0.30,
        0.40

    ),

    # Mixture fitting

    max_iter = 100,

    tol = 1e-6,

    nstart = 100

)

# Candidate k Grid

candidate_k_grid <- 10^seq(

    -16,
    -1,

    length.out = 30

)

# Helper function

simulation_id <- function(replicate,
                          dimension,
                          groups,
                          contamination,
                          mechanism,
                          initialization){

    paste(

        paste0(
            dimension[1],
            "x",
            dimension[2]
        ),

        paste0(
            "G",
            groups
        ),

        mechanism,

        paste0(
            contamination*100,
            "pct"
        ),

        initialization,

        paste0(
            "Rep",
            replicate
        ),

        sep = "_"

    )

}

# Matrix-variate Gaussian generator (core engine)

generate_matrix_observation <- function(mean_matrix,
                                        row_cov,
                                        col_cov,
                                        r,
                                        p) {

  Z <- matrix(rnorm(r * p), r, p)

  row_chol <- chol(make_spd(row_cov))
  col_chol <- chol(make_spd(col_cov))

  mean_matrix +
    row_chol %*% Z %*% t(col_chol)
}

# Clean mixture generator

simulate_matrix_mixture <- function(n,
                                     means,
                                     row_covs,
                                     col_covs,
                                     proportions) {

  g <- length(means)
  r <- nrow(means[[1]])
  p <- ncol(means[[1]])

  cluster_true <- sample(
    seq_len(g),
    size = n,
    replace = TRUE,
    prob = proportions
  )

  x_list <- vector("list", n)

  for (i in seq_len(n)) {
    k <- cluster_true[i]

    x_list[[i]] <- generate_matrix_observation(
      mean_matrix = means[[k]],
      row_cov = row_covs[[k]],
      col_cov = col_covs[[k]],
      r = r,
      p = p
    )
  }

  list(
    x_list = x_list,
    cluster_true = cluster_true
  )
}

# Matrix-level contamination

contam_matrix_level <- function(x_list,
                                proportion = 0.1,
                                r,
                                p,
                                low = -5,
                                high = 5) {

  n <- length(x_list)
  n_contam <- ceiling(proportion * n)

  idx <- sample(seq_len(n), n_contam)

  for (i in idx) {
    x_list[[i]] <- matrix(
      runif(r * p, low, high),
      r,
      p
    )
  }

  list(
    x_list = x_list,
    contam_idx = idx
  )
}

# Column contamination

contam_column <- function(x_list,
                          proportion = 0.1,
                          low = -10,
                          high = 10) {

  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  n_contam <- ceiling(proportion * n)
  idx <- sample(seq_len(n), n_contam)

  for (i in idx) {

    col_id <- sample(seq_len(p), 1)

    x_list[[i]][, col_id] <- runif(r, low, high)
  }

  list(
    x_list = x_list,
    contam_idx = idx
  )
}

# Element-wise contamination

contam_elementwise <- function(x_list,
                               proportion = 0.1,
                               entry_prob = 0.1,
                               low = -10,
                               high = 10) {

  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  n_contam <- ceiling(proportion * n)
  idx <- sample(seq_len(n), n_contam)

  for (i in idx) {

    mask <- matrix(
      rbinom(r * p, 1, entry_prob),
      r,
      p
    )

    noise <- matrix(runif(r * p, low, high), r, p)

    x_list[[i]] <- x_list[[i]] * (1 - mask) + noise * mask
  }

  list(
    x_list = x_list,
    contam_idx = idx
  )
}

# Permutation contamination

contam_permutation <- function(x_list,
                               proportion = 0.1) {

  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  n_contam <- ceiling(proportion * n)
  idx <- sample(seq_len(n), n_contam)

  for (i in idx) {
    x_list[[i]] <- matrix(
      sample(x_list[[i]]),
      r,
      p
    )
  }

  list(
    x_list = x_list,
    contam_idx = idx
  )
}

# Contamination dispatcher

apply_contamination <- function(x_list,
                                type = c("matrix",
                                         "column",
                                         "element",
                                         "permutation"),
                                proportion,
                                r,
                                p) {

  type <- match.arg(type)

  if (type == "matrix") {
    return(contam_matrix_level(x_list, proportion, r, p))

  } else if (type == "column") {
    return(contam_column(x_list, proportion))

  } else if (type == "element") {
    return(contam_elementwise(x_list, proportion))

  } else if (type == "permutation") {
    return(contam_permutation(x_list, proportion))
  }
}

################################
# Evaluation metrics
################################

compute_misclassification_rate <- function(true, pred) {
  mean(true != pred)
}

compute_ari <- function(true, pred) {
  mclust::adjustedRandIndex(true, pred)
}

################################
# Safe HC fitting wrapper
################################

run_hc_fit <- function(x_list,
                        g,
                        init_method,
                        noise_type = "hc",
                        estimate_k = TRUE,
                        noise_k = NULL,
                        max_iter = 100,
                        tol = 1e-6,
                        nstart = 100,
                        verbose = FALSE) {

  fit <- tryCatch(
    matrix_variate_noise_fit(
      x_list = x_list,
      g = g,
      noise_type = noise_type,
      init = init_method,
      estimate_k = estimate_k,
      noise_k = noise_k,
      max_iter = max_iter,
      tol = tol,
      nstart = nstart,
      verbose = verbose
    ),
    error = function(e) {
      return(NULL)
    }
  )

  fit
}

################################
# Extract estimated k
################################

extract_k <- function(fit) {

  if (is.null(fit)) return(NA_real_)

  if (!is.null(fit$k_selection)) {
    return(fit$k_selection$selected_k)
  }

  if (!is.null(fit$noise$k)) {
    return(fit$noise$k)
  }

  NA_real_
}

################################
# Single benchmark run
################################

run_benchmark_case <- function(replicate_id,
                               dimension,
                               g,
                               n,
                               init_method,
                               contamination_type,
                               contamination_level,
                               candidate_k_grid,
                               means,
                               row_covs,
                               col_covs,
                               proportions) {

  r <- dimension[1]
  p <- dimension[2]

  ##############################
  # 1. Generate clean data
  ##############################

  sim <- simulate_matrix_mixture(
    n = n,
    means = means,
    row_covs = row_covs,
    col_covs = col_covs,
    proportions = proportions
  )

  x_list <- sim$x_list
  true_cluster <- sim$cluster_true

  ##############################
  # 2. Apply contamination
  ##############################

  contam <- apply_contamination(
    x_list = x_list,
    type = contamination_type,
    proportion = contamination_level,
    r = r,
    p = p
  )

  x_list <- contam$x_list
  contam_idx <- contam$contam_idx

  ##############################
  # 3. Fit HC model (with k estimation)
  ##############################

  fit <- run_hc_fit(
    x_list = x_list,
    g = g,
    init_method = init_method,
    estimate_k = TRUE,
    verbose = FALSE
  )

  ##############################
  # 4. Extract results
  ##############################

  pred_cluster <- if (!is.null(fit)) fit$cluster else rep(NA, n)

  ari <- compute_ari(true_cluster, pred_cluster)
  misclass <- compute_misclassification_rate(true_cluster, pred_cluster)
  est_k <- extract_k(fit)

  ##############################
  # 5. Confusion table
  ##############################

  cluster_table <- if (!is.null(fit)) {
    table(Predicted = pred_cluster, True = true_cluster)
  } else {
    matrix(NA)
  }

  ##############################
  # 6. Summary output
  ##############################

  summary_row <- tibble::tibble(
    replicate = replicate_id,
    dimension_r = r,
    dimension_p = p,
    groups = g,
    sample_size = n,
    init = init_method,
    contamination_type = contamination_type,
    contamination_level = contamination_level,
    estimated_k = est_k,
    ARI = ari,
    misclassification = misclass,
    n_noise = length(contam_idx)
  )

  ##############################
  # 7. Full candidate log (placeholder structure)
  ##############################

  candidate_log <- tibble::tibble(
    replicate = replicate_id,
    init = init_method,
    contamination_type = contamination_type,
    contamination_level = contamination_level,
    estimated_k = est_k,
    cluster_table = paste(capture.output(print(cluster_table)), collapse = "\n")
  )

  list(
    summary = summary_row,
    candidates = candidate_log
  )
}

# Experimental design grid

design_grid <- expand.grid(
  dimension = seq_along(benchmark_settings$dimensions),
  groups = benchmark_settings$n_groups,
  n = benchmark_settings$sample_size,
  init = benchmark_settings$initialization,
  contamination_type = benchmark_settings$noise_mechanism,
  contamination_level = benchmark_settings$contamination,
  replicate = seq_len(benchmark_settings$n_replications),
  KEEP.OUT.ATTRS = FALSE
)

design_grid$dimension_value <- benchmark_settings$dimensions[design_grid$dimension]

# Output storage

summary_results <- list()
candidate_results <- list()

# Scenario-specific mixture generator

generate_scenario_mixture <- function(r, p, g) {

  means <- vector("list", g)
  row_covs <- vector("list", g)
  col_covs <- vector("list", g)

  for (k in seq_len(g)) {

    # structured but separable means
    M <- matrix(rnorm(r * p, sd = 1), r, p)
    means[[k]] <- M + k * 0.5

    row_covs[[k]] <- diag(runif(r, 0.5, 1.5))
    col_covs[[k]] <- diag(runif(p, 0.5, 1.5))
  }

  proportions <- rep(1 / g, g)

  list(
    means = means,
    row_covs = row_covs,
    col_covs = col_covs,
    proportions = proportions
  )
}

# Main benchmark execution loop

for (i in seq_len(nrow(design_grid))) {

  config <- design_grid[i, ]

  r <- config$dimension_value[[1]]
  p <- config$dimension_value[[2]]

  g <- config$groups
  n <- config$n
  init_method <- config$init
  contamination_type <- config$contamination_type
  contamination_level <- config$contamination_level
  rep <- config$replicate

  # Create simulation ID

  sim_id <- simulation_id(
    replicate = rep,
    dimension = c(r, p),
    groups = g,
    contamination = contamination_level,
    mechanism = contamination_type,
    initialization = init_method
  )

  cat("\nRunning:", sim_id, "\n")

  # Generate mixture model

  mix <- generate_scenario_mixture(r, p, g)

  # Run benchmark

  res <- run_benchmark_case(
    replicate_id = sim_id,
    dimension = c(r, p),
    g = g,
    n = n,
    init_method = init_method,
    contamination_type = contamination_type,
    contamination_level = contamination_level,
    candidate_k_grid = candidate_k_grid,
    means = mix$means,
    row_covs = mix$row_covs,
    col_covs = mix$col_covs,
    proportions = mix$proportions
  )

  # Store results

  summary_results[[i]] <- res$summary
  candidate_results[[i]] <- res$candidates

  if (i %% 50 == 0) {

    summary_df <- dplyr::bind_rows(summary_results)
    candidate_df <- dplyr::bind_rows(candidate_results)

    readr::write_csv(summary_df,
                     file.path(output_dir, "estimate_k_summary_partial.csv"))

    readr::write_csv(candidate_df,
                     file.path(output_dir, "estimate_k_candidates_partial.csv"))

    cat("\nCheckpoint saved at iteration", i, "\n")
  }
}

# Final output export

summary_df <- dplyr::bind_rows(summary_results)
candidate_df <- dplyr::bind_rows(candidate_results)

readr::write_csv(summary_df,
                 file.path(output_dir, "estimate_k_summary.csv"))

readr::write_csv(candidate_df,
                 file.path(output_dir, "estimate_k_candidates.csv"))
