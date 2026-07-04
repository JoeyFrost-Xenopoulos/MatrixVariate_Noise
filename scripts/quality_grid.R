###############################################################
## Grid Search Heuristic Validation (CLEAN VERSION)
##
## Tests whether candidate_k_grid recovers optimal noise k
## using explicit profile evaluation over noise_k
###############################################################

rm(list = ls())

library(Ampharos)
library(mclust)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)

set.seed(12345)

###############################################################
## Output directory
###############################################################

output_dir <- "grid_search_diagnostics"
if (!dir.exists(output_dir)) dir.create(output_dir)

###############################################################
## Candidate grid
###############################################################

candidate_k_grid <- 10^seq(-16, -1, length.out = 30)

###############################################################
## Simulation settings
###############################################################

test_settings <- list(
  n_replications = 30,
  dimensions = list(c(3,3), c(5,5)),
  n_groups = c(2,3),
  sample_size = c(150),
  contamination = c(0, 0.1, 0.3),
  noise_type = c("matrix", "column"),
  init = c("kmeans", "ecme")
)

###############################################################
## Storage
###############################################################

results_list <- list()

###############################################################
## Scenario generator (assumes already defined elsewhere)
###############################################################

# generate_scenario_mixture()
# simulate_matrix_mixture()
# apply_contamination()

###############################################################
## Fit model at fixed k
###############################################################

fit_fixed_k <- function(x_list, g, init_method, k_val) {

  tryCatch(
    matrix_variate_noise_fit(
      x_list = x_list,
      g = g,
      init = init_method,
      estimate_k = FALSE,
      noise_k = k_val,
      max_iter = 100,
      tol = 1e-6,
      nstart = 100,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
}

###############################################################
## Profile oracle over k grid
###############################################################

profile_k_grid <- function(x_list, g, init_method, grid) {

  scores <- numeric(length(grid))

  for (i in seq_along(grid)) {

    fit <- fit_fixed_k(x_list, g, init_method, grid[i])

    if (is.null(fit)) {
      scores[i] <- NA_real_
    } else if (!is.null(fit$loglik)) {
      scores[i] <- fit$loglik
    } else {
      scores[i] <- NA_real_
    }
  }

  best_idx <- which.max(scores)

  list(
    grid = grid,
    scores = scores,
    oracle_k = grid[best_idx],
    oracle_score = scores[best_idx]
  )
}

###############################################################
## Run single experiment
###############################################################

run_case <- function(dim, g, n, cont, type, init) {

  r <- dim[1]
  p <- dim[2]

  mix <- generate_scenario_mixture(r, p, g)

  sim <- simulate_matrix_mixture(
    n = n,
    means = mix$means,
    row_covs = mix$row_covs,
    col_covs = mix$col_covs,
    proportions = mix$proportions
  )

  contam <- apply_contamination(
    sim$x_list,
    type = type,
    proportion = cont,
    r = r,
    p = p
  )

  x_list <- contam$x_list

  ###########################################################
  ## Oracle via grid profiling
  ###########################################################

  oracle <- profile_k_grid(
    x_list = x_list,
    g = g,
    init_method = init,
    grid = candidate_k_grid
  )

  ###########################################################
  ## Automatic estimator
  ###########################################################

  fit_auto <- tryCatch(
    matrix_variate_noise_fit(
      x_list = x_list,
      g = g,
      init = init,
      estimate_k = TRUE,
      max_iter = 100,
      tol = 1e-6,
      nstart = 100,
      verbose = FALSE
    ),
    error = function(e) NULL
  )

  selected_k <- if (!is.null(fit_auto$k_selection)) {
    fit_auto$k_selection$selected_k
  } else {
    NA_real_
  }

  ###########################################################
  ## Metrics
  ###########################################################

  log_error <- abs(log10(selected_k) - log10(oracle$oracle_k))

  tibble(
    r = r,
    p = p,
    g = g,
    n = n,
    contamination = cont,
    type = type,
    init = init,
    oracle_k = oracle$oracle_k,
    selected_k = selected_k,
    log_error = log_error
  )
}

###############################################################
## Main loop
###############################################################

for (d in seq_along(test_settings$dimensions)) {
  for (g in test_settings$n_groups) {
    for (n in test_settings$sample_size) {
      for (cont in test_settings$contamination) {
        for (type in test_settings$noise_type) {
          for (init in test_settings$init) {

            dim <- test_settings$dimensions[[d]]

            cat("\nRunning:", dim[1], "x", dim[2],
                "G", g,
                "cont", cont,
                "type", type,
                "init", init, "\n")

            res <- run_case(dim, g, n, cont, type, init)

            results_list[[length(results_list) + 1]] <- res
          }
        }
      }
    }
  }
}

###############################################################
## Save outputs
###############################################################

results_df <- bind_rows(results_list)

write_csv(results_df,
          file.path(output_dir, "grid_search_diagnostics.csv"))

###############################################################
## Summary table
###############################################################

summary_df <- results_df %>%
  group_by(r, p, g, contamination, type, init) %>%
  summarise(
    mean_log_error = mean(log_error, na.rm = TRUE),
    sd_log_error = sd(log_error, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_df,
          file.path(output_dir, "grid_search_summary.csv"))