################################
# HC Noise Grid Search Quality Benchmark
#
# Goal:
#   Evaluate whether heuristic k-grid is sufficient to
#   capture optimal noise constant across scenarios
################################

rm(list = ls())

library(Ampharos)
library(mclust)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(tibble)

set.seed(123)

# Grid definition

k_grid <- 10^seq(-16, -1, length.out = 30)

# Helper: compute KS score for a given k

compute_k_score <- function(fit, x_list) {
  matrix_noise_ks_score(fit, x_list)$statistic
}

# Run grid search for a single dataset

evaluate_k_grid <- function(x_list, g, init_method) {

  scores <- numeric(length(k_grid))

  for (i in seq_along(k_grid)) {

    fit <- tryCatch(
      matrix_variate_noise_fit(
        x_list = x_list,
        g = g,
        noise_type = "hc",
        init = init_method,
        estimate_k = FALSE,
        noise_k = k_grid[i],
        verbose = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      scores[i] <- NA
    } else {
      scores[i] <- compute_k_score(fit, x_list)
    }
  }

  tibble(
    k = k_grid,
    log_k = log10(k_grid),
    score = scores
  )
}

# Extract optimal k from grid

extract_best_k <- function(grid_df) {

  grid_df <- grid_df %>% filter(!is.na(score))

  if (nrow(grid_df) == 0) return(NA_real_)

  grid_df$k[which.max(grid_df$score)]
}

# Grid diagnostics

grid_diagnostics <- function(grid_df, k_star) {

  k_min <- min(grid_df$k, na.rm = TRUE)
  k_max <- max(grid_df$k, na.rm = TRUE)

  tibble(
    k_star = k_star,
    hit_grid = k_star >= k_min & k_star <= k_max,
    at_lower_edge = abs(k_star - k_min) < 1e-12,
    at_upper_edge = abs(k_star - k_max) < 1e-12,
    log_error = min(abs(log10(grid_df$k) - log10(k_star)))
  )
}

# Single experiment wrapper

run_grid_experiment <- function(x_list, g, init_method) {

  grid_df <- evaluate_k_grid(x_list, g, init_method)

  k_star <- extract_best_k(grid_df)

  diag <- grid_diagnostics(grid_df, k_star)

  list(
    grid = grid_df,
    diagnostics = diag
  )
}

# Example simulation (replace with full benchmark loop later)

simulate_simple <- function(n = 100, r = 3, p = 3, g = 2) {

  M1 <- matrix(0, r, p)
  M2 <- matrix(2, r, p)

  x_list <- c(
    replicate(n/2, M1 + matrix(rnorm(r*p), r, p), simplify = FALSE),
    replicate(n/2, M2 + matrix(rnorm(r*p), r, p), simplify = FALSE)
  )

  x_list
}

# Run test

x_list <- simulate_simple()

res_kmeans <- run_grid_experiment(x_list, g = 2, init_method = "kmeans")
res_ecme   <- run_grid_experiment(x_list, g = 2, init_method = "ecme")

# Combine outputs

summary <- bind_rows(
  mutate(res_kmeans$diagnostics, init = "kmeans"),
  mutate(res_ecme$diagnostics, init = "ecme")
)

print(summary)

write_csv(summary, "grid_quality_summary.csv")