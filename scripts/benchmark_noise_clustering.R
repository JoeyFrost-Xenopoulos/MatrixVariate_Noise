#!/usr/bin/env Rscript
#' Benchmark: Noise Mixture Clustering
#'
#' Tests matrix_variate_noise_fit across combinations of initialization
#' schemes and KS score types. Records clustering accuracy, noise recovery,
#' convergence, and automatic k-selection diagnostics.
#'
#' Usage:
#'   Rscript scripts/benchmark_noise_clustering.R
#'
#' Adjust the parameters below before running.

# --- Source package code ---
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# ============================================================
# Configuration
# ============================================================

# Simulation settings
N_TRIALS        <- 10   # repetitions per scenario
N_RESTARTS      <- 5    # restarts per fit (best log-lik kept)

# Methods to sweep
INIT_METHODS    <- c("kmeans", "kmeans++", "random", "ecme")
KS_TYPES        <- c("onesample", "twosample")
NOISE_TYPES     <- c("hc", "br")

# Fixed fit parameters
MAX_ITER        <- 200
TOL             <- 1e-6

# Output
OUTPUT_CSV      <- "scripts/noise_clustering_results.csv"

# ============================================================
# Helpers
# ============================================================

#' Adjusted Rand Index
adjusted_rand_index <- function(labels_true, labels_pred) {
  n <- length(labels_true)
  if (n != length(labels_pred)) stop("Label vectors must have equal length.")
  tab <- table(labels_true, labels_pred)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2))
  c <- sum(choose(colSums(tab), 2))
  d <- choose(n, 2)
  if (d == 0) return(1)
  expected <- (b * c) / d
  max_idx <- (b + c) / 2
  if (max_idx == expected) return(1)
  (a - expected) / (max_idx - expected)
}

#' Simulate a matrix-variate mixture with known parameters
simulate_mixture <- function(n_per_group, g, r, p, separation, noise_sd = 0.5) {
  true_means <- lapply(seq_len(g), function(k) {
    matrix(rnorm(r * p, mean = separation * (k - (g + 1) / 2)), r, p)
  })
  x_list <- list()
  true_labels <- integer(0)
  for (k in seq_len(g)) {
    for (i in seq_len(n_per_group)) {
      noise <- matrix(rnorm(r * p, sd = noise_sd), r, p)
      x_list <- c(x_list, list(true_means[[k]] + noise))
      true_labels <- c(true_labels, k)
    }
  }
  list(x_list = x_list, true_labels = true_labels)
}

#' Add contiguous outliers (permuted entries) to simulate data-entry errors
add_permuted_outliers <- function(x_list, n_outliers = 15) {
  n <- length(x_list)
  out_idx <- sample(n, min(n_outliers, n))
  for (idx in out_idx) {
    x_list[[idx]] <- matrix(sample(x_list[[idx]]), nrow = nrow(x_list[[idx]]))
  }
  x_list
}

#' Add column-replacement outliers (Tomarchio-style)
add_column_outliers <- function(x_list, n_outliers = 10, p, r, range = 15) {
  n <- length(x_list)
  out_idx <- sample(n, min(n_outliers, n))
  for (idx in out_idx) {
    col_to_replace <- sample(p, 1)
    x_list[[idx]][, col_to_replace] <- runif(r, -range, range)
  }
  x_list
}

#' Run one fit and return key metrics as a data.frame row
run_fit <- function(x_list, g, true_labels, init, noise_type,
                    ks_type, noise_k = 1e-4, estimate_k = FALSE) {
  fit <- tryCatch(
    suppressWarnings(
      matrix_variate_noise_fit(
        x_list, g = g, noise_type = noise_type,
        init = init, ks_type = ks_type,
        noise_k = noise_k, estimate_k = estimate_k,
        max_iter = MAX_ITER, tol = TOL,
        nstart = N_RESTARTS, verbose = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.frame(
      scenario = NA, init = init, ks_type = ks_type,
      noise_type = noise_type, estimate_k = estimate_k,
      trial = NA, converged = FALSE, iterations = NA,
      loglik = NA_real_, ari = NA_real_,
      noise_pi = NA_real_, noise_detected = NA,
      n_noise = NA, selected_k = NA,
      grid_min = NA, grid_max = NA, grid_len = NA,
      ks_selected = NA,
      stringsAsFactors = FALSE
    ))
  }

  keep <- fit$cluster > 0
  ari <- tryCatch(
    adjusted_rand_index(true_labels[keep], fit$cluster[keep]),
    error = function(e) NA_real_
  )

  n_noise <- sum(fit$cluster == 0)

  ks_info <- if (!is.null(fit$k_selection)) {
    list(
      selected_k = fit$k_selection$selected_k,
      grid_min = if (!is.null(fit$k_selection$k_grid))
        min(fit$k_selection$k_grid) else NA,
      grid_max = if (!is.null(fit$k_selection$k_grid))
        max(fit$k_selection$k_grid) else NA,
      grid_len = length(fit$k_selection$k_grid),
      ks_selected = min(fit$k_selection$ks_scores, na.rm = TRUE)
    )
  } else {
    list(selected_k = NA, grid_min = NA, grid_max = NA,
         grid_len = NA, ks_selected = NA)
  }

  data.frame(
    scenario = NA, init = init, ks_type = ks_type,
    noise_type = noise_type, estimate_k = estimate_k,
    trial = NA, converged = fit$converged,
    iterations = fit$iterations,
    loglik = fit$logLik[length(fit$logLik)],
    ari = ari,
    noise_pi = fit$noise$pi,
    noise_detected = n_noise > 0,
    n_noise = n_noise,
    selected_k = ks_info$selected_k,
    grid_min = ks_info$grid_min,
    grid_max = ks_info$grid_max,
    grid_len = ks_info$grid_len,
    ks_selected = ks_info$ks_selected,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Scenarios
# ============================================================

scenarios <- list(
  list(
    name = "Easy (2x3, g=2, sep=4)",
    n = 20, g = 2, r = 2, p = 3, sep = 4, sd = 0.5,
    sim_type = "basic"
  ),
  list(
    name = "Medium (4x4, g=3, sep=2.5)",
    n = 25, g = 3, r = 4, p = 4, sep = 2.5, sd = 0.7,
    sim_type = "basic"
  ),
  list(
    name = "Hard (2x3, g=3, sep=1.5)",
    n = 25, g = 3, r = 2, p = 3, sep = 1.5, sd = 1.0,
    sim_type = "basic"
  ),
  list(
    name = "Permuted outliers (3x5, g=3, sep=3)",
    n = 20, g = 3, r = 3, p = 5, sep = 3, sd = 0.6,
    sim_type = "permuted_outliers", n_outliers = 15
  ),
  list(
    name = "Column outliers (2x4, g=2, sep=2.5)",
    n = 20, g = 2, r = 2, p = 4, sep = 2.5, sd = 0.6,
    sim_type = "column_outliers", n_outliers = 10
  )
)

# ============================================================
# Run benchmark
# ============================================================

cat(paste(rep("=", 72), collapse = ""), "\n")
cat("  NOISE MIXTURE CLUSTERING BENCHMARK\n")
cat("  Init methods: ", paste(INIT_METHODS, collapse = ", "), "\n", sep = "")
cat("  KS types:     ", paste(KS_TYPES, collapse = ", "), "\n", sep = "")
cat("  Noise types:  ", paste(NOISE_TYPES, collapse = ", "), "\n", sep = "")
cat("  Trials:       ", N_TRIALS, "\n", sep = "")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

all_results <- data.frame()

for (scenario in scenarios) {
  cat(sprintf("\n--- Scenario: %s ---\n", scenario$name))

  for (trial in seq_len(N_TRIALS)) {
    set.seed(10000 + trial)
    sim <- simulate_mixture(
      n_per_group = scenario$n, g = scenario$g,
      r = scenario$r, p = scenario$p,
      separation = scenario$sep, noise_sd = scenario$sd
    )

    if (scenario$sim_type == "permuted_outliers") {
      sim$x_list <- add_permuted_outliers(sim$x_list, scenario$n_outliers)
    } else if (scenario$sim_type == "column_outliers") {
      sim$x_list <- add_column_outliers(
        sim$x_list, scenario$n_outliers,
        p = scenario$p, r = scenario$r
      )
    }

    true_noise_prop <- 0
    if (scenario$sim_type %in% c("permuted_outliers", "column_outliers")) {
      true_noise_prop <- scenario$n_outliers / length(sim$x_list)
    }

    total_combos <- length(NOISE_TYPES) * length(INIT_METHODS) * length(KS_TYPES)
    combo_idx <- 1

    for (noise_type in NOISE_TYPES) {
      for (init in INIT_METHODS) {
        for (ks_type in KS_TYPES) {
          cat(sprintf("  Trial %d/%d: %s | %s | %s | %s\r",
                      trial, N_TRIALS, scenario$name,
                      noise_type, init, ks_type))

          result <- run_fit(
            x_list = sim$x_list, g = scenario$g,
            true_labels = sim$true_labels,
            init = init, noise_type = noise_type,
            ks_type = ks_type,
            noise_k = 1e-4,
            estimate_k = TRUE
          )
          result$scenario <- scenario$name
          result$trial <- trial
          result$true_noise_prop <- true_noise_prop
          result$n <- length(sim$x_list)
          result$true_g <- scenario$g
          result$r <- scenario$r
          result$p <- scenario$p
          all_results <- rbind(all_results, result)
          combo_idx <- combo_idx + 1
        }
      }
    }
  }
}

# ============================================================
# Summary tables
# ============================================================

cat("\n\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("  SCENARIO SUMMARY: ARI (mean +/- sd)\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

for (scenario in scenarios) {
  cat(sprintf("\n  %s\n", scenario$name))
  sdata <- all_results[all_results$scenario == scenario$name, ]
  header <- sprintf(
    "    %-10s | %-12s | %-8s | %10s | %10s | %10s | %10s",
    "Noise", "Init", "KS", "Conv %", "ARI", "NoisePi", "Iters"
  )
  cat(header, "\n")
  cat("    ", paste(rep("-", nchar(header) - 2), collapse = ""), "\n", sep = "")

  for (noise_type in NOISE_TYPES) {
    for (init in INIT_METHODS) {
      for (ks_type in KS_TYPES) {
        sub <- sdata[sdata$noise_type == noise_type &
                     sdata$init == init &
                     sdata$ks_type == ks_type, ]
        if (nrow(sub) == 0) next
        conv_rate <- mean(sub$converged, na.rm = TRUE) * 100
        ari_m <- mean(sub$ari, na.rm = TRUE)
        ari_s <- sd(sub$ari, na.rm = TRUE)
        noise_pi_m <- mean(sub$noise_pi, na.rm = TRUE)
        noise_pi_s <- sd(sub$noise_pi, na.rm = TRUE)
        iters_m <- mean(sub$iterations, na.rm = TRUE)
        iters_s <- sd(sub$iterations, na.rm = TRUE)
        cat(sprintf(
          "    %-10s | %-12s | %-8s | %8.0f%% | %5.3f+-%.3f | %5.3f+-%.3f | %5.1f+-%.1f\n",
          noise_type, init, ks_type,
          conv_rate, ari_m, ari_s, noise_pi_m, noise_pi_s,
          iters_m, iters_s
        ))
      }
    }
  }
}

# ============================================================
# Aggregate ranking
# ============================================================

cat("\n\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("  AGGREGATE RANKING (by mean ARI)\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

cat(sprintf(
  "  %-10s | %-8s | %-12s | %8s | %10s | %10s | %10s\n",
  "Noise", "KS", "Init", "Conv %", "Mean ARI", "Mean NoisePi", "Mean Iters"
))
cat("  ", paste(rep("-", 68), collapse = ""), "\n", sep = "")

for (noise_type in NOISE_TYPES) {
  for (ks_type in KS_TYPES) {
    for (init in INIT_METHODS) {
      sub <- all_results[all_results$noise_type == noise_type &
                          all_results$ks_type == ks_type &
                          all_results$init == init, ]
      if (nrow(sub) == 0) next
      conv_rate <- mean(sub$converged, na.rm = TRUE) * 100
      ari_m <- mean(sub$ari, na.rm = TRUE)
      noise_pi_m <- mean(sub$noise_pi, na.rm = TRUE)
      iters_m <- mean(sub$iterations, na.rm = TRUE)
      cat(sprintf(
        "  %-10s | %-8s | %-12s | %7.0f%% | %10.3f | %10.3f | %10.1f\n",
        noise_type, ks_type, init, conv_rate, ari_m, noise_pi_m, iters_m
      ))
    }
  }
}

# ============================================================
# KS selection diagnostics
# ============================================================

cat("\n\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("  KS SELECTION DIAGNOSTICS (estimate_k = TRUE)\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

ks_data <- all_results[all_results$estimate_k == TRUE &
                        !is.na(all_results$selected_k), ]
if (nrow(ks_data) > 0) {
  cat(sprintf(
    "  %-10s | %-8s | %-12s | %8s | %12s | %12s\n",
    "Noise", "KS", "Init", "Conv %",
    "Mean Sel k", "Mean KS stat"
  ))
  cat("  ", paste(rep("-", 68), collapse = ""), "\n", sep = "")

  for (noise_type in NOISE_TYPES) {
    for (ks_type in KS_TYPES) {
      for (init in INIT_METHODS) {
        sub <- ks_data[ks_data$noise_type == noise_type &
                        ks_data$ks_type == ks_type &
                        ks_data$init == init, ]
        if (nrow(sub) == 0) next
        conv_rate <- mean(sub$converged, na.rm = TRUE) * 100
        sel_k_m <- mean(sub$selected_k, na.rm = TRUE)
        ks_m <- mean(sub$ks_selected, na.rm = TRUE)
        cat(sprintf(
          "  %-10s | %-8s | %-12s | %7.0f%% | %12.4e | %12.4f\n",
          noise_type, ks_type, init, conv_rate, sel_k_m, ks_m
        ))
      }
    }
  }
} else {
  cat("  No KS selection results available.\n")
}

# ============================================================
# Save
# ============================================================

write.csv(all_results, OUTPUT_CSV, row.names = FALSE)
cat("\n\n  Results saved to:", OUTPUT_CSV, "\n")
