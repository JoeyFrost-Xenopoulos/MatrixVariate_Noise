#!/usr/bin/env Rscript
#' Benchmark Comparison of Initialization Schemes
#'
#' Compares k-means, k-means++, random, and ECME initialization across multiple
#' simulated scenarios measuring:
#'   - Final log-likelihood (model quality)
#'   - Number of EM iterations to convergence
#'   - Clustering accuracy (Adjusted Rand Index)
#'   - Wall-clock time
#'
#' Usage:
#'   Rscript scripts/compare_initializations.R

# --- Source package code ---
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# --- Helpers ---

#' Adjusted Rand Index
#'
#' Computes ARI between two integer label vectors.
adjusted_rand_index <- function(labels_true, labels_pred) {
  n <- length(labels_true)
  if (n != length(labels_pred)) stop("Label vectors must have equal length.")

  # Contingency table
  tab <- table(labels_true, labels_pred)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2))
  c <- sum(choose(colSums(tab), 2))
  d <- choose(n, 2)

  expected <- (b * c) / d
  max_idx <- (b + c) / 2
  if (max_idx == expected) return(1)
  (a - expected) / (max_idx - expected)
}

#' Simulate matrix-variate data from a known mixture
#'
#' @param n_per_group Number of observations per group.
#' @param g Number of groups.
#' @param r Number of rows per matrix.
#' @param p Number of columns per matrix.
#' @param separation Mean separation between clusters.
#' @param noise_sd Noise standard deviation.
#' @return List with x_list and true_labels.
simulate_mixture_data <- function(n_per_group, g, r, p, separation, noise_sd = 0.5) {
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

#' Run a single benchmark trial
#'
#' @param x_list List of matrices.
#' @param g Number of components.
#' @param true_labels True cluster assignments.
#' @param init Initialization method.
#' @param max_iter Maximum EM iterations.
#' @param n_restarts Number of restarts (best result kept).
#' @return Data frame row with metrics.
run_trial <- function(x_list, g, true_labels, init, max_iter = 100, n_restarts = 5) {
  best_ll <- -Inf
  best_fit <- NULL
  total_time <- 0

  for (restart in seq_len(n_restarts)) {
    t0 <- proc.time()["elapsed"]
    fit <- tryCatch(
      suppressWarnings(
        matrix_variate_mixture_fit(
          x_list, g = g, max_iter = max_iter, init = init, verbose = FALSE
        )
      ),
      error = function(e) NULL
    )
    t1 <- proc.time()["elapsed"]
    total_time <- total_time + (t1 - t0)

    if (!is.null(fit)) {
      final_ll <- fit$logLik[length(fit$logLik)]
      if (final_ll > best_ll) {
        best_ll <- final_ll
        best_fit <- fit
      }
    }
  }

  if (is.null(best_fit)) {
    return(data.frame(
      init = init,
      loglik = NA_real_,
      iterations = NA_integer_,
      ari = NA_real_,
      time_sec = total_time,
      converged = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  ari <- tryCatch(
    adjusted_rand_index(true_labels, best_fit$cluster),
    error = function(e) NA_real_
  )

  data.frame(
    init = init,
    loglik = best_ll,
    iterations = best_fit$iterations,
    ari = ari,
    time_sec = total_time,
    converged = best_fit$converged,
    stringsAsFactors = FALSE
  )
}

# --- Benchmark Scenarios ---

scenarios <- list(
  list(name = "Easy (2x3, g=2, sep=4)", n = 20, g = 2, r = 2, p = 3, sep = 4, sd = 0.5),
  list(name = "Medium (3x4, g=3, sep=2.5)", n = 20, g = 3, r = 3, p = 4, sep = 2.5, sd = 0.7),
  list(name = "Hard (2x3, g=3, sep=1.5)", n = 25, g = 3, r = 2, p = 3, sep = 1.5, sd = 1.0),
  list(name = "Large matrices (4x5, g=2, sep=3)", n = 15, g = 2, r = 4, p = 5, sep = 3, sd = 0.6)
)

init_methods <- c("kmeans", "kmeans++", "random", "ecme")
n_trials <- 10  # trials per scenario x init combination

# --- Run Benchmark ---

cat(paste(rep("=", 72), collapse = ""), "\n")
cat("  INITIALIZATION SCHEME BENCHMARK\n")
cat("  Comparing: kmeans | kmeans++ | random | ecme\n")
cat("  Trials per condition:", n_trials, "\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

all_results <- data.frame()

for (scenario in scenarios) {
  cat(sprintf("\n--- Scenario: %s ---\n", scenario$name))
  cat(sprintf("    n_per_group=%d, g=%d, r=%d, p=%d, separation=%.1f, noise_sd=%.1f\n",
              scenario$n, scenario$g, scenario$r, scenario$p, scenario$sep, scenario$sd))

  for (trial in seq_len(n_trials)) {
    set.seed(1000 + trial)
    sim <- simulate_mixture_data(
      n_per_group = scenario$n,
      g = scenario$g,
      r = scenario$r,
      p = scenario$p,
      separation = scenario$sep,
      noise_sd = scenario$sd
    )

    for (init in init_methods) {
      result <- run_trial(sim$x_list, scenario$g, sim$true_labels, init)
      result$scenario <- scenario$name
      result$trial <- trial
      all_results <- rbind(all_results, result)
    }
  }

  # Print scenario summary
  scenario_data <- all_results[all_results$scenario == scenario$name, ]
  cat("\n    Summary (mean ± sd across", n_trials, "trials):\n")
  cat(sprintf("    %-8s  %12s  %10s  %8s  %10s\n",
              "Init", "Log-lik", "ARI", "Iters", "Time (s)"))
  cat("    ", paste(rep("-", 58), collapse = ""), "\n")

  for (init in init_methods) {
    d <- scenario_data[scenario_data$init == init, ]
    cat(sprintf("    %-8s  %6.1f ± %4.1f  %5.3f ± %4.3f  %4.1f ± %3.1f  %5.3f ± %4.3f\n",
                init,
                mean(d$loglik, na.rm = TRUE), sd(d$loglik, na.rm = TRUE),
                mean(d$ari, na.rm = TRUE), sd(d$ari, na.rm = TRUE),
                mean(d$iterations, na.rm = TRUE), sd(d$iterations, na.rm = TRUE),
                mean(d$time_sec, na.rm = TRUE), sd(d$time_sec, na.rm = TRUE)))
  }
}

# --- Final Aggregate Summary ---

cat("\n\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("  AGGREGATE RESULTS\n")
cat(paste(rep("=", 72), collapse = ""), "\n\n")

cat(sprintf("  %-8s  %12s  %10s  %10s  %12s  %10s\n",
            "Init", "Mean LL", "Mean ARI", "Mean Iter", "Mean Time", "Conv Rate"))
cat("  ", paste(rep("-", 68), collapse = ""), "\n")

for (init in init_methods) {
  d <- all_results[all_results$init == init, ]
  conv_rate <- mean(d$converged, na.rm = TRUE)
  cat(sprintf("  %-8s  %12.1f  %10.3f  %10.1f  %10.3f s  %9.1f%%\n",
              init,
              mean(d$loglik, na.rm = TRUE),
              mean(d$ari, na.rm = TRUE),
              mean(d$iterations, na.rm = TRUE),
              mean(d$time_sec, na.rm = TRUE),
              conv_rate * 100))
}

cat("\n")
cat("  Best init by ARI:       ", init_methods[which.max(
  sapply(init_methods, function(m) mean(all_results$ari[all_results$init == m], na.rm = TRUE))
)], "\n")
cat("  Best init by log-lik:   ", init_methods[which.max(
  sapply(init_methods, function(m) mean(all_results$loglik[all_results$init == m], na.rm = TRUE))
)], "\n")
cat("  Fastest init:           ", init_methods[which.min(
  sapply(init_methods, function(m) mean(all_results$time_sec[all_results$init == m], na.rm = TRUE))
)], "\n")

# --- Save results to CSV ---
output_file <- "scripts/init_benchmark_results.csv"
write.csv(all_results, output_file, row.names = FALSE)
cat(sprintf("\n  Results saved to: %s\n", output_file))
