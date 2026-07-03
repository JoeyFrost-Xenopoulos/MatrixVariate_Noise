#!/usr/bin/env Rscript

# --- Source package code ---
for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# ============================================================
# Configuration
# ============================================================

N_TRIALS   <- 10
N_RESTARTS <- 5

INIT_METHODS <- c("kmeans", "kmeans++", "random", "ecme")

NOISE_TYPE <- "hc"

MAX_ITER <- 200
TOL      <- 1e-6

OUTPUT_CSV <- "scripts/noise_clustering_hc_estk_results.csv"

# ============================================================
# Metrics
# ============================================================

adjusted_rand_index <- function(labels_true, labels_pred) {
  n <- length(labels_true)
  tab <- table(labels_true, labels_pred)

  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2))
  c <- sum(choose(colSums(tab), 2))
  d <- choose(n, 2)

  if (d == 0) return(1)

  expected <- (b * c) / d
  max_idx  <- (b + c) / 2

  if (max_idx == expected) return(1)
  (a - expected) / (max_idx - expected)
}

# ============================================================
# Simulation WITH TRUE NOISE CLASS
# ============================================================

simulate_mixture <- function(n_per_group, g, r, p, separation,
                             noise_sd = 0.5, noise_pi = 0.1) {

  true_means <- lapply(seq_len(g), function(k) {
    matrix(
      rnorm(r * p, mean = separation * (k - (g + 1) / 2)),
      r, p
    )
  })

  x_list <- list()
  true_labels <- integer(0)

  total_n <- n_per_group * g
  n_noise <- round(total_n * noise_pi)

  # ---- signal points ----
  for (k in seq_len(g)) {
    for (i in seq_len(n_per_group)) {
      noise <- matrix(rnorm(r * p, sd = noise_sd), r, p)
      x_list <- c(x_list, list(true_means[[k]] + noise))
      true_labels <- c(true_labels, k)
    }
  }

  # ---- TRUE noise points (important change) ----
  for (i in seq_len(n_noise)) {
    noise_matrix <- matrix(rnorm(r * p, sd = 3 * noise_sd), r, p)
    x_list <- c(x_list, list(noise_matrix))
    true_labels <- c(true_labels, 0)
  }

  list(
    x_list = x_list,
    true_labels = true_labels,
    noise_pi_true = noise_pi
  )
}

# ============================================================
# Fit runner
# ============================================================

run_fit <- function(x_list, g, true_labels, init, noise_pi_true) {

  fit <- tryCatch(
    suppressWarnings(
      matrix_variate_noise_fit(
        x_list,
        g = g,
        noise_type = NOISE_TYPE,
        init = init,
        estimate_k = TRUE,
        max_iter = MAX_ITER,
        tol = TOL,
        nstart = N_RESTARTS,
        verbose = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.frame(
      init = init,
      converged = FALSE,
      ari = NA_real_,
      misclassification_rate = NA_real_,
      noise_pi_true = noise_pi_true,
      noise_pi_est = NA_real_,
      noise_precision = NA_real_,
      noise_recall = NA_real_,
      noise_f1 = NA_real_,
      iterations = NA_real_
    ))
  }

  pred <- fit$cluster
  true <- true_labels

  signal_idx <- true != 0
  noise_idx   <- true == 0

  # ---- clustering quality (signal only) ----
  ari <- tryCatch(
    adjusted_rand_index(true[signal_idx], pred[signal_idx]),
    error = function(e) NA_real_
  )

  misclassification_rate <- mean(pred[signal_idx] != true[signal_idx])

  # ---- noise detection metrics ----
  pred_noise <- pred == 0

  tp <- sum(pred_noise & noise_idx)
  fp <- sum(pred_noise & signal_idx)
  fn <- sum(!pred_noise & noise_idx)

  noise_precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  noise_recall    <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)

  noise_f1 <- if (is.na(noise_precision) || is.na(noise_recall) ||
                  (noise_precision + noise_recall) == 0) {
    NA_real_
  } else {
    2 * noise_precision * noise_recall / (noise_precision + noise_recall)
  }

  data.frame(
    init = init,
    converged = fit$converged,
    ari = ari,
    misclassification_rate = misclassification_rate,
    noise_pi_true = noise_pi_true,
    noise_pi_est = if (!is.null(fit$noise$pi)) fit$noise$pi else NA_real_,
    noise_precision = noise_precision,
    noise_recall = noise_recall,
    noise_f1 = noise_f1,
    iterations = fit$iterations
  )
}

# ============================================================
# Scenarios
# ============================================================

scenarios <- list(
  list(name="Easy",   n=20, g=2, r=2, p=3, sep=4,   sd=0.5, pi=0.05),
  list(name="Medium", n=25, g=3, r=4, p=4, sep=2.5, sd=0.7, pi=0.10),
  list(name="Hard",   n=25, g=3, r=2, p=3, sep=1.5, sd=1.0, pi=0.20)
)

# ============================================================
# Run benchmark
# ============================================================

all_results <- data.frame()

for (scenario in scenarios) {

  cat("\nScenario:", scenario$name, "\n")

  for (trial in seq_len(N_TRIALS)) {

    set.seed(10000 + trial)

    sim <- simulate_mixture(
      n_per_group = scenario$n,
      g = scenario$g,
      r = scenario$r,
      p = scenario$p,
      separation = scenario$sep,
      noise_sd = scenario$sd,
      noise_pi = scenario$pi
    )

    for (init in INIT_METHODS) {

      cat(sprintf(
        "  Trial %d | %s | %s\r",
        trial, scenario$name, init
      ))

      result <- run_fit(
        x_list = sim$x_list,
        g = scenario$g,
        true_labels = sim$true_labels,
        init = init,
        noise_pi_true = sim$noise_pi_true
      )

      result$scenario <- scenario$name
      result$trial <- trial

      all_results <- rbind(all_results, result)
    }
  }
}

# ============================================================
# Summary
# ============================================================

cat("\n\n===== SUMMARY =====\n")

for (scenario in scenarios) {

  cat("\n", scenario$name, "\n")

  sdata <- all_results[all_results$scenario == scenario$name, ]

  for (init in INIT_METHODS) {

    sub <- sdata[sdata$init == init, ]
    if (nrow(sub) == 0) next

    cat(sprintf(
      "  %-10s ARI=%.3f  miss=%.3f  F1=%.3f  pi_err=%.3f  conv=%.2f%%\n",
      init,
      mean(sub$ari, na.rm = TRUE),
      mean(sub$misclassification_rate, na.rm = TRUE),
      mean(sub$noise_f1, na.rm = TRUE),
      mean(abs(sub$noise_pi_est - sub$noise_pi_true), na.rm = TRUE),
      mean(sub$converged, na.rm = TRUE) * 100
    ))
  }
}

# ============================================================
# Save
# ============================================================

write.csv(all_results, OUTPUT_CSV, row.names = FALSE)
cat("\nSaved to:", OUTPUT_CSV, "\n")