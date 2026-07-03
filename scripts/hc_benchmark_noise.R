#!/usr/bin/env Rscript
# Benchmark: Noise Mixture Clustering (hc only, estimate_k only, no KS sweep)

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
# Helpers
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

simulate_mixture <- function(n_per_group, g, r, p, separation, noise_sd = 0.5) {

  true_means <- lapply(seq_len(g), function(k) {
    matrix(
      rnorm(r * p, mean = separation * (k - (g + 1) / 2)),
      r, p
    )
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

run_fit <- function(x_list, g, true_labels, init) {

  fit <- tryCatch(
    suppressWarnings(
      matrix_variate_noise_fit(
        x_list,
        g = g,
        noise_type = NOISE_TYPE,
        init = init,
        ks_type = "onesample",   # fixed placeholder (required arg)
        noise_k = 1e-4,
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
      noise_pi = NA_real_,
      iterations = NA_real_
    ))
  }

  keep <- fit$cluster > 0

  ari <- tryCatch(
    adjusted_rand_index(true_labels[keep], fit$cluster[keep]),
    error = function(e) NA_real_
  )

  data.frame(
    init = init,
    converged = fit$converged,
    ari = ari,
    noise_pi = fit$noise$pi,
    iterations = fit$iterations
  )
}

# ============================================================
# Scenarios
# ============================================================

scenarios <- list(
  list(name="Easy",   n=20, g=2, r=2, p=3, sep=4,   sd=0.5),
  list(name="Medium", n=25, g=3, r=4, p=4, sep=2.5, sd=0.7),
  list(name="Hard",   n=25, g=3, r=2, p=3, sep=1.5, sd=1.0)
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
      noise_sd = scenario$sd
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
        init = init
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

cat("\n\n===== SUMMARY (Mean ARI) =====\n")

for (scenario in scenarios) {

  cat("\n", scenario$name, "\n")

  sdata <- all_results[all_results$scenario == scenario$name, ]

  for (init in INIT_METHODS) {

    sub <- sdata[sdata$init == init, ]
    if (nrow(sub) == 0) next

    cat(sprintf(
      "  %-10s ARI=%.3f  conv=%.2f%%  pi=%.3f  iters=%.1f\n",
      init,
      mean(sub$ari, na.rm = TRUE),
      mean(sub$converged, na.rm = TRUE) * 100,
      mean(sub$noise_pi, na.rm = TRUE),
      mean(sub$iterations, na.rm = TRUE)
    ))
  }
}

# ============================================================
# Save
# ============================================================

write.csv(all_results, OUTPUT_CSV, row.names = FALSE)
cat("\nSaved to:", OUTPUT_CSV, "\n")