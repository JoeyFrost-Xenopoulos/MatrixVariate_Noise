#' Log-likelihood trend across estimate_k candidates (Viroli simulation)
#'
#' Log-likelihood trend across estimate_k candidates (Viroli simulation)
#'
#' Demonstrates how the observed-data log-likelihood of the HC fit changes
#' across the k-grid used by estimate_k. This is meant as a diagnostic for
#' the selection path, not as a separate post-cleaning refit analysis.
#'
#' Uses the same 3-group Viroli simulation as in Noise_Tests.qmd.
#'
#' @param n_iter  Number of Monte Carlo replicates
#' @return list(data.table plot_df)
#' @export
library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)

set.seed(42)

r_viroli  <- 3
p_viroli  <- 5
n_viroli  <- 300

M1_viroli <- matrix(c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = r_viroli, ncol = p_viroli, byrow = FALSE)
M2_viroli <- matrix(0, r_viroli, p_viroli)
M3_viroli <- matrix(c(-0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = r_viroli, ncol = p_viroli, byrow = FALSE)

U1_viroli <- rcorrmatrix(r_viroli); V1_viroli <- rcorrmatrix(p_viroli)
U2_viroli <- rcorrmatrix(r_viroli); V2_viroli <- rcorrmatrix(p_viroli)
U3_viroli <- rcorrmatrix(r_viroli); V3_viroli <- rcorrmatrix(p_viroli)

n1_v <- round(0.3 * n_viroli)
n2_v <- round(0.4 * n_viroli)
n3_v <- n_viroli - n1_v - n2_v

simulate_matrix_group <- function(n, mean_mat, row_cov, col_cov) {
  lapply(seq_len(n), function(i) {
    mean_mat + row_cov %*% matrix(rnorm(r_viroli * p_viroli),
                                  r_viroli, p_viroli) %*% col_cov
  })
}

x_list_viroli <- c(
  simulate_matrix_group(n1_v, M1_viroli, U1_viroli, V1_viroli),
  simulate_matrix_group(n2_v, M2_viroli, U2_viroli, V2_viroli),
  simulate_matrix_group(n3_v, M3_viroli, U3_viroli, V3_viroli)
)

outlier_idx_viroli <- sample(length(x_list_viroli), 15)
for (idx in outlier_idx_viroli) {
  x_list_viroli[[idx]] <- matrix(sample(x_list_viroli[[idx]]),
                                  r_viroli, p_viroli)
}

cat(sprintf("Viroli simulation: n = %d, groups = %d %d %d\n",
            length(x_list_viroli), n1_v, n2_v, n3_v))

# ──────────────────────────────────────────────────────────────────────────────
# Internal evaluation with final log-likelihood capture
# ──────────────────────────────────────────────────────────────────────────────
evaluate_k_candidate_with_loglik <- function(idx, x_list, g, max_iter, tol,
                                              k_grid, nstart, init,
                                              noise_pi_init, verbose) {
  current_k <- k_grid[idx]

  fit_noise <- Ampharos:::mv_noise_fit_impl(
    x_list        = x_list,
    g             = g,
    noise_type    = "hc",
    max_iter      = max_iter,
    tol           = tol,
    nstart        = nstart,
    noise_k       = current_k,
    noise_pi_init = noise_pi_init,
    init          = init,
    verbose       = FALSE
  )

  keep_idx <- fit_noise$cluster != 0
  x_clean  <- x_list[keep_idx]

  if (length(x_clean) <= g) {
    return(list(
      k             = current_k,
      logLik        = tail(fit_noise$logLik, 1),
      ks_statistic  = Inf,
      ks_p.value    = NA_real_,
      n_used        = length(x_clean)
    ))
  }

  ks_result <- suppressWarnings(
    tryCatch(
      Ampharos:::mv_noise_ks_score(fit_noise, x_clean),
      error = function(e) list(statistic = Inf, p.value = NA_real_,
                                n_used = length(x_clean))
    )
  )

  list(
    k             = current_k,
    logLik        = tail(fit_noise$logLik, 1),
    ks_statistic  = ks_result$statistic,
    ks_p.value    = ks_result$p.value,
    n_used        = ks_result$n_used
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Run estimate_k on Viroli data
# ──────────────────────────────────────────────────────────────────────────────
cat("Running mv_noise_fit with estimate_k = TRUE on Viroli data...\n")
fit_hc_viroli <- Ampharos::mv_noise_fit(
  x_list       = x_list_viroli,
  g            = 3,
  noise_type   = "hc",
  nstart       = 10,
  estimate_k   = TRUE,
  init         = "kmeans",
  verbose      = TRUE,
  use_parallel = FALSE
)

k_selection <- fit_hc_viroli$k_selection
best_k <- k_selection$selected_k

cat(sprintf("\nSelected k  = %.4e  (KS = %.4f)\n",
            best_k, k_selection$ks_scores[which.min(k_selection$ks_scores)]))

# ──────────────────────────────────────────────────────────────────────────────
# Re-run grid with log-likelihood capture
# ──────────────────────────────────────────────────────────────────────────────
cat("\nCollecting final log-likelihoods across the k_grid ...\n")

grid_results <- lapply(
  seq_along(k_selection$k_grid),
  evaluate_k_candidate_with_loglik,
  x_list        = x_list_viroli,
  g             = 3,
  max_iter      = 100,
  tol           = 1e-6,
  k_grid        = k_selection$k_grid,
  nstart        = 10,
  init          = "kmeans",
  noise_pi_init = 0.05,
  verbose       = FALSE
)

plot_df <- rbindlist(grid_results)

# Annotate rows
ok_loglik <- complete.cases(plot_df$logLik)
if (any(ok_loglik)) {
  max_loglik_k <- plot_df$k[which.max(plot_df$logLik[ok_loglik])]
  plot_df$type[plot_df$k == max_loglik_k & ok_loglik] <- "Max log-likelihood"
}
plot_df$ks_flag <- plot_df$k == k_selection$k_grid[which.min(k_selection$ks_scores)]

cat(sprintf("Grid size: %d candidates\n", nrow(plot_df)))
cat(sprintf("Best KS:   k = %.4e\n", k_selection$k_grid[which.min(k_selection$ks_scores)]))
if (any(plot_df$type == "Max log-likelihood")) {
  cat(sprintf("Max logLik k = %.4e\n", plot_df$k[plot_df$type == "Max log-likelihood"][1]))
} else {
  cat("Max logLik k = (no valid log-likelihoods)\n")
}

# ──────────────────────────────────────────────────────────────────────────────
# Plot 1: log-likelihood vs k (linear y)
# ──────────────────────────────────────────────────────────────────────────────
p_loglik <- ggplot(plot_df, aes(x = k, y = logLik)) +
  geom_vline(xintercept = best_k, color = "black", linetype = "dashed",
             linewidth = 0.8) +
  geom_vline(xintercept = k_selection$k_grid[which.min(k_selection$ks_scores)],
             color = "steelblue", linetype = "dotdash", linewidth = 0.8) +
  geom_line(color = "darkorange2", linewidth = 1) +
  geom_point(data = subset(plot_df, type != "Other"),
             aes(shape = type), color = "darkorange2", size = 3) +
  scale_x_log10() +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    title = "Final log-likelihood vs noise_k across estimate_k grid",
    subtitle = "Viroli simulation (3 groups, 15 outliers, n = 300, 3x5 matrices)",
    x = expression(k~("noise height, log scale")),
    y = expression( Final~log-likelihood ),
    shape = ""
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank()
  )

# ──────────────────────────────────────────────────────────────────────────────
# Plot 2: log-likelihood vs k (log y) to accentuate spread
# ──────────────────────────────────────────────────────────────────────────────
p_loglik_logy <- p_loglik +
  labs(
    title = "Final log-likelihood vs noise_k (alternate linear-scale view)",
    y = expression( Final~log-likelihood~"(linear scale)" )
  ) +
  scale_y_continuous(labels = scales::comma_format())

# ──────────────────────────────────────────────────────────────────────────────
# Plot 3: combined view — log-likelihood (left) + KS score (right)
# ──────────────────────────────────────────────────────────────────────────────
plot_df_long <- melt(
  plot_df[, .(k, logLik, ks_statistic)],
  id.vars = "k", variable.name = "metric", value.name = "value"
)
plot_df_long$metric <- factor(plot_df_long$metric,
  levels = c("logLik", "ks_statistic"),
  labels = c("Log-likelihood", "KS statistic")
)

p_combined <- ggplot(plot_df_long, aes(x = k, y = value)) +
  geom_vline(xintercept = best_k, color = "black", linetype = "dashed",
             linewidth = 0.8) +
  geom_line(color = "darkorange2", linewidth = 1) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_x_log10() +
  labs(
    title    = "estimate_k diagnostics across the k-grid (Viroli simulation)",
    subtitle = "Black dashed line = selected k; top = log-likelihood, bottom = KS statistic",
    x        = expression(k~("noise height, log scale")),
    y        = "Value"
  ) +
  theme_minimal()

# ──────────────────────────────────────────────────────────────────────────────
# Print diagnostics and save
# ──────────────────────────────────────────────────────────────────────────────
cat("\n--- trend summary ---\n")
ok <- complete.cases(plot_df$logLik)
if (any(ok)) {
  cor_val <- cor(plot_df$k[ok], plot_df$logLik[ok], method = "spearman")
  cat(sprintf("Spearman rank correlation (k vs log-likelihood): %.4f\n", cor_val))
}
cat(sprintf("KS-selected k:                   %.4e\n",
            k_selection$k_grid[which.min(k_selection$ks_scores)]))
if (any(ok_loglik)) {
  cat(sprintf("Log-likelihood-maximising k:     %.4e\n",
              plot_df$k[which.max(plot_df$logLik[ok_loglik])]))
} else {
  cat("Log-likelihood-maximising k:     (no valid values)\n")
}

print(p_loglik)
print(p_loglik_logy)
print(p_combined)
