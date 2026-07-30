library(ggplot2)
library(clusterGeneration)
library(data.table)

# ──────────────────────────────────────────────────────────────────────────────
# Shared helper — evaluate one k-grid candidate and return log-likelihood
# ──────────────────────────────────────────────────────────────────────────────
evaluate_k_candidate_with_loglik <- function(idx, x_list, r, p, g,
                                              max_iter, tol,
                                              k_grid, nstart, init,
                                              noise_pi_init, verbose) {
  current_k <- k_grid[idx]

  fit_noise <- mv_noise_fit_impl(
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
      k            = current_k,
      logLik       = NA_real_,
      ks_statistic = Inf,
      ks_p.value   = NA_real_,
      n_used       = length(x_clean)
    ))
  }

  fit_clean <- tryCatch(
    mv_mixture_fit(x_list = x_clean, g = g,
                   max_iter = max_iter, tol = tol, verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit_clean)) {
    return(list(
      k            = current_k,
      logLik       = NA_real_,
      ks_statistic = Inf,
      ks_p.value   = NA_real_,
      n_used       = length(x_clean)
    ))
  }

  ks_result <- suppressWarnings(
    tryCatch(
      mv_noise_ks_score(fit_clean, x_clean),
      error = function(e) list(statistic = Inf, p.value = NA_real_,
                                n_used = length(x_clean))
    )
  )

  list(
    k            = current_k,
    logLik       = tail(fit_clean$logLik, 1),
    ks_statistic = ks_result$statistic,
    ks_p.value   = ks_result$p.value,
    n_used       = ks_result$n_used
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Scenario helpers
# ──────────────────────────────────────────────────────────────────────────────
viroli_simulation <- function(r, p, n, n1, n2, n3, M1, M2, M3,
                               U1, V1, U2, V2, U3, V3, n_outliers) {
  simulate_matrix_group <- function(n, mean_mat, row_cov, col_cov) {
    lapply(seq_len(n), function(i) {
      mean_mat + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov
    })
  }
  x_list <- c(
    simulate_matrix_group(n1, M1, U1, V1),
    simulate_matrix_group(n2, M2, U2, V2),
    simulate_matrix_group(n3, M3, U3, V3)
  )
  outlier_idx <- sample(length(x_list), n_outliers)
  for (idx in outlier_idx) {
    x_list[[idx]] <- matrix(sample(x_list[[idx]]), r, p)
  }
  list(x_list = x_list, outlier_idx = outlier_idx)
}

two_group_uniform_simulation <- function(r, p, n1, n2, M1, M2,
                                         n_outliers, n_outliers_perm) {
  simulate_group <- function(n, mean_matrix, row_sd = 0.5, col_sd = 0.5) {
    row_cov <- diag(row_sd, nrow(mean_matrix))
    col_cov <- diag(col_sd, ncol(mean_matrix))
    lapply(seq_len(n), function(i)
      mean_matrix + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov)
  }
  x_list <- c(simulate_group(n1, M1), simulate_group(n2, M2))

  # Uniform entry outliers
  if (n_outliers > 0) {
    contam <- lapply(seq_len(n_outliers),
                     function(i) matrix(runif(r * p, -3, 3), r, p))
    x_list <- c(x_list, contam)
  }
  # Permuted-entry outliers
  if (n_outliers_perm > 0) {
    perm_idx <- sample(seq_along(x_list), n_outliers_perm)
    for (idx in perm_idx) {
      x_list[[idx]] <- matrix(sample(x_list[[idx]]), r, p)
    }
  }
  list(x_list = x_list)
}

# ──────────────────────────────────────────────────────────────────────────────
# Scenario runner
# ──────────────────────────────────────────────────────────────────────────────
run_loglik_scenario <- function(x_list, g, scenario_name,
                                subtitle, r, p,
                                nstart = 10, max_iter = 100, tol = 1e-6,
                                use_kmeans = TRUE, init = "kmeans",
                                noise_pi_init = 0.05,
                                save_plots = TRUE, plots_dir = "loglik/plots") {
  cat(sprintf("\n===== Scenario: %s =====\n", scenario_name))

  cat(sprintf("  Fitting estimate_k with g = %d ...\n", g))
  fit_noise <- mv_noise_fit(
    x_list       = x_list,
    g            = g,
    noise_type   = "hc",
    nstart       = nstart,
    estimate_k   = TRUE,
    init         = init,
    verbose      = TRUE,
    use_parallel = FALSE
  )

  k_selection <- fit_noise$k_selection
  best_k      <- k_selection$selected_k
  best_ks_idx <- which.min(k_selection$ks_scores)

  cat(sprintf("  KS-selected k = %.4e (KS = %.4f)\n",
              best_k, k_selection$ks_scores[best_ks_idx]))

  cat("  Collecting log-likelihoods across k_grid ...\n")
  grid_results <- lapply(
    seq_along(k_selection$k_grid),
    evaluate_k_candidate_with_loglik,
    x_list        = x_list,
    r             = r,
    p             = p,
    g             = g,
    max_iter      = max_iter,
    tol           = tol,
    k_grid        = k_selection$k_grid,
    nstart        = nstart,
    init          = init,
    noise_pi_init = noise_pi_init,
    verbose       = FALSE
  )

  plot_df <- rbindlist(grid_results)

  ok_loglik <- complete.cases(plot_df$logLik)
  if (any(ok_loglik)) {
    max_loglik_k <- plot_df$k[which.max(plot_df$logLik[ok_loglik])]
    plot_df$type  <- rep("Other", nrow(plot_df))
    plot_df$type[plot_df$k == max_loglik_k] <- "Max log-likelihood"
    plot_df$type[plot_df$k == best_k]        <- "Selected k"
  } else {
    plot_df$type <- rep("Other", nrow(plot_df))
    plot_df$type[plot_df$k == best_k] <- "Selected k"
  }

  plot_df$ks_flag <- plot_df$k == k_selection$k_grid[best_ks_idx]

  # Diagnostics
  cat(sprintf("  Grid size:           %d candidates\n", nrow(plot_df)))
  if (any(plot_df$type == "Max log-likelihood")) {
    cat(sprintf("  Max log-likelihood k = %.4e\n",
                plot_df$k[plot_df$type == "Max log-likelihood"][1]))
  } else {
    cat("  Max log-likelihood k = (no valid values)\n")
  }

  cor_val <- NA_real_
  if (any(ok_loglik)) {
    cor_val <- cor(plot_df$k[ok_loglik], plot_df$logLik[ok_loglik],
                   method = "spearman")
    cat(sprintf("  Spearman (k vs logLik): %.4f\n", cor_val))
  }

  # ── Plot 1: log-likelihood vs k (linear y) ─────────────────────────────────
  p1 <- ggplot(plot_df, aes(x = k, y = logLik)) +
    geom_vline(xintercept = best_k, color = "black", linetype = "dashed",
               linewidth = 0.8) +
    geom_vline(xintercept = k_selection$k_grid[best_ks_idx],
               color = "steelblue", linetype = "dotdash", linewidth = 0.8) +
    geom_line(color = "darkorange2", linewidth = 1) +
    geom_point(data = subset(plot_df, type != "Other"),
               aes(shape = type), color = "darkorange2", size = 3) +
    scale_x_log10() +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(
      title    = sprintf("Scenario: %s", scenario_name),
      subtitle = subtitle,
      x        = expression(k~("noise height, log scale")),
      y        = expression( Final~log-likelihood ),
      shape    = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank())

  # ── Plot 2: log-likelihood with log y ──────────────────────────────────────
  p2 <- p1 +
    labs(
      title = sprintf("Scenario: %s (log y)", scenario_name),
      y     = expression( Final~log-likelihood~"(log scale)" )
    ) +
    scale_y_continuous(labels = scales::comma_format())

  # ── Plot 3: log-likelihood + KS combined ───────────────────────────────────
  plot_df_long <- melt(
    plot_df[, .(k, logLik, ks_statistic)],
    id.vars = "k", variable.name = "metric", value.name = "value"
  )
  plot_df_long$metric <- factor(plot_df_long$metric,
    levels = c("logLik", "ks_statistic"),
    labels = c("Log-likelihood", "KS statistic")
  )

  p3 <- ggplot(plot_df_long, aes(x = k, y = value)) +
    geom_vline(xintercept = best_k, color = "black", linetype = "dashed",
               linewidth = 0.8) +
    geom_line(color = "darkorange2", linewidth = 1) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    scale_x_log10() +
    labs(
      title    = sprintf("Scenario: %s — estimate_k diagnostics", scenario_name),
      subtitle = "Black dashed = selected k; top = log-likelihood, bottom = KS",
      x        = expression(k~("noise height, log scale")),
      y        = "Value"
    ) +
    theme_minimal()

  # ── Save plots ─────────────────────────────────────────────────────────────
  if (save_plots) {
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

    tag <- paste0(scenario_name, collapse = "_")
    tag <- gsub("[^a-zA-Z0-9_]+", "_", tag)

    ggsave(
      filename = file.path(plots_dir, sprintf("%s_loglik_trend.png", tag)),
      plot     = p1, width = 9, height = 5.5, dpi = 150
    )
    cat(sprintf("  Saved plots to %s/\n", plots_dir))
  }

  list(
    plot_df  = plot_df,
    p_loglik = p1,
    p_loglik_logy = p2,
    p_combined    = p3,
    best_k        = best_k,
    k_grid        = k_selection$k_grid,
    ks_scores     = k_selection$ks_scores,
    correlation   = cor_val
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# SCENARIO DEFINITIONS
# ══════════════════════════════════════════════════════════════════════════════

set.seed(42)

# ─── Viroli-style scenarios (3×5, 3 groups, permutation outliers) ─────────────
v_r <- 3; v_p <- 5; v_n <- 300

M1_v_base <- matrix(c( 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = v_r, ncol = v_p, byrow = FALSE)
M2_v_base <- matrix(0, v_r, v_p)
M3_v_base <- matrix(c(-0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = v_r, ncol = v_p, byrow = FALSE)

U1v <- rcorrmatrix(v_r); V1v <- rcorrmatrix(v_p)
U2v <- rcorrmatrix(v_r); V2v <- rcorrmatrix(v_p)
U3v <- rcorrmatrix(v_r); V3v <- rcorrmatrix(v_p)

# Scenario V1: small signal (first-column means 0.05, 0 −0.05)
sim_v1 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = matrix(c(0.05, 0.05, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  M2 = M2_v_base, M3 = matrix(c(-0.05, 0.05, 0, rep(0, 12)),
                               nrow = v_r, ncol = v_p, byrow = FALSE),
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V2: reduced covariance (row/col noise scaled to 0.3)
U1v2 <- 0.3 * rcorrmatrix(v_r); V1v2 <- 0.3 * rcorrmatrix(v_p)
U2v2 <- 0.3 * rcorrmatrix(v_r); V2v2 <- 0.3 * rcorrmatrix(v_p)
U3v2 <- 0.3 * rcorrmatrix(v_r); V3v2 <- 0.3 * rcorrmatrix(v_p)
sim_v2 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v2, V1 = V1v2, U2 = U2v2, V2 = V2v2, U3 = U3v2, V3 = V3v2,
  n_outliers = 15
)

# Scenario V3: fewer permutation outliers (5 instead of 15)
set.seed(42 + 1)
sim_v3 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 5
)

# ─── Two-group contamination scenarios (4×4, M1=1, M2=−1) ─────────────────────
c_r <- 4; c_p <- 4
M1_c <- matrix(1,    c_r, c_p)
M2_c <- matrix(-1,   c_r, c_p)

# Scenario C1: low uniform contamination (5 out of 65)
sim_c1 <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = M1_c, M2 = M2_c, n_outliers = 5, n_outliers_perm = 0
)

# Scenario C2: medium positive contamination matching vignette (15 ut of 75)
sim_c2 <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = M1_c, M2 = M2_c, n_outliers = 15, n_outliers_perm = 0
)

# Scenario C3: high contamination with both contamination types (15 uniform + 15 permuted)
sim_c3 <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = M1_c, M2 = M2_c, n_outliers = 15, n_outliers_perm = 15
)

# ══════════════════════════════════════════════════════════════════════════════
# RUN ALL SCENARIOS
# ══════════════════════════════════════════════════════════════════════════════

results <- list()

results$V1_small_signal <- run_loglik_scenario(
  x_list      = sim_v1$x_list,
  g           = 3,
  scenario_name = "Viroli_small_signal",
  subtitle    = "Viroli 3×5 | 3 groups | first-col means 0.05,0,−0.05 | 15 outliers",
  r           = v_r, p = v_p
)

results$V2_low_covariance <- run_loglik_scenario(
  x_list      = sim_v2$x_list,
  g           = 3,
  scenario_name = "Viroli_low_cov",
  subtitle    = "Viroli 3×5 | 3 groups | row/col SD = 0.3 instead of 0.5 | 15 outliers",
  r           = v_r, p = v_p
)

results$V3_few_outliers <- run_loglik_scenario(
  x_list      = sim_v3$x_list,
  g           = 3,
  scenario_name = "Viroli_few_outliers",
  subtitle    = "Viroli 3×5 | 3 groups | 5 permuted-outliers instead of 15",
  r           = v_r, p = v_p
)

results$C1_low_uniform_contam <- run_loglik_scenario(
  x_list      = sim_c1$x_list,
  g           = 2,
  scenario_name = "Contam_low_uniform",
  subtitle    = "Two-group 4×4 | M1=1, M2=-1 | 5 uniform outliers (Unif(-3,3)) | n=65",
  r           = c_r, p = c_p
)

results$C2_medium_uniform_contam <- run_loglik_scenario(
  x_list      = sim_c2$x_list,
  g           = 2,
  scenario_name = "Contam_medium_uniform",
  subtitle    = "Two-group 4×4 | M1=1, M2=-1 | 15 uniform outliers (Unif(-3,3)) | n=75",
  r           = c_r, p = c_p
)

results$C3_high_mixed_contam <- run_loglik_scenario(
  x_list      = sim_c3$x_list,
  g           = 2,
  scenario_name = "Contam_high_mixed",
  subtitle    = "Two-group 4×4 | M1=1, M2=-1 | 15 uniform + 15 permuted outliers | n=90",
  r           = c_r, p = c_p
)

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════════════════

summary_data <- rbindlist(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.table(
    Scenario      = nm,
    Selected_k    = r$best_k,
    Max_logLik_k  = ifelse(any(r$plot_df$logLik[complete.cases(r$plot_df$logLik)] > -Inf),
                          r$plot_df$k[complete.cases(r$plot_df$logLik)][
                            which.max(r$plot_df$logLik[complete.cases(r$plot_df$logLik)])
                          ], NA_real_),
    Spearman_rho  = r$correlation,
    N_candidates  = nrow(r$plot_df),
    N_valid_logLik= sum(complete.cases(r$plot_df$logLik))
  )
}), fill = TRUE)

colnames(summary_data)[2] <- "Selected_k"
colnames(summary_data)[3] <- "Max_logLik_k"

cat("\n===== Cross-scenario summary =====\n")
print(summary_data)

cat(sprintf("\nAll 18 plot files saved under %s/\n", normalizePath("loglik/plots")))
