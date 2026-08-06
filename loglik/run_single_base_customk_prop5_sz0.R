library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)

# ──────────────────────────────────────────────────────────────────────────────
# Custom k-grid (exact match for base_customk scenarios)
# ──────────────────────────────────────────────────────────────────────────────
custom_k_grid <- 10^seq(-100, -4, length.out = 60)
cat(sprintf("Custom k_grid: [%.0e, %.0e] | %d candidates\n\n",
            min(custom_k_grid), max(custom_k_grid), length(custom_k_grid)))

# ──────────────────────────────────────────────────────────────────────────────
# Simulation helpers (from run_all_base.R)
# ──────────────────────────────────────────────────────────────────────────────
viroli_simulation <- function(r, p, n, n1, n2, n3, M1, M2, M3,
                              U1, V1, U2, V2, U3, V3, n_outliers,
                              row_cov_scale = 1, col_cov_scale = 1) {
  U1 <- row_cov_scale * U1; V1 <- col_cov_scale * V1
  U2 <- row_cov_scale * U2; V2 <- col_cov_scale * V2
  U3 <- row_cov_scale * U3; V3 <- col_cov_scale * V3
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

scaled_viroli_simulation <- function(r, p, n, n1, n2, n3, n_outliers,
                                      signal_strength = 0.5, cov_scale = 1,
                                      outlier_type = c("perm", "row_spike", "col_out"),
                                      seed = NULL) {
  outlier_type <- match.arg(outlier_type)
  if (!is.null(seed)) set.seed(seed)
  M1 <- matrix(0, r, p)
  M2 <- matrix(0, r, p)
  M3 <- matrix(0, r, p)
  signal_col <- min(2, p)
  for (j in seq_len(signal_col)) {
    M1[1, j] <- signal_strength
    M3[1, j] <- -signal_strength
  }
  U1 <- cov_scale * rcorrmatrix(r)
  U2 <- cov_scale * rcorrmatrix(r)
  U3 <- cov_scale * rcorrmatrix(r)
  V1 <- cov_scale * rcorrmatrix(p)
  V2 <- cov_scale * rcorrmatrix(p)
  V3 <- cov_scale * rcorrmatrix(p)
  result <- viroli_simulation(r, p, n, n1, n2, n3,
                              M1, M2, M3, U1, V1, U2, V2, U3, V3,
                              n_outliers = n_outliers)
  if (outlier_type == "row_spike") {
    spike_idx <- sample(seq_along(result$x_list), n_outliers)
    for (idx in spike_idx) {
      row_id <- sample.int(r, 1)
      result$x_list[[idx]][row_id, ] <- runif(p, -8, 8)
    }
  }
  if (outlier_type == "col_out") {
    col_out_idx <- sample(seq_along(result$x_list), n_outliers)
    for (idx in col_out_idx) {
      col_id <- sample.int(p, 1)
      result$x_list[[idx]][, col_id] <- runif(r, -10, 10)
    }
  }
  list(x_list = result$x_list, outlier_idx = result$outlier_idx)
}

# ──────────────────────────────────────────────────────────────────────────────
# Evaluate one k-grid candidate and return log-likelihood + noise-recovery flag
# ──────────────────────────────────────────────────────────────────────────────
evaluate_k_candidate_with_loglik <- function(idx, x_list, r, p, g,
                                               max_iter, tol,
                                               k_grid, nstart, init,
                                               noise_pi_init, verbose,
                                               outlier_idx) {
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
    verbose       = FALSE,
    use_parallel  = FALSE
  )

  keep_idx <- fit_noise$cluster != 0
  x_clean  <- x_list[keep_idx]

  detected_noise <- which(fit_noise$cluster == 0)
  correct_noise <- length(detected_noise) == length(outlier_idx) &&
    length(intersect(detected_noise, outlier_idx)) == length(outlier_idx)

  if (length(x_clean) <= g) {
    return(list(
      k            = current_k,
      logLik       = tail(fit_noise$logLik, 1),
      ks_statistic = Inf,
      ks_p.value   = NA_real_,
      n_used       = length(x_clean),
      correct_noise = correct_noise
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
    k            = current_k,
    logLik       = tail(fit_noise$logLik, 1),
    ks_statistic = ks_result$statistic,
    ks_p.value   = ks_result$p.value,
    n_used       = ks_result$n_used,
    correct_noise = correct_noise
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Run a single base scenario (no replicates, no parallel)
# ──────────────────────────────────────────────────────────────────────────────
run_base_scenario <- function(x_list, g, scenario_name,
                               subtitle, r, p,
                               true_k_noise = NA_integer_,
                               outlier_idx = NULL,
                               nstart = 10, max_iter = 100, tol = 1e-6,
                               init = "kmeans",
                               noise_pi_init = 0.05,
                               save_plots = TRUE, plots_dir = "loglik/plots",
                               seed = 42, k_grid = custom_k_grid) {
  if (!is.null(seed)) set.seed(seed)

  cat(sprintf("\n===== Scenario: %s =====\n", scenario_name))
  cat(sprintf("  Subtitle : %s\n", subtitle))
  cat(sprintf("  r=%d, p=%d, g=%d, n=%d\n\n", r, p, g, length(x_list)))

  # ── Step 1: KS-based k selection ──────────────────────────────────────────
  cat("  [Step 1] Fitting estimate_k with k_grid ...\n")
  fit_noise <- mv_noise_fit(
    x_list       = x_list,
    g            = g,
    noise_type   = "hc",
    nstart       = nstart,
    estimate_k   = TRUE,
    k_grid       = k_grid,
    init         = init,
    verbose      = TRUE,
    use_parallel = FALSE
  )

  k_selection <- fit_noise$k_selection
  best_k      <- k_selection$selected_k
  best_ks_idx <- which.min(k_selection$ks_scores)

  cat(sprintf("\n  KS-selected k = %.4e (KS = %.4f, p = %.4f)\n",
              best_k,
              k_selection$ks_scores[best_ks_idx],
              k_selection$ks_pvalues[best_ks_idx]))

  # ── Step 2: Evaluate log-likelihood across the full k_grid ────────────────
  cat("\n  [Step 2] Evaluating log-likelihood across k_grid ...\n\n")
  active_k_grid <- k_selection$k_grid

  grid_results <- lapply(
    seq_along(active_k_grid),
    evaluate_k_candidate_with_loglik,
    x_list        = x_list,
    r             = r,
    p             = p,
    g             = g,
    max_iter      = max_iter,
    tol           = tol,
    k_grid        = active_k_grid,
    nstart        = nstart,
    init          = init,
    noise_pi_init = noise_pi_init,
    verbose       = FALSE,
    outlier_idx   = outlier_idx
  )

  plot_df <- rbindlist(grid_results)
  plot_df$scenario_name <- scenario_name

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

  # ── Step 3: Final fit with selected k ─────────────────────────────────────
  cat("\n  [Step 3] Final fit with selected k ...\n")
  fit_final <- Ampharos:::mv_noise_fit_impl(
    x_list        = x_list,
    g             = g,
    noise_type    = "hc",
    max_iter      = max_iter,
    tol           = tol,
    nstart        = nstart,
    noise_k       = best_k,
    noise_pi_init = noise_pi_init,
    init          = init,
    verbose       = FALSE,
    use_parallel  = FALSE
  )

  n_total <- length(x_list)
  detected_noise_idx <- which(fit_final$cluster == 0)
  tp <- length(intersect(detected_noise_idx, outlier_idx))
  fp <- length(setdiff(detected_noise_idx, outlier_idx))
  fn <- length(setdiff(outlier_idx, detected_noise_idx))
  tn <- n_total - tp - fp - fn

  accuracy_rate <- if (n_total > 0) (tp + tn) / n_total else NA_real_
  fpr <- if ((fp + tn) > 0) fp / (fp + tn) else NA_real_
  fnr <- if ((fn + tp) > 0) fn / (fn + tp) else NA_real_

  logLik_best <- plot_df$logLik[plot_df$k == best_k][1]
  logLik_max  <- if (any(ok_loglik)) max(plot_df$logLik[ok_loglik]) else NA_real_
  logLik_true <- if (!is.na(true_k_noise) && true_k_noise %in% plot_df$k)
    plot_df$logLik[plot_df$k == true_k_noise][1] else NA_real_

  loglik_gap_vs_max  <- if (!is.na(logLik_best) && !is.na(logLik_max))
    logLik_max - logLik_best else NA_real_
  loglik_gap_vs_true <- if (!is.na(logLik_best) && !is.na(logLik_true))
    logLik_true - logLik_best else NA_real_

  metrics <- data.table(
    scenario             = scenario_name,
    n_total              = n_total,
    true_k_noise         = ifelse(!is.na(true_k_noise), true_k_noise, NA_integer_),
    best_k               = best_k,
    best_ks_k            = k_selection$k_grid[best_ks_idx],
    detected_noise_count = length(detected_noise_idx),
    tp = tp, fp = fp, fn = fn, tn = tn,
    accuracy_rate        = accuracy_rate,
    fpr                  = fpr,
    fnr                  = fnr,
    logLik_best_k        = logLik_best,
    logLik_max           = logLik_max,
    logLik_true_k        = logLik_true,
    loglik_gap_vs_max    = loglik_gap_vs_max,
    loglik_gap_vs_true   = loglik_gap_vs_true
  )

  if (!is.na(true_k_noise)) {
    selected_idx <- match(best_k, plot_df$k)
    correct_noise_detected <- !is.na(selected_idx) &&
      isTRUE(plot_df$correct_noise[selected_idx])
    cat(sprintf("  True noise k = %.4e | Correctly identified = %s\n\n",
                true_k_noise, if (correct_noise_detected) "YES" else "NO"))
    metrics$correct_noise_detected <- correct_noise_detected
  }

  cor_val <- NA_real_
  if (any(ok_loglik)) {
    cor_val <- cor(plot_df$k[ok_loglik], plot_df$logLik[ok_loglik], method = "spearman")
    cat(sprintf("  Spearman (k vs logLik): %.4f\n", cor_val))
  }

  # ── Step 4: Plot ───────────────────────────────────────────────────────────
  p1 <- ggplot() +
    {
      ok_rows <- plot_df[complete.cases(plot_df$logLik), ]
      if (nrow(ok_rows) > 1) {
        ok_rows <- ok_rows[order(ok_rows$k), ]
        seg_df <- data.frame(
          x            = head(ok_rows$k, -1),
          xend         = tail(ok_rows$k, -1),
          y            = head(ok_rows$logLik, -1),
          yend         = tail(ok_rows$logLik, -1),
          correct_noise = factor(head(ok_rows$correct_noise, -1),
                                 levels = c(FALSE, TRUE))
        )
        list(geom_segment(
          data = seg_df,
          aes(x = x, xend = xend, y = y, yend = yend, color = correct_noise),
          linewidth = 1
        ))
      } else NULL
    } +
    geom_point(data = subset(plot_df, type != "Other"),
               aes(x = k, y = logLik, shape = type), color = "darkorange2", size = 3) +
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "green"),
                       labels = c("FALSE" = "Incorrect recovery",
                                  "TRUE"  = "Correct noise recovery"),
                       name = "Noise recovery") +
    scale_x_log10() +
    scale_y_continuous(
      name = expression(Final~log-likelihood),
      labels = scales::comma_format()
    ) +
    labs(
      title    = sprintf("Scenario: %s", scenario_name),
      subtitle = subtitle,
      x        = expression(k~("noise height, log scale")),
      shape    = "",
      colour   = "Noise recovery"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_blank())

  if (save_plots) {
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
    ggsave(filename = file.path(plots_dir, paste0(scenario_name, ".png")),
           plot = p1, width = 10, height = 6, dpi = 150)
    cat(sprintf("  Plot saved to %s\n", file.path(plots_dir, paste0(scenario_name, ".png"))))
  }

  cat("\n===== Done =====\n\n")
  print(metrics)

  list(
    plot_df    = plot_df,
    p_loglik   = p1,
    best_k     = best_k,
    k_grid     = active_k_grid,
    ks_scores  = k_selection$ks_scores,
    correlation = cor_val,
    metrics    = metrics,
    plots_dir  = plots_dir
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Exact scenario: base_customk_prop5_sz0
#   size_grid index 0 → r=2, p=3
#   5% of n=150 → 8 permuted outliers
#   simulation seed = 42 + 700 + s=1 + noise_prop=5 = 748
# ──────────────────────────────────────────────────────────────────────────────
rg       <- 2
pg       <- 3
n_total  <- 150
noise_prop <- 5
n_outliers <- min(max(0, round(n_total * noise_prop / 100)), n_total - 3)
sim_seed   <- 42 + 700 + 1 + noise_prop   # 748

cat(sprintf("Scenario: base_customk_prop5_sz0\n"))
cat(sprintf("  r=%d, p=%d, g=3, n=%d, outliers=%d (%.0f%%)\n",
            rg, pg, n_total, n_outliers, noise_prop))
cat(sprintf("  Simulation seed = %d\n\n", sim_seed))

set.seed(sim_seed)
sim_viroli_prop <- scaled_viroli_simulation(
  r = rg, p = pg, n = n_total,
  n1 = round(0.3 * n_total), n2 = round(0.4 * n_total),
  n3 = n_total - round(0.3 * n_total) - round(0.4 * n_total),
  n_outliers = n_outliers, signal_strength = 0.5, cov_scale = 1,
  outlier_type = "perm"
)

scenario_name <- "base_customk_prop5_sz0"
plots_dir     <- file.path("loglik", "plots_base_custom_k", "size_2x3")

result <- run_base_scenario(
  x_list       = sim_viroli_prop$x_list,
  g            = 3,
  scenario_name = scenario_name,
  subtitle     = sprintf("Base Viroli %dx%d | 3 groups | %d permuted outliers (%.0f%% of n=%d)",
                         rg, pg, n_outliers, noise_prop, n_total),
  r            = rg,
  p            = pg,
  true_k_noise = n_outliers,
  outlier_idx  = sim_viroli_prop$outlier_idx,
  plots_dir    = plots_dir,
  k_grid       = custom_k_grid,
  seed         = sim_seed,
  nstart       = 10,
  max_iter     = 100,
  tol          = 1e-6,
  init         = "kmeans",
  noise_pi_init = 0.05
)

cat(sprintf("\nFinal result:\n"))
cat(sprintf("  best_k             = %.4e\n", result$best_k))
cat(sprintf("  correlation (k/logLik) = %.4f\n", result$correlation))
cat(sprintf("  accuracy_rate      = %.4f\n", result$metrics$accuracy_rate))
cat(sprintf("  logLik_gap_vs_true = %.4f\n", result$metrics$loglik_gap_vs_true))
