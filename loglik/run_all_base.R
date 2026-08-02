library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)
library(future)
library(future.apply)

plan(multisession, workers = 8)

# ──────────────────────────────────────────────────────────────────────────────
# Simulation helpers
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
# Shared helper — evaluate one k-grid candidate and return log-likelihood
# ──────────────────────────────────────────────────────────────────────────────
evaluate_k_candidate_with_loglik <- function(idx, x_list, r, p, g,
                                               max_iter, tol,
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
      k      = current_k,
      logLik = tail(fit_noise$logLik, 1),
      ks_statistic = Inf,
      ks_p.value   = NA_real_,
      n_used       = length(x_clean)
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
    n_used       = ks_result$n_used
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Run a single base scenario
# ──────────────────────────────────────────────────────────────────────────────
run_base_scenario <- function(x_list, g, scenario_name,
                               subtitle, r, p,
                               true_k_noise = NA_integer_,
                               outlier_idx = NULL,
                               nstart = 10, max_iter = 100, tol = 1e-6,
                               use_kmeans = TRUE, init = "kmeans",
                               noise_pi_init = 0.05,
                               save_plots = TRUE, plots_dir = "loglik/plots",
                               seed = 42) {
  if (!is.null(seed)) set.seed(seed)

  for (pkg in c("Ampharos", "ggplot2", "clusterGeneration", "data.table")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package ", pkg, " required")
  }

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
    verbose       = FALSE
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
    verbose       = FALSE
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
  logLik_max <- if (any(ok_loglik)) max(plot_df$logLik[ok_loglik]) else NA_real_
  logLik_true <- if (!is.na(true_k_noise) && true_k_noise %in% plot_df$k)
                   plot_df$logLik[plot_df$k == true_k_noise][1] else NA_real_

  loglik_gap_vs_max <- if (!is.na(logLik_best) && !is.na(logLik_max))
                         logLik_max - logLik_best else NA_real_
  loglik_gap_vs_true <- if (!is.na(logLik_best) && !is.na(logLik_true))
                          logLik_true - logLik_best else NA_real_

  metrics <- data.table(
    scenario = scenario_name,
    n_total = n_total,
    true_k_noise = ifelse(!is.na(true_k_noise), true_k_noise, NA_integer_),
    best_k = best_k,
    best_ks_k = k_selection$k_grid[best_ks_idx],
    detected_noise_count = length(detected_noise_idx),
    tp = tp, fp = fp, fn = fn, tn = tn,
    accuracy_rate = accuracy_rate,
    fpr = fpr,
    fnr = fnr,
    logLik_best_k = logLik_best,
    logLik_max = logLik_max,
    logLik_true_k = logLik_true,
    loglik_gap_vs_max = loglik_gap_vs_max,
    loglik_gap_vs_true = loglik_gap_vs_true
  )

  if (!is.na(true_k_noise)) {
    correct_noise_detected <- best_k == true_k_noise
    cat(sprintf("  True noise count = %.4e | Correctly identified = %s\n",
                true_k_noise, if (correct_noise_detected) "YES" else "NO"))
    metrics$correct_noise_detected <- correct_noise_detected
  }

  cor_val <- NA_real_
  if (any(ok_loglik)) {
    cor_val <- cor(plot_df$k[ok_loglik], plot_df$logLik[ok_loglik],
                   method = "spearman")
    cat(sprintf("  Spearman (k vs logLik): %.4f\n", cor_val))
  }

  ok_ks <- plot_df[is.finite(plot_df$ks_statistic), ]
  ks_scale <- 1
  ks_offset <- 0
  if (nrow(ok_ks) > 0 && any(ok_loglik)) {
    logLik_range <- range(plot_df$logLik[ok_loglik], na.rm = TRUE)
    ks_range <- range(plot_df$ks_statistic[ok_ks], na.rm = TRUE)
    ks_scale <- diff(logLik_range) / diff(ks_range)
    ks_offset <- logLik_range[1] - ks_range[1] * ks_scale
  }

  p1 <- ggplot() +
    {
      ok_rows <- plot_df[complete.cases(plot_df$logLik), ]
      if (nrow(ok_rows) > 1) {
        ok_rows <- ok_rows[order(plot_df$k), ]
        seg_df <- data.frame(
          x      = head(plot_df$k, -1),
          xend   = tail(plot_df$k, -1),
          y      = head(plot_df$logLik, -1),
          yend   = tail(plot_df$logLik, -1)
        )
        list(geom_segment(
          data = seg_df,
          aes(x = x, xend = xend, y = y, yend = yend),
          color = "grey60", linewidth = 1
        ))
      } else NULL
    } +
    geom_point(data = subset(plot_df, type != "Other"),
               aes(x = k, y = logLik, shape = type), color = "darkorange2", size = 3) +
    {
      if (nrow(ok_ks) > 0 && any(ok_loglik)) {
        list(
          geom_line(data = ok_ks, aes(x = k, y = ks_statistic * ks_scale + ks_offset),
                    color = "steelblue", linewidth = 1),
          geom_point(data = ok_ks, aes(x = k, y = ks_statistic * ks_scale + ks_offset),
                     color = "steelblue", size = 2)
        )
      } else NULL
    } +
    scale_color_manual(values = c("FALSE" = "red", "TRUE" = "green"),
                       labels = c("FALSE" = "Incorrect recovery",
                                  "TRUE"  = "Correct noise recovery"),
                       name = "Noise recovery") +
    scale_x_log10() +
    scale_y_continuous(
      name = expression( Final~log-likelihood ),
      labels = scales::comma_format(),
      sec.axis = if (nrow(ok_ks) > 0 && any(ok_loglik))
        sec_axis(~(. - ks_offset) / ks_scale, name = "KS statistic")
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
  }

  list(
    plot_df   = plot_df,
    p_loglik  = p1,
    best_k    = best_k,
    k_grid    = active_k_grid,
    ks_scores = k_selection$ks_scores,
    correlation = cor_val,
    metrics   = metrics,
    plots_dir = plots_dir
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Run replicates of a base scenario
# ──────────────────────────────────────────────────────────────────────────────
run_base_replicates <- function(x_list, g, B = 10, base_seed = 42, ...) {
  if (B < 2) {
    res <- run_base_scenario(x_list = x_list, g = g, seed = base_seed, ...)
    res$B <- 1L
    return(res)
  }

  replicate_results <- future_lapply(seq_len(B), function(b) {
    run_base_scenario(
      x_list = x_list,
      g = g,
      seed = base_seed + b * 1000,
      ...
    )
  }, future.seed = TRUE)

  all_plot_df <- rbindlist(lapply(seq_along(replicate_results), function(b) {
    df <- replicate_results[[b]]$plot_df
    df$replicate <- b
    df
  }))

  agg <- all_plot_df[, .(
    logLik_mean = mean(logLik, na.rm = TRUE),
    logLik_se = ifelse(sum(!is.na(logLik)) > 1,
                       sd(logLik, na.rm = TRUE) / sqrt(sum(!is.na(logLik))),
                       NA_real_),
    ks_statistic_mean = mean(ks_statistic, na.rm = TRUE),
    ks_statistic_se = ifelse(sum(!is.na(ks_statistic) & is.finite(ks_statistic)) > 1,
                             sd(ks_statistic, na.rm = TRUE) / sqrt(sum(!is.na(ks_statistic) & is.finite(ks_statistic))),
                             NA_real_),
    n_used_mean = mean(n_used, na.rm = TRUE),
    scenario_name = scenario_name[1],
    .N
  ), by = k]

  ok_loglik <- !is.na(agg$logLik_mean) & is.finite(agg$logLik_mean)
  best_k_overall <- if (any(ok_loglik)) agg$k[which.max(agg$logLik_mean[ok_loglik])] else NA_real_

  p <- ggplot() +
    geom_line(data = all_plot_df, aes(x = k, y = logLik, group = replicate),
              color = "gray60", alpha = 0.4, linewidth = 0.5) +
    {
      if (any(ok_loglik)) {
        list(
          geom_ribbon(data = agg[ok_loglik, ],
                      aes(x = k, ymin = logLik_mean - logLik_se, ymax = logLik_mean + logLik_se),
                      fill = "darkorange2", alpha = 0.3),
          geom_line(data = agg[ok_loglik, ], aes(x = k, y = logLik_mean),
                    color = "darkorange2", linewidth = 1.2)
        )
      } else NULL
    } +
    geom_vline(xintercept = best_k_overall, color = "black", linetype = "dashed", linewidth = 0.8) +
    scale_x_log10() +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(
      title = sprintf("Scenario: %s (%d replicates)", scenario_name, B),
      subtitle = subtitle,
      x = expression(k~("noise height, log scale")),
      y = expression( Final~log-likelihood )
    ) +
    theme_minimal()

  metrics_list <- lapply(replicate_results, function(r) r$metrics)
  metrics_combined <- rbindlist(metrics_list)

  non_numeric <- sapply(metrics_combined, function(x) !is.numeric(x) || all(is.na(x)))
  numeric_cols <- names(metrics_combined)[!non_numeric]

  if (length(numeric_cols) > 0) {
    metrics_summary <- metrics_combined[, c(
      as.list(colMeans(.SD[, ..numeric_cols], na.rm = TRUE)),
      as.list(sapply(.SD[, ..numeric_cols], function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))))
    ), .SDcols = numeric_cols]
    setnames(metrics_summary,
             c(numeric_cols, paste0(numeric_cols, "_se")),
             c(numeric_cols, paste0(numeric_cols, "_se")))
    for (col in names(metrics_combined)[non_numeric]) {
      metrics_summary[[col]] <- metrics_combined[[col]][1]
    }
  } else {
    metrics_summary <- metrics_combined[1, ]
  }

  list(
    plot_df   = agg,
    p_loglik  = p,
    best_k    = mean(sapply(replicate_results, function(r) r$best_k), na.rm = TRUE),
    k_grid    = replicate_results[[1]]$k_grid,
    ks_scores = replicate_results[[1]]$ks_scores,
    correlation = mean(sapply(replicate_results, function(r) r$correlation), na.rm = TRUE),
    metrics   = metrics_summary,
    B         = B,
    plots_dir = plots_dir
  )
}
