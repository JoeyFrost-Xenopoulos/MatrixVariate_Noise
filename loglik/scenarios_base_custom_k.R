library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)
library(future)
library(future.apply)

plan(multisession, workers = 8)

source("loglik/run_all_base.R")

# ──────────────────────────────────────────────────────────────────────────────
# Custom k-grid
# ──────────────────────────────────────────────────────────────────────────────
custom_k_grid <- 10^seq(-100, -4, length.out = 96)

cat(sprintf("Using custom k_grid: [%.0e, %.0e] with %d candidates\n",
            min(custom_k_grid), max(custom_k_grid), length(custom_k_grid)))

# ──────────────────────────────────────────────────────────────────────────────
# Modified run_base_scenario that accepts a custom k_grid
# ──────────────────────────────────────────────────────────────────────────────
run_base_scenario_custom <- function(x_list, g, scenario_name,
                                      subtitle, r, p,
                                      true_k_noise = NA_integer_,
                                      outlier_idx = NULL,
                                      nstart = 10, max_iter = 100, tol = 1e-6,
                                      use_kmeans = TRUE, init = "kmeans",
                                      noise_pi_init = 0.05,
                                      save_plots = TRUE, plots_dir = "loglik/plots",
                                      seed = 42, k_grid = custom_k_grid) {
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
    k_grid       = k_grid,
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
    selected_idx <- match(best_k, plot_df$k)
    correct_noise_detected <- !is.na(selected_idx) &&
      isTRUE(plot_df$correct_noise[selected_idx])
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

  p1 <- ggplot() +
    {
      ok_rows <- plot_df[complete.cases(plot_df$logLik), ]
      if (nrow(ok_rows) > 1) {
        ok_rows <- ok_rows[order(ok_rows$k), ]
        seg_df <- data.frame(
          x      = head(ok_rows$k, -1),
          xend   = tail(ok_rows$k, -1),
          y      = head(ok_rows$logLik, -1),
          yend   = tail(ok_rows$logLik, -1),
          correct_noise = factor(
            head(ok_rows$correct_noise, -1),
            levels = c(FALSE, TRUE)
          )
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
      name = expression( Final~log-likelihood ),
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
# Run replicates of a base scenario with custom k_grid
# ──────────────────────────────────────────────────────────────────────────────
run_base_replicates_custom <- function(x_list, g, B = 10, base_seed = 42, k_grid = custom_k_grid, ...) {
  if (B < 2) {
    res <- run_base_scenario_custom(x_list = x_list, g = g, seed = base_seed,
                                     k_grid = k_grid, ...)
    res$B <- 1L
    return(res)
  }

  replicate_results <- future_lapply(seq_len(B), function(b) {
    run_base_scenario_custom(
      x_list = x_list,
      g = g,
      seed = base_seed + b * 1000,
      k_grid = k_grid,
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
    k_grid    = k_grid,
    ks_scores = replicate_results[[1]]$ks_scores,
    correlation = mean(sapply(replicate_results, function(r) r$correlation), na.rm = TRUE),
    metrics   = metrics_summary,
    B         = B,
    plots_dir = plots_dir
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Matrix size grid
# ──────────────────────────────────────────────────────────────────────────────
size_grid <- list(
  c(4, 4), c(4, 5), c(4, 6), c(5, 5), c(5, 6),
  c(5, 7), c(5, 8), c(6, 9), c(7, 10),
  c(8, 11), c(9, 12), c(10, 13)
)

base_plots_dir <- file.path("loglik", "plots_base_custom_k")

cat("\n===== Running base Viroli across size grid with custom k_grid =====\n")

results <- list()

for (s in seq_along(size_grid)) {
  rg <- size_grid[[s]][1]
  pg <- size_grid[[s]][2]
  n_total <- 150

  plots_dir <- file.path(base_plots_dir, sprintf("size_%dx%d", rg, pg))

  cat(sprintf("\n  >>> Size %d: r=%d, p=%d (n=%d) -> %s\n",
              s - 1, rg, pg, n_total, plots_dir))

  for (noise_prop in seq(0, 30, by = 1)) {
    prop <- noise_prop / 100
    n_outliers <- min(max(0, round(n_total * prop)), n_total - 3)

    set.seed(42 + 700 + s + noise_prop)
    sim_viroli_prop <- scaled_viroli_simulation(
      r = rg, p = pg, n = n_total,
      n1 = round(0.3 * n_total), n2 = round(0.4 * n_total),
      n3 = n_total - round(0.3 * n_total) - round(0.4 * n_total),
      n_outliers = n_outliers, signal_strength = 0.5, cov_scale = 1,
      outlier_type = "perm"
    )

    results[[sprintf("base_customk_prop%d_sz%d", noise_prop, s - 1)]] <- future(run_base_replicates_custom(
      x_list = sim_viroli_prop$x_list, g = 3,
      scenario_name = sprintf("base_customk_prop%d_sz%d", noise_prop, s - 1),
      subtitle = sprintf("Base Viroli %dx%d | 3 groups | %d permuted outliers (%d%% of n=%d)",
                         rg, pg, n_outliers, noise_prop, n_total),
      r = rg, p = pg, true_k_noise = n_outliers,
      outlier_idx = sim_viroli_prop$outlier_idx,
      plots_dir = plots_dir,
      B = 10
    ), seed = TRUE)
  }
}

cat("\n  Resolving all futures ...\n")
for (nm in names(results)) {
  results[[nm]] <- value(results[[nm]])
}

# ──────────────────────────────────────────────────────────────────────────────
# Save metrics
# ──────────────────────────────────────────────────────────────────────────────
all_metrics <- rbindlist(lapply(seq_along(results), function(i) {
  nm <- names(results)[[i]]
  r  <- results[[nm]]
  if (!is.null(r$metrics)) r$metrics[, scenario_idx := i]
}), fill = TRUE)

if (!dir.exists(base_plots_dir)) dir.create(base_plots_dir, recursive = TRUE)
fwrite(all_metrics, file.path(base_plots_dir, "summary_metrics.csv"), row.names = FALSE)
cat(sprintf("Saved base metrics to %s\n", file.path(base_plots_dir, "summary_metrics.csv")))

for (s in seq_along(size_grid)) {
  rg <- size_grid[[s]][1]; pg <- size_grid[[s]][2]
  sz_metrics <- all_metrics[grepl(sprintf("_sz%d$", s - 1), scenario)]
  if (nrow(sz_metrics) > 0) {
    plots_dir <- file.path(base_plots_dir, sprintf("size_%dx%d", rg, pg))
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
    fwrite(sz_metrics, file.path(plots_dir, "metrics.csv"), row.names = FALSE)
    cat(sprintf("Saved metrics for size %dx%d to %s\n", rg, pg, file.path(plots_dir, "metrics.csv")))
  }
}

# ──────────────────────────────────────────────────────────────────────────────
# Combined plots per size folder
# ──────────────────────────────────────────────────────────────────────────────
cat("\n===== Generating combined plots per size folder =====\n")

results_by_dir <- split(results, sapply(results, `[[`, "plots_dir"))

for (dir in names(results_by_dir)) {
  dir_results <- results_by_dir[[dir]]
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  combined_df <- rbindlist(lapply(dir_results, function(r) r$plot_df))

  has_se <- "logLik_se" %in% names(combined_df) && any(!is.na(combined_df$logLik_se))

  p_combined <- ggplot(combined_df, aes(x = k, y = ifelse(has_se, logLik_mean, logLik))) +
    {
      if (has_se) {
        geom_ribbon(aes(ymin = logLik_mean - logLik_se, ymax = logLik_mean + logLik_se),
                    fill = "darkorange2", alpha = 0.2)
      } else NULL
    } +
    geom_line(color = "darkorange2", linewidth = 0.8) +
    facet_wrap(~scenario_name, scales = "free_y") +
    scale_x_log10() +
    scale_y_continuous(labels = scales::comma_format()) +
    labs(
      title = sprintf("Base Viroli: combined log-likelihood trends (%s)", basename(dir)),
      x = expression(k~("noise height, log scale")),
      y = ifelse(has_se, expression(Mean~Final~log-likelihood), expression(Final~log-likelihood))
    ) +
    theme_minimal() +
    theme(strip.text.x = element_text(size = 7))

  combined_path <- file.path(dir, "combined.png")
  ggsave(filename = combined_path, plot = p_combined, width = 16, height = 12, dpi = 150)
  cat(sprintf("Saved combined plot to %s\n", combined_path))
}

cat(sprintf("\nDone. %d scenarios saved under %s/\n", length(results), normalizePath(base_plots_dir)))
