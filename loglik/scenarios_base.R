library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)
library(future)
library(future.apply)

plan(multisession, workers = 8)

source("loglik/run_all_base.R")

# ──────────────────────────────────────────────────────────────────────────────
# Matrix size grid
# ──────────────────────────────────────────────────────────────────────────────
size_grid <- list(
  c(2, 3), c(2, 4), c(3, 3), c(3, 4), c(3, 5),
  c(4, 6), c(5, 7), c(5, 8), c(6, 9), c(7, 10),
  c(8, 11), c(9, 12), c(10, 13)
)

base_plots_dir <- file.path("loglik", "plots_base")

cat("\n===== Running base Viroli across size grid with increasing contamination =====\n")

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

    results[[sprintf("base_prop%d_sz%d", noise_prop, s - 1)]] <- future(run_base_replicates(
      x_list = sim_viroli_prop$x_list, g = 3,
      scenario_name = sprintf("base_prop%d_sz%d", noise_prop, s - 1),
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
