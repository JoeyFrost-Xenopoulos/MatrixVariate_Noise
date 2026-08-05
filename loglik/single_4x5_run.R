library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)
library(future)
library(future.apply)

plan(multisession, workers = 8)

source("loglik/run_all_base.R")

r <- 4
p <- 5
n_total <- 20
plots_dir <- file.path("loglik", "plots_base", sprintf("size_%dx%d", r, p))

cat(sprintf("\n===== Running base n20 | 4x54 | ONE run =====\n"))

set.seed(42)
sim_viroli <- scaled_viroli_simulation(
  r = r, p = p, n = n_total,
  n1 = round(0.3 * n_total), n2 = round(0.4 * n_total),
  n3 = n_total - round(0.3 * n_total) - round(0.4 * n_total),
  n_outliers = 5, signal_strength = 0.5, cov_scale = 1,
  outlier_type = "perm"
)

res <- run_base_replicates(
  x_list = sim_viroli$x_list, g = 3,
  scenario_name = "base_n20_4x54",
  subtitle = sprintf("Base Viroli 4x54 | 3 groups | 5 permuted outliers (n=20)"),
  r = r, p = p, true_k_noise = 5,
  outlier_idx = sim_viroli$outlier_idx,
  plots_dir = plots_dir,
  B = 1
)

if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
fwrite(res$metrics, file.path(plots_dir, "base_n20_4x54_metrics.csv"), row.names = FALSE)
cat(sprintf("Saved metrics to %s\n", file.path(plots_dir, "base_n20_4x54_metrics.csv")))
