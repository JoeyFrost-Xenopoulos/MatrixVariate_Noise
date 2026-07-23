# Run and save a single HC benchmark scenario.
#
# Usage:
#   devtools::load_all()
#   source("scripts/run_benchmark_estimate_k_hc.R")
#
# This script sources the benchmark helpers, runs one benchmark configuration,
# and writes CSV outputs to a timestamped folder under `results/`.

source(file.path("scripts", "benchmark_estimate_k_hc.R"))

benchmark_output_dir <- file.path(
  "results",
  paste0("hc_benchmark_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)

benchmark_k_grid <- mv_hc_benchmark_grid_for_dims(
  row_count = 2,
  col_count = 3,
  n_points = 120L
)

benchmark_run <- mv_hc_benchmark_run(
  n_per_group = c(10, 10),
  contamination_levels = c(0, 0.05, 0.10, 0.20),
  initializations = c("kmeans", "emrefine", "dbscan"),
  g = 2,
  replicates = 5,
  row_sd = 0.35,
  col_sd = 0.35,
  contamination_type = "column_replace",
  contamination_range = c(-15, 15),
  max_iter = 50,
  tol = 1e-06,
  nstart = 25,
  adaptive_grid = TRUE,
  k_grid = benchmark_k_grid,
  noise_pi_init = 0.05,
  use_parallel = FALSE,
  seed = 123,
  verbose = FALSE
)

saved_files <- mv_hc_benchmark_save_results(
  results = benchmark_run,
  output_dir = benchmark_output_dir,
  prefix = "hc_noise"
)

message("Benchmark finished. Results saved to: ", benchmark_output_dir)
message("Auto-selection CSV: ", saved_files$auto_selection)
if (!is.null(saved_files$exhaustive)) {
  message("Exhaustive CSV: ", saved_files$exhaustive)
}
if (!is.null(saved_files$best_k_ranges)) {
  message("Best-k ranges CSV: ", saved_files$best_k_ranges)
}
