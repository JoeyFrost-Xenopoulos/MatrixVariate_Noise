#!/usr/bin/env Rscript
# Run all noise test simulations (Tomarchio + Viroli)
# Usage: Rscript R/tests/run_all_noise_tests.R

# Resolve project root from script location so this works from any working directory
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[1])
} else {
  "run_all_noise_tests.R"
}
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))

find_project_root <- function(start_dir) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(6)) {
    if (dir.exists(file.path(current, "R"))) {
      return(current)
    }
    if (basename(current) == "R" && file.exists(file.path(dirname(current), "DESCRIPTION"))) {
      return(dirname(current))
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      break
    }
    current <- parent
  }
  NULL
}

project_root <- find_project_root(script_dir)
if (is.null(project_root)) {
  project_root <- if (dir.exists("R")) normalizePath(getwd(), winslash = "/", mustWork = TRUE) else NULL
}
if (is.null(project_root)) {
  stop("Could not locate project root containing R/ directory.")
}

setwd(project_root)

# Source package R files (except this runner) so helpers are available
r_files <- list.files(file.path(project_root, "R"), pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
runner_path <- normalizePath(file.path(project_root, "R", "tests", "run_all_noise_tests.R"), winslash = "/", mustWork = FALSE)
helpers <- setdiff(r_files, runner_path)
for (f in helpers) {
  tryCatch(source(f), error = function(e) stop("Error sourcing ", f, ": ", conditionMessage(e)))
}

# Ensure temp directory exists
temp_dir <- file.path(project_root, "temp")
if (!dir.exists(temp_dir)) dir.create(temp_dir)

cat("Starting all noise simulations...\n")

cat("1) Tomarchio — Clean (N_sim=50, n=100)\n")
tom_clean <- run_tomarchio_noise_simulation(N_sim = 50, n = 100, contam_n = 0, verbose = TRUE, methods = "hc")
saveRDS(tom_clean, file = file.path(temp_dir, "tomarchio_noise_clean.rds"))

cat("2) Tomarchio — Contaminated (N_sim=50, n=100, contam_n=30)\n")
tom_contam <- run_tomarchio_noise_simulation(N_sim = 10, n = 100, contam_n = 30, verbose = TRUE, methods = "hc")
saveRDS(tom_contam, file = file.path(temp_dir, "tomarchio_noise_contam.rds"))

cat("3) Viroli — Clean (N_sim=50, n=100)\n")
viro_clean <- run_viroli_noise_simulation(N_sim = 50, n = 100, contam_n = 0, verbose = TRUE, methods = "hc")
saveRDS(viro_clean, file = file.path(temp_dir, "viroli_noise_clean.rds"))

cat("4) Viroli — Contaminated (N_sim=50, n=100, contam_n=30)\n")
viro_contam <- run_viroli_noise_simulation(N_sim = 50, n = 100, contam_n = 30, verbose = TRUE, methods = "hc")
saveRDS(viro_contam, file = file.path(temp_dir, "viroli_noise_contam.rds"))

# Combine summaries
scenario_summary <- function(results, model_name, scenario_name) {
  df <- as.data.frame(results$summary, stringsAsFactors = FALSE)
  if ("sd_noise_rate" %in% names(df)) df$sd_noise_rate <- NULL
  df$method <- rownames(df)
  df$model <- model_name
  df$scenario <- scenario_name
  df <- df[, c("model", "scenario", "method", setdiff(names(df), c("model", "scenario", "method")))]
  rownames(df) <- NULL
  df
}

final_noise_summary <- rbind(
  scenario_summary(tom_clean, "Tomarchio", "Clean"),
  scenario_summary(tom_contam, "Tomarchio", "Contaminated"),
  scenario_summary(viro_clean, "Viroli", "Clean"),
  scenario_summary(viro_contam, "Viroli", "Contaminated")
)

print(final_noise_summary)
write.csv(final_noise_summary, file = file.path(temp_dir, "final_noise_summary.csv"), row.names = FALSE)

# Save all results
saveRDS(list(
  tomarchio_clean = tom_clean,
  tomarchio_contam = tom_contam,
  viroli_clean = viro_clean,
  viroli_contam = viro_contam,
  final_summary = final_noise_summary
), file = file.path(temp_dir, "all_noise_results.rds"))

cat("All simulations complete. Results saved to temp/\n")
