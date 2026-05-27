#!/usr/bin/env Rscript
# Run all noise test simulations (Tomarchio + Viroli)
# Usage: Rscript R/run_all_noise_tests.R

# Ensure we're in project root where R/ exists
if (!dir.exists("R")) stop("R/ directory not found. Run this script from the project root.")

# Source package R files (except this runner) so helpers are available
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
runner_path <- normalizePath("R/run_all_noise_tests.R", winslash = "/", mustWork = FALSE)
helpers <- setdiff(r_files, runner_path)
for (f in helpers) {
  tryCatch(source(f), error = function(e) stop("Error sourcing ", f, ": ", conditionMessage(e)))
}

# Ensure temp directory exists
if (!dir.exists("temp")) dir.create("temp")

cat("Starting all noise simulations...\n")

cat("1) Tomarchio — Clean (N_sim=50, n=100)\n")
tom_clean <- run_tomarchio_noise_simulation(N_sim = 50, n = 100, contam_n = 0, verbose = TRUE)
saveRDS(tom_clean, file = "temp/tomarchio_noise_clean.rds")

cat("2) Tomarchio — Contaminated (N_sim=10, n=100, contam_n=30)\n")
tom_contam <- run_tomarchio_noise_simulation(N_sim = 10, n = 100, contam_n = 30, verbose = TRUE)
saveRDS(tom_contam, file = "temp/tomarchio_noise_contam.rds")

cat("3) Viroli — Clean (N_sim=10, n=100)\n")
viro_clean <- run_viroli_noise_simulation(N_sim = 10, n = 100, contam_n = 0, verbose = TRUE)
saveRDS(viro_clean, file = "temp/viroli_noise_clean.rds")

cat("4) Viroli — Contaminated (N_sim=10, n=30, contam_n=15)\n")
viro_contam <- run_viroli_noise_simulation(N_sim = 10, n = 30, contam_n = 15, verbose = TRUE)
saveRDS(viro_contam, file = "temp/viroli_noise_contam.rds")

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
write.csv(final_noise_summary, file = "temp/final_noise_summary.csv", row.names = FALSE)

# Save all results
saveRDS(list(
  tomarchio_clean = tom_clean,
  tomarchio_contam = tom_contam,
  viroli_clean = viro_clean,
  viroli_contam = viro_contam,
  final_summary = final_noise_summary
), file = "temp/all_noise_results.rds")

cat("All simulations complete. Results saved to temp/\n")
