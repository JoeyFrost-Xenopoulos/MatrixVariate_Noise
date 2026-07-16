# Source all R files in order (dependencies first)
pkg_root <- normalizePath(file.path(dirname(dirname(getwd()))))
r_dir <- file.path(pkg_root, "R")

if (!dir.exists(r_dir)) {
  # Fallback: try relative to test file location
  r_dir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
}

source_order <- c(
  "matrix_utils.R",
  "matrix_base.R",
  "matrix_init_whiten.R",
  "matrix_init_kmeans.R",
  "matrix_init_dbscan.R",
  "matrix_init_emrefine.R",
  "matrix_noise_br.R",
  "ks_score.R",
  "matrix_parallel.R",
  "matrix_noise.R"
)

for (f in source_order) {
  fpath <- file.path(r_dir, f)
  if (file.exists(fpath)) source(fpath)
}
