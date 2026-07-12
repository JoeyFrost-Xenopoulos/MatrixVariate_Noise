# Source all R files in order (dependencies first)
pkg_root <- normalizePath(file.path(dirname(dirname(getwd()))))
r_dir <- file.path(pkg_root, "R")

if (!dir.exists(r_dir)) {
  # Fallback: try relative to test file location
  r_dir <- normalizePath(file.path("..", "..", "R"), mustWork = FALSE)
}

source_order <- c(
  "Utils.R",
  "Matrix_Base.R",
  "Init_Whiten.R",
  "Init_Kmeans.R",
  "Init_DBSCAN.R",
  "Init_EMRefine.R",
  "Matrix_Noise_BR.R",
  "KS_Score.R",
  "Matrix_Noise.R"
)

for (f in source_order) {
  fpath <- file.path(r_dir, f)
  if (file.exists(fpath)) source(fpath)
}
