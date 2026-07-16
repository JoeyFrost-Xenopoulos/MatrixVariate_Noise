## ---- benchmark: mv_noise_fit n_cores scaling ----

## Goal: time mv_noise_fit across different n_cores settings using
## kmeans initialization with automatic k selection (estimate_k = TRUE,
## no k_grid supplied -> adaptive heuristic grid + nstart restarts).

## Run (after load_all in the package root, or an installed library):
##   devtools::load_all()
##   source("scripts/benchmark_ncores_mv_noise_fit.R")

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
})

## Ensure the Ampharos namespace is available (via devtools::load_all() or an
## installed library() call) without forcing a double-load.
if (!requireNamespace("Ampharos", quietly = TRUE)) {
  stop("Ampharos is not available. Run devtools::load_all() in the package root",
       " or install.packages() / library(Ampharos) first.")
}

set.seed(20260716)

## ---------------------------------------------------------------- data setup
## Two well-separated 2x3 matrix groups so k selection / mixture fitting
## has real work to do. Scale n_each to make timing meaningful.
make_two_group <- function(seed, n_each = 40, sd = 0.3) {
  set.seed(seed)
  r <- 2; p <- 3
  m1 <- matrix(c(2, 1.8, 1.5, 1.7, 1.6, 1.9), r, p)
  m2 <- matrix(c(-2, -1.8, -1.5, -1.7, -1.6, -1.9), r, p)
  mk <- function(n, m) lapply(seq_len(n), function(i)
    m + matrix(rnorm(r * p, sd = sd), r, p))
  c(mk(n_each, m1), mk(n_each, m2))
}

x_list <- make_two_group(1, n_each = 40)

g <- 2
max_cores <- future::availableCores()

## n_cores configurations to evaluate:
##   - "sequential": use_parallel = FALSE (single core, no future overhead)
##   - "default"   : use_parallel = TRUE, n_cores = NULL (optimal default)
##   - 1..max_cores: explicit core counts
core_configs <- c("sequential", "default", as.character(seq_len(max_cores)))

cat(sprintf("\nData: %d matrices, %d cores available.\n", length(x_list), max_cores))
cat("Configs:", paste(core_configs, collapse = ", "), "\n\n")

## ---------------------------------------------------------------- timing loop
run_once <- function(use_parallel, n_cores) {
  t0 <- system.time({
    fit <- mv_noise_fit(
      x_list, g = g, noise_type = "hc",
      max_iter = 50, nstart = 25,
      estimate_k = TRUE, k_grid = NULL,   # adaptive heuristic grid
      init = "kmeans",
      verbose = FALSE,
      use_parallel = use_parallel,
      n_cores = n_cores
    )
  })
  list(
    elapsed = t0[["elapsed"]],
    selected_k = fit$k_selection$selected_k,
    logLik = fit$logLik
  )
}

results <- data.frame(
  config = character(0),
  n_cores = integer(0),
  elapsed = numeric(0),
  selected_k = numeric(0),
  stringsAsFactors = FALSE
)

reps <- 3
for (cfg in core_configs) {
  if (cfg == "sequential") {
    up <- FALSE; nc <- NULL
  } else if (cfg == "default") {
    up <- TRUE; nc <- NULL
  } else {
    up <- TRUE; nc <- as.integer(cfg)
  }

  cat(sprintf("Running config '%s' (use_parallel=%s, n_cores=%s) x%d ... ",
              cfg, up, if (is.null(nc)) "NULL" else nc, reps))
  times <- numeric(reps)
  sel_k <- NA
  for (i in seq_len(reps)) {
    out <- run_once(up, nc)
    times[i] <- out$elapsed
    sel_k <- out$selected_k
  }
  med <- median(times)
  cat(sprintf("median %.2fs (reps: %s)\n",
              med, paste(sprintf("%.2f", times), collapse = ", ")))

  results <- rbind(results, data.frame(
    config = cfg,
    n_cores = if (is.null(nc)) NA_integer_ else nc,
    elapsed = med,
    selected_k = sel_k,
    stringsAsFactors = FALSE
  ))
}

## ---------------------------------------------------------------- summary
cat("\n=== Timing summary (median elapsed seconds) ===\n")
print(results, row.names = FALSE)

seq_time <- results$elapsed[results$config == "sequential"]
if (length(seq_time) == 1 && seq_time > 0) {
  results$speedup <- seq_time / results$elapsed
  cat("\n=== Speedup vs sequential ===\n")
  print(results[, c("config", "n_cores", "elapsed", "speedup")], row.names = FALSE)
}

## ---------------------------------------------------------------- residuals
plan(sequential)
invisible(gc())
cat("\nDone.\n")
