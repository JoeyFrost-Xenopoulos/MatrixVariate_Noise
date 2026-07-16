## ---- benchmark: estimate_k robustness to contamination + noise_k search ----

## Goal:
##   Evaluate mv_noise_fit(estimate_k = TRUE) under varying contamination
##   levels and initialization methods, and find the noise_k at which the
##   algorithm consistently recovers the GROUND-TRUTH clustering.
##
##   For every (contamination_level, init) scenario we first let the built-in
##   KS grid SELECT k automatically, then POST-HOC check whether the returned
##   clustering matches the ground truth. If it does not, we walk a noise_k
##   ladder (coarse -> fine) calling mv_noise_fit with estimate_k = FALSE until
##   we find a noise_k that recovers the ground truth (or exhaust the ladder).
##
##   "Perfect recovery" = every non-contaminated matrix lands in its true
##   component AND every contaminated matrix is flagged as noise (cluster == 0),
##   up to a permutation of the component labels.
##
## Run (after devtools::load_all() in the package root, or install + library):
##   devtools::load_all()
##   source("scripts/benchmark_estimate_k_contamination.R")

## Ensure the Ampharos namespace is attached so bare `mv_noise_fit` resolves,
## whether it was brought in by devtools::load_all() or an installed library().
if (!requireNamespace("Ampharos", quietly = TRUE)) {
  stop("Ampharos is not available. Run devtools::load_all() in the package root",
       " or install.packages() / library(Ampharos) first.")
}
suppressPackageStartupMessages(library("Ampharos"))

set.seed(20260716)

## ============================================================ tunable params
r <- 2; p <- 3                       # matrix dimensions
g <- 2                               # number of true Gaussian components
n_clean_per_group <- 20             # clean matrices per true component
max_iter <- 40
nstart <- 20
init_methods <- c("kmeans", "dbscan")         # initialization strategies to sweep
                                  # (add "emrefine", "dbscan" if desired)
reps <- 4                           # replicate datasets per scenario

## contamination levels: fraction of matrices injected as outliers/noise
contamination_levels <- c(0, 0.05, 0.10, 0.20, 0.30, 0.40)

## noise_k ladder used for the post-hoc search when auto-selection fails.
## Coarse geometric sweep first, then a finer sweep around the first success.
noise_k_coarse <- 10^seq(-14, -1, length.out = 28)
noise_k_fine_template <- function(around) {
  10^seq(log10(around) - 1.5, log10(around) + 1.5, length.out = 20)
}

## ============================================================ data generation
## Two well-separated component means; clean matrices are perturbed by small SD.
## Contaminated matrices are drawn from a wide uniform cloud (clear outliers).
make_dataset <- function(seed, contam_frac, n_per_group) {
  set.seed(seed)
  m1 <- matrix(c(2, 1.8, 1.5, 1.7, 1.6, 1.9), r, p)
  m2 <- matrix(c(-2, -1.8, -1.5, -1.7, -1.6, -1.9), r, p)

  mk_clean <- function(n, m, sd = 0.3) {
    lapply(seq_len(n), function(i) m + matrix(rnorm(r * p, sd = sd), r, p))
  }
  clean1 <- mk_clean(n_per_group, m1)
  clean2 <- mk_clean(n_per_group, m2)

  n_contam <- round((n_per_group * 2) * contam_frac)
  contam <- if (n_contam > 0) {
    lapply(seq_len(n_contam), function(i)
      matrix(runif(r * p, min = -8, max = 8), r, p))
  } else list()

  ## interleave clean groups, then append contaminants at the end
  x_list <- c(clean1, clean2, contam)
  true_component <- c(rep(1, n_per_group), rep(2, n_per_group),
                      rep(NA_integer_, n_contam))   # NA -> ground-truth noise
  list(x_list = x_list, true_component = true_component)
}

## ============================================================ recovery check
## Relabel predicted component labels (0 = noise) to best match true component
## labels (NA = true noise), via greedy one-to-one assignment on the cross-tab.
relabel_to_true <- function(true_component, pred_cluster) {
  tc <- true_component
  pc <- as.integer(pred_cluster)
  true_levels <- sort(unique(tc[!is.na(tc)]))
  pred_levels <- sort(unique(pc[pc != 0]))

  used <- integer(0)
  mapping <- integer(length(pred_levels))
  names(mapping) <- as.character(pred_levels)

  for (t_lvl in true_levels) {
    best_par <- NA_integer_
    best_count <- -1L
    for (p_lvl in pred_levels) {
      if (p_lvl %in% used) next
      cnt <- sum(tc == t_lvl & pc == p_lvl, na.rm = TRUE)
      if (cnt > best_count) { best_count <- cnt; best_par <- p_lvl }
    }
    used <- c(used, best_par)
    mapping[as.character(best_par)] <- t_lvl
  }
  out <- pc
  for (p_lvl in pred_levels) {
    out[pc == p_lvl] <- mapping[as.character(p_lvl)]
  }
  out
}

## Perfect recovery: every non-noise true matrix in its true component, and
## every true-noise matrix flagged as cluster == 0.
is_perfect <- function(true_component, pred_cluster) {
  pc <- as.integer(pred_cluster)
  tc <- true_component
  ## true-noise must be classified as noise
  noise_ok <- all(pc[is.na(tc)] == 0)
  ## clean matrices must match their true component (after relabeling)
  rel <- relabel_to_true(tc, pc)
  clean_ok <- all(rel[!is.na(tc)] == tc[!is.na(tc)])
  noise_ok && clean_ok
}

## ============================================================ run one dataset
run_one <- function(x_list, true_component, init, contam_frac) {
  ## 1) automatic k selection (no k_grid supplied)
  auto_fit <- mv_noise_fit(
    x_list, g = g, noise_type = "hc",
    max_iter = max_iter, nstart = nstart,
    estimate_k = TRUE, k_grid = NULL,
    init = init, verbose = FALSE
  )
  auto_k <- auto_fit$k_selection$selected_k
  auto_ok <- is_perfect(true_component, auto_fit$cluster)

  ## 2) post-hoc search for a noise_k that recovers ground truth, if needed
  found_k <- if (auto_ok) auto_k else NA_real_
  if (!auto_ok) {
    for (k in noise_k_coarse) {
      fit <- mv_noise_fit(
        x_list, g = g, noise_type = "hc",
        max_iter = max_iter, nstart = nstart,
        estimate_k = FALSE, noise_k = k,
        init = init, verbose = FALSE
      )
      if (is_perfect(true_component, fit$cluster)) { found_k <- k; break }
    }
    ## finer sweep around first coarse success
    if (!is.na(found_k)) {
      fine <- noise_k_fine_template(found_k)
      for (k in fine) {
        fit <- mv_noise_fit(
          x_list, g = g, noise_type = "hc",
          max_iter = max_iter, nstart = nstart,
          estimate_k = FALSE, noise_k = k,
          init = init, verbose = FALSE
        )
        if (is_perfect(true_component, fit$cluster)) { found_k <- k; break }
      }
    }
  }

  list(
    auto_k = auto_k,
    auto_ok = auto_ok,
    found_k = found_k,
    cluster = auto_fit$cluster,
    logLik = auto_fit$logLik
  )
}

## ============================================================ main sweep
cat(sprintf(
  "Benchmark: r=%dx%d, g=%d, clean/group=%d, reps=%d, inits=[%s]\n",
  r, p, g, n_clean_per_group, reps, paste(init_methods, collapse = ", ")
))

results <- list()
details <- data.frame(
  init = character(0), contamination = numeric(0), rep = integer(0),
  auto_k = numeric(0), auto_ok = logical(0), recovered_k = numeric(0),
  stringsAsFactors = FALSE
)

for (init in init_methods) {
  for (cl in contamination_levels) {
    cat(sprintf("\n--- init=%s, contamination=%.2f ---\n", init, cl))
    auto_success <- 0L
    found_ks <- numeric(0)
    auto_ks <- numeric(0)

    for (rep_i in seq_len(reps)) {
      seed <- 1000 * as.integer(cl * 100) + rep_i
      ds <- make_dataset(seed, cl, n_clean_per_group)
      out <- run_one(ds$x_list, ds$true_component, init, cl)

      auto_ks <- c(auto_ks, out$auto_k)
      if (out$auto_ok) auto_success <- auto_success + 1L
      if (!is.na(out$found_k)) found_ks <- c(found_ks, out$found_k)

      details <- rbind(details, data.frame(
        init = init, contamination = cl, rep = rep_i,
        auto_k = out$auto_k, auto_ok = out$auto_ok,
        recovered_k = if (is.na(out$found_k)) NA_real_ else out$found_k,
        stringsAsFactors = FALSE
      ))

      cat(sprintf(
        "  rep %d: auto_k=%.3e auto_ok=%s recovered_k=%s\n",
        rep_i, out$auto_k,
        out$auto_ok,
        if (is.na(out$found_k)) "none" else sprintf("%.3e", out$found_k)
      ))
    }

    ## noise threshold for consistent recovery across replicates:
    ## if post-hoc search recovered in ALL reps, take the max needed noise_k
    ## (the worst-case requirement); if auto already perfect, threshold = auto_k.
    if (auto_success == reps) {
      threshold_k <- max(auto_ks)
      recovered <- TRUE
    } else if (length(found_ks) == reps) {
      threshold_k <- max(found_ks)
      recovered <- TRUE
    } else {
      threshold_k <- NA_real_
      recovered <- FALSE
    }

    results[[sprintf("%s_%.2f", init, cl)]] <- list(
      init = init,
      contamination = cl,
      auto_success_rate = auto_success / reps,
      auto_k_median = median(auto_ks),
      recovered_all = recovered,
      noise_threshold = threshold_k
    )

    cat(sprintf(
      "  >> auto success %.0f%%, recovery threshold noise_k = %s\n",
      auto_success / reps * 100,
      if (is.na(threshold_k)) "NOT ACHIEVABLE" else sprintf("%.3e", threshold_k)
    ))
  }
}

## ============================================================ summary table
cat("\n=== Summary: contamination vs noise_k threshold for consistent recovery ===\n")
summ <- data.frame(
  init = character(0), contamination = numeric(0),
  auto_success_rate = numeric(0), auto_k_median = numeric(0),
  recovered_all = logical(0), noise_threshold = numeric(0),
  stringsAsFactors = FALSE
)
for (nm in names(results)) {
  rr <- results[[nm]]
  summ <- rbind(summ, data.frame(
    init = rr$init,
    contamination = rr$contamination,
    auto_success_rate = rr$auto_success_rate,
    auto_k_median = rr$auto_k_median,
    recovered_all = rr$recovered_all,
    noise_threshold = rr$noise_threshold,
    stringsAsFactors = FALSE
  ))
}
print(summ, row.names = FALSE)

## ============================================================ write CSVs
## Save both the per-replicate detail and the per-scenario summary to CSV so the
## results can be inspected / plotted outside R.
out_dir <- file.path("scripts", "benchmark_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

detail_csv <- file.path(out_dir, "estimate_k_contamination_detail.csv")
summary_csv <- file.path(out_dir, "estimate_k_contamination_summary.csv")

utils::write.csv(details, detail_csv, row.names = FALSE)
utils::write.csv(summ, summary_csv, row.names = FALSE)

cat(sprintf("\nWrote detail  -> %s\n", detail_csv))
cat(sprintf("Wrote summary -> %s\n", summary_csv))

cat("\nDone.\n")
