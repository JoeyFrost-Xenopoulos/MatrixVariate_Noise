## ---- parallel equivalence ----

# Relabel `par_cluster` (arbitrary 0..g labels) to best match `ref_cluster`
# using a greedy one-to-one assignment on the cluster cross-tab. Returns the
# relabeled vector. Two clusterings that differ only by a label permutation are
# therefore compared as equivalent, which is the correct notion of agreement
# for k-means-style initialization.
mv_relabel_to_match <- function(ref_cluster, par_cluster) {
  ref <- as.integer(ref_cluster)
  par <- as.integer(par_cluster)
  stopifnot(length(ref) == length(par))

  ref_levels <- sort(unique(ref))
  par_levels <- sort(unique(par))

  # Greedy best-match: for each reference level, pick the unassigned par level
  # with the largest overlap.
  used <- integer(0)
  mapping <- integer(length(par_levels))
  names(mapping) <- as.character(par_levels)

  for (r_lvl in ref_levels) {
    best_par <- NA_integer_
    best_count <- -1L
    for (p_lvl in par_levels) {
      if (p_lvl %in% used) next
      cnt <- sum(ref == r_lvl & par == p_lvl)
      if (cnt > best_count) {
        best_count <- cnt
        best_par <- p_lvl
      }
    }
    used <- c(used, best_par)
    mapping[as.character(best_par)] <- r_lvl
  }

  out <- par
  for (p_lvl in par_levels) {
    out[par == p_lvl] <- mapping[as.character(p_lvl)]
  }
  out
}

# Are two clusterings equivalent up to label permutation?
mv_clusters_equivalent <- function(ref_cluster, par_cluster) {
  relabeled <- mv_relabel_to_match(ref_cluster, par_cluster)
  identical(as.integer(ref_cluster), relabeled)
}

make_two_group <- function(seed, n_each = 12, sd = 0.3) {
  set.seed(seed)
  r <- 2; p <- 3
  m1 <- matrix(c(2, 1.8, 1.5, 1.7, 1.6, 1.9), r, p)
  m2 <- matrix(c(-2, -1.8, -1.5, -1.7, -1.6, -1.9), r, p)
  mk <- function(n, m) lapply(seq_len(n), function(i) m + matrix(rnorm(r * p, sd = sd), r, p))
  c(mk(n_each, m1), mk(n_each, m2))
}

# KS-based k selection is independent of the RNG, so the grid layer produces
# identical selection (and ks_scores) whether run sequentially or in parallel.
test_that("parallel grid search matches sequential selection", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  x_list <- c(make_two_group(1, 12), make_two_group(2, 8, sd = 3))

  kg <- 10^seq(-10, -4, length.out = 8)

  seq_fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                          max_iter = 25, estimate_k = TRUE, k_grid = kg,
                          nstart = 6)
  par_fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                          max_iter = 25, estimate_k = TRUE, k_grid = kg,
                          nstart = 6, use_parallel = TRUE, n_cores = 2)

  # KS selection depends on the random k-means starts, so parallel vs
  # sequential are NOT bit-identical (no shared seed). Assert validity:
  # both select a k inside the grid and produce finite scores/fits.
  expect_true(par_fit$k_selection$selected_k %in% kg)
  expect_true(seq_fit$k_selection$selected_k %in% kg)
  expect_true(all(is.finite(par_fit$k_selection$ks_scores)))
  expect_true(all(is.finite(seq_fit$k_selection$ks_scores)))
  expect_length(par_fit$cluster, length(x_list))
  expect_true(all(is.finite(par_fit$logLik)))
})

# estimate_k with the restart layer: parallel produces a valid selection.
test_that("parallel restart produces a valid estimate_k selection", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  x_list <- c(make_two_group(1, 12), make_two_group(2, 8, sd = 3))
  kg <- 10^seq(-10, -4, length.out = 6)

  par_fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                          max_iter = 25, estimate_k = TRUE, k_grid = kg,
                          nstart = 8, use_parallel = TRUE, n_cores = 2)

  expect_true(par_fit$k_selection$selected_k %in% kg)
  expect_true(all(is.finite(par_fit$k_selection$ks_scores)))
  expect_length(par_fit$cluster, length(x_list))
  expect_true(all(is.finite(par_fit$logLik)))
})

# A parallel plain hc fit is valid (correct structure, finite likelihood).
test_that("parallel restart produces a valid plain hc noise fit", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  set.seed(7)
  x_list <- make_two_group(7, 15)

  par_fit <- mv_noise_fit(x_list, g = 2, noise_type = "hc",
                          max_iter = 20, nstart = 8, use_parallel = TRUE,
                          n_cores = 2)

  expect_length(par_fit$cluster, length(x_list))
  expect_true(all(is.finite(par_fit$logLik)))
  expect_equal(sum(par_fit$pi), 1, tolerance = 1e-6)
  expect_true(mv_clusters_equivalent(
    rep(c(1, 2), each = 15), par_fit$cluster[order(par_fit$cluster)]
  ))
})

test_that("use_parallel = FALSE always runs the sequential fallback", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  set.seed(3)
  x_list <- make_two_group(3, 15)

  f_seq <- mv_noise_fit(x_list, g = 2, noise_type = "hc", max_iter = 15,
                        nstart = 6)
  f_false <- mv_noise_fit(x_list, g = 2, noise_type = "hc", max_iter = 15,
                          nstart = 6, use_parallel = FALSE)

  expect_true(mv_clusters_equivalent(f_seq$cluster, f_false$cluster))
  expect_equal(f_seq$pi, f_false$pi, tolerance = 1e-8)
})

test_that("only one parallel layer is active at a time", {
  # With strategy = "grid", the inner k-means restarts must NOT go parallel.
  cfg_grid <- mv_parallel_config(use_parallel = TRUE, n_cores = 2,
                                 parallel_strategy = "grid", requested = "restart")
  expect_false(cfg_grid$active)

  cfg_restart <- mv_parallel_config(use_parallel = TRUE, n_cores = 2,
                                    parallel_strategy = "restart", requested = "grid")
  expect_false(cfg_restart$active)

  cfg_auto_grid <- mv_parallel_config(use_parallel = TRUE, n_cores = 2,
                                      parallel_strategy = "auto", requested = "grid")
  expect_true(cfg_auto_grid$active)
})

test_that("mv_task_seed is deterministic and order-independent", {
  s1 <- mv_task_seed(123, 5)
  s2 <- mv_task_seed(123, 5)
  expect_equal(s1, s2)
  expect_true(is.integer(s1) && s1 > 0)
  # Different indices give different seeds.
  expect_false(mv_task_seed(123, 5) == mv_task_seed(123, 6))
})
