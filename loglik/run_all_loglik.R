library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)
library(future)
library(future.apply)

plan(multisession, workers = availableCores() - 1,
     packages = c("Ampharos", "ggplot2", "clusterGeneration", "data.table"))

# ──────────────────────────────────────────────────────────────────────────────
# Shared helper — evaluate one k-grid candidate and return log-likelihood
# ──────────────────────────────────────────────────────────────────────────────
evaluate_k_candidate_with_loglik <- function(idx, x_list, r, p, g,
                                              max_iter, tol,
                                              k_grid, nstart, init,
                                              noise_pi_init, verbose) {
  current_k <- k_grid[idx]

  fit_noise <-  mv_noise_fit_impl(
    x_list        = x_list,
    g             = g,
    noise_type    = "hc",
    max_iter      = max_iter,
    tol           = tol,
    nstart        = nstart,
    noise_k       = current_k,
    noise_pi_init = noise_pi_init,
    init          = init,
    verbose       = FALSE
  )

  keep_idx <- fit_noise$cluster != 0
  x_clean  <- x_list[keep_idx]

  if (length(x_clean) <= g) {
    return(list(
      k            = current_k,
      logLik       = NA_real_,
      ks_statistic = Inf,
      ks_p.value   = NA_real_,
      n_used       = length(x_clean)
    ))
  }

  fit_clean <- tryCatch(
    mv_mixture_fit(x_list = x_clean, g = g,
                   max_iter = max_iter, tol = tol, verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit_clean)) {
    return(list(
      k            = current_k,
      logLik       = NA_real_,
      ks_statistic = Inf,
      ks_p.value   = NA_real_,
      n_used       = length(x_clean)
    ))
  }

  ks_result <- suppressWarnings(
    tryCatch(
      mv_noise_ks_score(fit_clean, x_clean),
      error = function(e) list(statistic = Inf, p.value = NA_real_,
                                n_used = length(x_clean))
    )
  )

  list(
    k            = current_k,
    logLik       = tail(fit_clean$logLik, 1),
    ks_statistic = ks_result$statistic,
    ks_p.value   = ks_result$p.value,
    n_used       = ks_result$n_used
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# Scenario helpers
# ──────────────────────────────────────────────────────────────────────────────
viroli_simulation <- function(r, p, n, n1, n2, n3, M1, M2, M3,
                               U1, V1, U2, V2, U3, V3, n_outliers,
                               row_cov_scale = 1, col_cov_scale = 1) {
  U1 <- row_cov_scale * U1; V1 <- col_cov_scale * V1
  U2 <- row_cov_scale * U2; V2 <- col_cov_scale * V2
  U3 <- row_cov_scale * U3; V3 <- col_cov_scale * V3
  simulate_matrix_group <- function(n, mean_mat, row_cov, col_cov) {
    lapply(seq_len(n), function(i) {
      mean_mat + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov
    })
  }
  x_list <- c(
    simulate_matrix_group(n1, M1, U1, V1),
    simulate_matrix_group(n2, M2, U2, V2),
    simulate_matrix_group(n3, M3, U3, V3)
  )
  outlier_idx <- sample(length(x_list), n_outliers)
  for (idx in outlier_idx) {
    x_list[[idx]] <- matrix(sample(x_list[[idx]]), r, p)
  }
  list(x_list = x_list, outlier_idx = outlier_idx)
}

two_group_uniform_simulation <- function(r, p, n1, n2, M1, M2,
                                         n_outliers, n_outliers_perm,
                                         n_row_spike = 0, n_col_out = 0,
                                         col_out_range = c(-15, 15)) {
  simulate_group <- function(n, mean_matrix, row_sd = 0.5, col_sd = 0.5) {
    row_cov <- diag(row_sd, nrow(mean_matrix))
    col_cov <- diag(col_sd, ncol(mean_matrix))
    lapply(seq_len(n), function(i)
      mean_matrix + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov)
  }
  x_list <- c(simulate_group(n1, M1), simulate_group(n2, M2))

  if (n_outliers > 0) {
    contam <- lapply(seq_len(n_outliers),
                     function(i) matrix(runif(r * p, -3, 3), r, p))
    x_list <- c(x_list, contam)
  }
  if (n_outliers_perm > 0) {
    perm_idx <- sample(seq_along(x_list), n_outliers_perm)
    for (idx in perm_idx) {
      x_list[[idx]] <- matrix(sample(x_list[[idx]]), r, p)
    }
  }
  if (n_row_spike > 0) {
    spike_idx <- sample(seq_along(x_list), n_row_spike)
    for (idx in spike_idx) {
      row_id <- sample.int(r, 1)
      x_list[[idx]][row_id, ] <- runif(p, -8, 8)
    }
  }
  if (n_col_out > 0) {
    col_idx <- sample(seq_along(x_list), n_col_out)
    for (idx in col_idx) {
      col_id <- sample.int(p, 1)
      x_list[[idx]][, col_id] <- runif(r, col_out_range[1], col_out_range[2])
    }
  }
  list(x_list = x_list)
}

tomarchio_simulation <- function(n, n_outliers = 10) {
  r <- 2; p <- 4
  M1 <- matrix(c(-2.60, -1.10, -0.50, -0.20, 1.30, 0.60, 0.30, 0.10),
               nrow = r, ncol = p, byrow = FALSE)
  M2 <- matrix(c(1.50, 1.70, 1.90, 2.20, -3.70, -2.70, -2.00, -1.50),
               nrow = r, ncol = p, byrow = FALSE)
  U1 <- matrix(c(2, 0, 0, 1), nrow = r, ncol = r)
  U2 <- matrix(c(1.70, 0.5, 0.5, 1.30), nrow = r, ncol = r)
  V_common <- matrix(c(1.00, 0.50, 0.25, 0.13,
                       0.50, 1.00, 0.50, 0.25,
                       0.25, 0.50, 1.00, 0.50,
                       0.13, 0.25, 0.50, 1.00), nrow = p, ncol = p)

  n1 <- round(0.5 * n); n2 <- n - n1

  simulate_tomarchio_group <- function(n, mean_mat, row_cov, col_cov) {
    lapply(seq_len(n), function(i) {
      mean_mat + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov
    })
  }

  x_list <- c(
    simulate_tomarchio_group(n1, M1, U1, V_common),
    simulate_tomarchio_group(n2, M2, U2, V_common)
  )

  if (n_outliers > 0) {
    out_idx <- sample(seq_along(x_list), n_outliers)
    for (idx in out_idx) {
      col_to_replace <- sample.int(p, 1)
      x_list[[idx]][, col_to_replace] <- runif(r, -15, 15)
    }
  }
  list(x_list = x_list)
}

# Scaled simulation helpers for matrix size scaling

scaled_viroli_simulation <- function(r, p, n, n1, n2, n3, n_outliers,
                                     signal_strength = 0.5, cov_scale = 1,
                                     outlier_type = c("perm", "row_spike", "col_out"),
                                     seed = NULL) {
  outlier_type <- match.arg(outlier_type)
  if (!is.null(seed)) set.seed(seed)
  M1 <- matrix(0, r, p)
  M2 <- matrix(0, r, p)
  M3 <- matrix(0, r, p)
  signal_col <- min(2, p)
  for (j in seq_len(signal_col)) {
    M1[1, j] <- signal_strength
    M3[1, j] <- -signal_strength
  }
  U1 <- cov_scale * rcorrmatrix(r)
  U2 <- cov_scale * rcorrmatrix(r)
  U3 <- cov_scale * rcorrmatrix(r)
  V1 <- cov_scale * rcorrmatrix(p)
  V2 <- cov_scale * rcorrmatrix(p)
  V3 <- cov_scale * rcorrmatrix(p)
  result <- viroli_simulation(r, p, n, n1, n2, n3,
                              M1, M2, M3, U1, V1, U2, V2, U3, V3,
                              n_outliers = n_outliers)
  if (outlier_type == "row_spike") {
    spike_idx <- sample(seq_along(result$x_list), n_outliers)
    for (idx in spike_idx) {
      row_id <- sample.int(r, 1)
      result$x_list[[idx]][row_id, ] <- runif(p, -8, 8)
    }
  }
  if (outlier_type == "col_out") {
    col_out_idx <- sample(seq_along(result$x_list), n_outliers)
    for (idx in col_out_idx) {
      col_id <- sample.int(p, 1)
      result$x_list[[idx]][, col_id] <- runif(r, -10, 10)
    }
  }
  list(x_list = result$x_list, outlier_idx = result$outlier_idx)
}

scaled_two_group_simulation <- function(r, p, n1, n2, n_outliers = 0,
                                        n_outliers_perm = 0,
                                        n_row_spike = 0, n_col_out = 0,
                                        mean_diff = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  M1 <- matrix(mean_diff / 2, r, p)
  M2 <- matrix(-mean_diff / 2, r, p)
  result <- two_group_uniform_simulation(
    r = r, p = p, n1 = n1, n2 = n2,
    M1 = M1, M2 = M2,
    n_outliers = n_outliers,
    n_outliers_perm = n_outliers_perm,
    n_row_spike = n_row_spike,
    n_col_out = n_col_out
  )
  list(x_list = result$x_list)
}

scaled_tomarchio_simulation <- function(r, p, n, n_outliers = 10, seed = NULL) {
  stopifnot(r >= 2, p >= 4)
  if (!is.null(seed)) set.seed(seed)
  cols_pattern <- c(-2.60, -1.10, -0.50, -0.20,  1.30, 0.60, 0.30, 0.10)
  col2_pattern <- c( 1.50,  1.70,  1.90,  2.20, -3.70,-2.70,-2.00,-1.50)
  M1 <- matrix(0, r, p)
  M2 <- matrix(0, r, p)
  for (j in seq_len(min(p, 4))) {
    idx <- (j - 1) * min(r, 2) + seq_len(min(r, 2))
    if (length(idx) == 1) {
      M1[1, j] <- cols_pattern[idx]
      M2[1, j] <- col2_pattern[idx]
    } else {
      M1[seq_len(min(r, 2)), j] <- cols_pattern[idx]
      M2[seq_len(min(r, 2)), j] <- col2_pattern[idx]
    }
  }
  if (p > 4 && r >= 2) {
    M1[1, 5] <- 0.6; M1[2, 5] <- 1.1
    M2[1, 5] <- -1.4; M2[2, 5] <- -1.8
  }
  U1 <- matrix(2, r, r)
  U2 <- matrix(1.7, r, r)
  diag(U1) <- c(2, rep(1, r - 1))
  diag(U2) <- c(1.7, rep(1.3, r - 1))
  idx <- row(U1); jdx <- col(U1)
  U1[idx != jdx] <- 0.5
  U2[idx != jdx] <- 0.5
  V <- matrix(0, p, p)
  if (p >= 4) {
    V0 <- matrix(c(1.00, 0.50, 0.25, 0.13,
                   0.50, 1.00, 0.50, 0.25,
                   0.25, 0.50, 1.00, 0.50,
                   0.13, 0.25, 0.50, 1.00), nrow = 4, ncol = 4)
    V[seq(4), seq(4)] <- V0
  }
  diag(V) <- 1
  V0 <- V
  n1 <- round(0.5 * n); n2 <- n - n1
  simulate_tomarchio_group <- function(n, mean_mat, row_cov, col_cov) {
    lapply(seq_len(n), function(i) {
      mean_mat + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov
    })
  }
  x_list <- c(
    simulate_tomarchio_group(n1, M1, U1, V),
    simulate_tomarchio_group(n2, M2, U2, V)
  )
  if (n_outliers > 0) {
    out_idx <- sample(seq_along(x_list), n_outliers)
    for (idx in out_idx) {
      col_to_replace <- sample.int(p, 1)
      x_list[[idx]][, col_to_replace] <- runif(r, -15, 15)
    }
  }
  list(x_list = x_list)
}
# ──────────────────────────────────────────────────────────────────────────────
run_loglik_scenario <- function(x_list, g, scenario_name,
                                 subtitle, r, p,
                                 true_k_noise = NA_integer_,
                                 nstart = 10, max_iter = 100, tol = 1e-6,
                                 use_kmeans = TRUE, init = "kmeans",
                                 noise_pi_init = 0.05,
                                 save_plots = TRUE, plots_dir = "loglik/plots") {
  cat(sprintf("\n===== Scenario: %s =====\n", scenario_name))

  cat(sprintf("  Fitting estimate_k with g = %d ...\n", g))
  fit_noise <- mv_noise_fit(
    x_list       = x_list,
    g            = g,
    noise_type   = "hc",
    nstart       = nstart,
    estimate_k   = TRUE,
    init         = init,
    verbose      = TRUE,
    use_parallel = FALSE
  )

   k_selection <- fit_noise$k_selection
  best_k      <- k_selection$selected_k
  best_ks_idx <- which.min(k_selection$ks_scores)

  cat(sprintf("  KS-selected k = %.4e (KS = %.4f)\n",
              best_k, k_selection$ks_scores[best_ks_idx]))

  cat("  Collecting log-likelihoods across k_grid ...\n")
  grid_results <- lapply(
    seq_along(k_selection$k_grid),
    evaluate_k_candidate_with_loglik,
    x_list        = x_list,
    r             = r,
    p             = p,
    g             = g,
    max_iter      = max_iter,
    tol           = tol,
    k_grid        = k_selection$k_grid,
    nstart        = nstart,
    init          = init,
    noise_pi_init = noise_pi_init,
    verbose       = FALSE
  )

  plot_df <- rbindlist(grid_results)

  ok_loglik <- complete.cases(plot_df$logLik)
  if (any(ok_loglik)) {
    max_loglik_k <- plot_df$k[which.max(plot_df$logLik[ok_loglik])]
    plot_df$type  <- rep("Other", nrow(plot_df))
    plot_df$type[plot_df$k == max_loglik_k] <- "Max log-likelihood"
    plot_df$type[plot_df$k == best_k]        <- "Selected k"
  } else {
    plot_df$type <- rep("Other", nrow(plot_df))
    plot_df$type[plot_df$k == best_k] <- "Selected k"
  }

  plot_df$ks_flag <- plot_df$k == k_selection$k_grid[best_ks_idx]
  correct_noise_detected <- !is.na(true_k_noise) && best_k == true_k_noise

   selected_k_color <- if (correct_noise_detected) "green" else "black"

   if (!is.na(true_k_noise)) {
     cat(sprintf("  True noise count = %.4e | Correctly identified = %s\n",
                 true_k_noise, if (correct_noise_detected) "YES" else "NO"))
   }

   # Diagnostics
   if (any(plot_df$type == "Max log-likelihood")) {
     cat(sprintf("  Max log-likelihood k = %.4e\n",
                 plot_df$k[plot_df$type == "Max log-likelihood"][1]))
   } else {
     cat("  Max log-likelihood k = (no valid values)\n")
   }

   cor_val <- NA_real_
  if (any(ok_loglik)) {
    cor_val <- cor(plot_df$k[ok_loglik], plot_df$logLik[ok_loglik],
                   method = "spearman")
    cat(sprintf("  Spearman (k vs logLik): %.4f\n", cor_val))
  }

   # ── Plot 1: log-likelihood vs k (linear y) ─────────────────────────────────
   p1 <- ggplot(plot_df, aes(x = k, y = logLik)) +
     geom_vline(xintercept = best_k, color = selected_k_color, linetype = "dashed",
                linewidth = 0.8) +
     geom_vline(xintercept = k_selection$k_grid[best_ks_idx],
                color = "steelblue", linetype = "dotdash", linewidth = 0.8)
   if (!is.na(true_k_noise)) {
     p1 <- p1 + geom_vline(xintercept = true_k_noise, color = "red", linetype = "twodash", linewidth = 0.8)
   }
   p1 <- p1 +
     geom_line(color = "darkorange2", linewidth = 1) +
     geom_point(data = subset(plot_df, type != "Other"),
                aes(shape = type), color = "darkorange2", size = 3) +
     scale_x_log10() +
     scale_y_continuous(labels = scales::comma_format()) +
     labs(
       title    = sprintf("Scenario: %s", scenario_name),
       subtitle = subtitle,
       x        = expression(k~("noise height, log scale")),
       y        = expression( Final~log-likelihood ),
       shape    = "",
       colour   = ""
     ) +
     theme_minimal() +
     theme(legend.position = "bottom", legend.title = element_blank())

  # ── Plot 2: log-likelihood with log y ──────────────────────────────────────
  p2 <- p1 +
    labs(
      title = sprintf("Scenario: %s (log y)", scenario_name),
      y     = expression( Final~log-likelihood~"(log scale)" )
    ) +
    scale_y_continuous(labels = scales::comma_format())

  # ── Plot 3: log-likelihood + KS combined ───────────────────────────────────
  plot_df_long <- melt(
    plot_df[, .(k, logLik, ks_statistic)],
    id.vars = "k", variable.name = "metric", value.name = "value"
  )
  plot_df_long$metric <- factor(plot_df_long$metric,
    levels = c("logLik", "ks_statistic"),
    labels = c("Log-likelihood", "KS statistic")
  )

    p3 <- ggplot(plot_df_long, aes(x = k, y = value)) +
      geom_vline(xintercept = best_k, color = selected_k_color, linetype = "dashed",
                 linewidth = 0.8)
    if (!is.na(true_k_noise)) {
      p3 <- p3 + geom_vline(xintercept = true_k_noise, color = "red", linetype = "twodash", linewidth = 0.8)
    }
    p3 <- p3 +
      geom_line(color = "darkorange2", linewidth = 1) +
      facet_wrap(~metric, scales = "free_y", ncol = 1) +
      scale_x_log10() +
      labs(
        title    = sprintf("Scenario: %s — estimate_k diagnostics", scenario_name),
        subtitle = if (is.na(true_k_noise)) "Black dashed = selected k; top = log-likelihood, bottom = KS"
                   else sprintf("Green dashed = correct selected k; Red dashed = true k=%d; top = log-likelihood, bottom = KS", true_k_noise),
        x        = expression(k~("noise height, log scale")),
        y        = "Value"
      ) +
      theme_minimal()

  # ── Save plots ─────────────────────────────────────────────────────────────
  if (save_plots) {
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

    tag <- paste0(scenario_name, collapse = "_")
    tag <- gsub("[^a-zA-Z0-9_]+", "_", tag)

    ggsave(
      filename = file.path(plots_dir, sprintf("%s_loglik_trend.png", tag)),
      plot     = p1, width = 9, height = 5.5, dpi = 150
    )
    cat(sprintf("  Saved plots to %s/\n", plots_dir))
  }

  list(
    plot_df  = plot_df,
    p_loglik = p1,
    p_loglik_logy = p2,
    p_combined    = p3,
    best_k        = best_k,
    k_grid        = k_selection$k_grid,
    ks_scores     = k_selection$ks_scores,
    correlation   = cor_val
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# SCENARIO DEFINITIONS
# ══════════════════════════════════════════════════════════════════════════════

set.seed(42)

# ─── Viroli-style scenarios (3×5, 3 groups, permutation outliers) ─────────────
v_r <- 3; v_p <- 5; v_n <- 300

M1_v_base <- matrix(c( 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = v_r, ncol = v_p, byrow = FALSE)
M2_v_base <- matrix(0, v_r, v_p)
M3_v_base <- matrix(c(-0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = v_r, ncol = v_p, byrow = FALSE)

U1v <- rcorrmatrix(v_r); V1v <- rcorrmatrix(v_p)
U2v <- rcorrmatrix(v_r); V2v <- rcorrmatrix(v_p)
U3v <- rcorrmatrix(v_r); V3v <- rcorrmatrix(v_p)

# Scenario V1: small signal (first-column means 0.05, 0 −0.05)
sim_v1 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = matrix(c(0.05, 0.05, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  M2 = M2_v_base, M3 = matrix(c(-0.05, 0.05, 0, rep(0, 12)),
                               nrow = v_r, ncol = v_p, byrow = FALSE),
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V2: reduced covariance (row/col noise scaled to 0.3)
U1v2 <- 0.3 * rcorrmatrix(v_r); V1v2 <- 0.3 * rcorrmatrix(v_p)
U2v2 <- 0.3 * rcorrmatrix(v_r); V2v2 <- 0.3 * rcorrmatrix(v_p)
U3v2 <- 0.3 * rcorrmatrix(v_r); V3v2 <- 0.3 * rcorrmatrix(v_p)
sim_v2 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v2, V1 = V1v2, U2 = U2v2, V2 = V2v2, U3 = U3v2, V3 = V3v2,
  n_outliers = 15
)

# Scenario V3: fewer permutation outliers (5 instead of 15)
set.seed(42 + 1)
sim_v3 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 5
)

# Scenario V4: higher overlap (smaller means + larger covariance)
U1v4 <- 1.2 * rcorrmatrix(v_r); V1v4 <- 1.2 * rcorrmatrix(v_p)
U2v4 <- 1.2 * rcorrmatrix(v_r); V2v4 <- 1.2 * rcorrmatrix(v_p)
U3v4 <- 1.2 * rcorrmatrix(v_r); V3v4 <- 1.2 * rcorrmatrix(v_p)
sim_v4 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3 * v_n) - round(0.4 * v_n),
  M1 = matrix(c(0.15, 0.15, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  M2 = M2_v_base,
  M3 = matrix(c(-0.15, 0.15, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  U1 = U1v4, V1 = V1v4, U2 = U2v4, V2 = V2v4, U3 = U3v4, V3 = V3v4,
  n_outliers = 15
)

# ─── Two-group contamination scenarios (4×4, M1=1, M2=−1) ─────────────────────
c_r <- 4; c_p <- 4
M1_c <- matrix(1,    c_r, c_p)
M2_c <- matrix(-1,   c_r, c_p)

# Scenario C1: low uniform contamination (5 out of 65)
sim_c1 <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = M1_c, M2 = M2_c, n_outliers = 5, n_outliers_perm = 0
)

# Scenario C2: medium positive contamination matching vignette (15 ut of 75)
sim_c2 <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = M1_c, M2 = M2_c, n_outliers = 15, n_outliers_perm = 0
)

# Scenario C3: high contamination with both contamination types (15 uniform + 15 permuted)
sim_c3 <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = M1_c, M2 = M2_c, n_outliers = 15, n_outliers_perm = 15
)

# Scenario C4: structured contamination (uniform + row spikes)
sim_c4 <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = M1_c, M2 = M2_c,
  n_outliers = 10,
  n_outliers_perm = 0,
  n_row_spike = 20
)

# ══════════════════════════════════════════════════════════════════════════════
# RUN ALL SCENARIOS
# ══════════════════════════════════════════════════════════════════════════════

results <- list()

results$V1_small_signal <- future(run_loglik_scenario(
  x_list      = sim_v1$x_list,
  g           = 3,
  scenario_name = "Viroli_small_signal",
  subtitle    = "Viroli 3×5 | 3 groups | first-col means 0.05,0,−0.05 | 15 outliers",
  r           = v_r, p = v_p
  ))

results$V2_low_covariance <- future(run_loglik_scenario(
  x_list      = sim_v2$x_list,
  g           = 3,
  scenario_name = "Viroli_low_cov",
  subtitle    = "Viroli 3×5 | 3 groups | row/col SD = 0.3 instead of 0.5 | 15 outliers",
  r           = v_r, p = v_p
  ))

results$V3_few_outliers <- future(run_loglik_scenario(
  x_list      = sim_v3$x_list,
  g           = 3,
  scenario_name = "Viroli_few_outliers",
  subtitle    = "Viroli 3×5 | 3 groups | 5 permuted-outliers instead of 15",
  r           = v_r, p = v_p
  ))

results$V4_high_overlap <- future(run_loglik_scenario(
  x_list      = sim_v4$x_list,
  g           = 3,
  scenario_name = "Viroli_high_overlap",
  subtitle    = "Viroli 3×5 | weaker means (±0.15) + larger covariance (x1.2) | 15 outliers",
  r           = v_r, p = v_p
  ))

results$C1_low_uniform_contam <- future(run_loglik_scenario(
  x_list      = sim_c1$x_list,
  g           = 2,
  scenario_name = "Contam_low_uniform",
  subtitle    = "Two-group 4×4 | M1=1, M2=-1 | 5 uniform outliers (Unif(-3,3)) | n=65",
  r           = c_r, p = c_p
  ))

results$C2_medium_uniform_contam <- future(run_loglik_scenario(
  x_list      = sim_c2$x_list,
  g           = 2,
  scenario_name = "Contam_medium_uniform",
  subtitle    = "Two-group 4×4 | M1=1, M2=-1 | 15 uniform outliers (Unif(-3,3)) | n=75",
  r           = c_r, p = c_p
  ))

results$C3_high_mixed_contam <- future(run_loglik_scenario(
  x_list      = sim_c3$x_list,
  g           = 2,
  scenario_name = "Contam_high_mixed",
  subtitle    = "Two-group 4×4 | M1=1, M2=-1 | 15 uniform + 15 permuted outliers | n=90",
  r           = c_r, p = c_p
  ))

results$C4_structured_spike_contam <- future(run_loglik_scenario(
  x_list      = sim_c4$x_list,
  g           = 2,
  scenario_name = "Contam_structured_spikes",
  subtitle    = "Two-group 4×4 | M1=1, M2=-1 | 10 uniform + 20 row-spike outliers | n=70",
  r           = c_r, p = c_p
  ))

# ══════════════════════════════════════════════════════════════════════════════
# ADDITIONAL SCENARIOS
# ══════════════════════════════════════════════════════════════════════════════

# ─── Tomarchio-style scenarios (r=2, p=4, column-replacement outliers) ───────
t_r <- 2; t_p <- 4; t_n <- 200
M1_t <- matrix(c(-2.60, -1.10, -0.50, -0.20, 1.30, 0.60, 0.30, 0.10),
               nrow = t_r, ncol = t_p, byrow = FALSE)
M2_t <- matrix(c(1.50, 1.70, 1.90, 2.20, -3.70, -2.70, -2.00, -1.50),
               nrow = t_r, ncol = t_p, byrow = FALSE)
U1_t <- matrix(c(2, 0, 0, 1), nrow = t_r, ncol = t_r)
U2_t <- matrix(c(1.70, 0.5, 0.5, 1.30), nrow = t_r, ncol = t_r)
V_t  <- matrix(c(1.00, 0.50, 0.25, 0.13,
                 0.50, 1.00, 0.50, 0.25,
                 0.25, 0.50, 1.00, 0.50,
                 0.13, 0.25, 0.50, 1.00), nrow = t_p, ncol = t_p)

# Scenario T1: Tomarchio base (equal weight, 10 column outliers)
set.seed(42 + 100)
sim_t1 <- tomarchio_simulation(n = t_n, n_outliers = 10)

# Scenario T2: Tomarchio with additional 15 permuted-entry outliers
set.seed(42 + 101)
sim_t2 <- tomarchio_simulation(n = t_n, n_outliers = 10)
perm_idx <- sample(seq_along(sim_t2$x_list), 15)
for (idx in perm_idx) {
  sim_t2$x_list[[idx]] <- matrix(sample(sim_t2$x_list[[idx]]), t_r, t_p)
}

# Scenario T3: Tomarchio heavier column contamination (20 column outliers)
set.seed(42 + 102)
sim_t3 <- tomarchio_simulation(n = t_n, n_outliers = 20)

# ─── Viroli-style: highly unequal 3-group proportions ────────────────────────
set.seed(42 + 200)
v_n_unequal <- 300
n1_uneq <- round(0.10 * v_n_unequal)   # 30
n2_uneq <- round(0.10 * v_n_unequal)   # 30
n3_uneq <- v_n_unequal - n1_uneq - n2_uneq   # 240
sim_v5 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n_unequal,
  n1 = n1_uneq, n2 = n2_uneq, n3 = n3_uneq,
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V6: Viroli with row-spike outliers (instead of permutation)
set.seed(42 + 201)
sim_v6 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3 * v_n) - round(0.4 * v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 0
)
spike_idx <- sample(seq_along(sim_v6$x_list), 15)
for (idx in spike_idx) {
  row_id <- sample.int(v_r, 1)
  sim_v6$x_list[[idx]][row_id, ] <- runif(v_p, -8, 8)
}

# Scenario V7: Viroli with column outliers (one column replaced per outlier)
set.seed(42 + 202)
sim_v7 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3 * v_n) - round(0.4 * v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 0
)
col_out_idx <- sample(seq_along(sim_v7$x_list), 15)
for (idx in col_out_idx) {
  col_id <- sample.int(v_p, 1)
  sim_v7$x_list[[idx]][, col_id] <- runif(v_r, -10, 10)
}

# ─── Two-group: highly imbalanced proportions ───────────────────────────────
# Scenario C5: very imbalanced 10:60 (n=70)
set.seed(42 + 300)
sim_c5 <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 10, n2 = 60,
  M1 = M1_c, M2 = M2_c, n_outliers = 5, n_outliers_perm = 0
)

# ─── Two-group: strong vs. weak separation ───────────────────────────────────
# Scenario C6_strong: well-separated means (M1=2, M2=-2)
set.seed(42 + 301)
sim_c6_strong <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = matrix(2, c_r, c_p), M2 = matrix(-2, c_r, c_p),
  n_outliers = 10, n_outliers_perm = 0
)

# Scenario C7_weak: weakly separated means (M1=0.25, M2=-0.25)
set.seed(42 + 302)
sim_c7_weak <- two_group_uniform_simulation(
  r = c_r, p = c_p, n1 = 40, n2 = 20,
  M1 = matrix(0.25, c_r, c_p), M2 = matrix(-0.25, c_r, c_p),
  n_outliers = 10, n_outliers_perm = 0
)

results$T1_tomarchio_base <- future(run_loglik_scenario(
  x_list      = sim_t1$x_list,
  g           = 2,
  scenario_name = "Tomarchio_base",
  subtitle    = "Tomarchio 2x4 | 2 groups | column-replacement outliers (Unif(-15,15)) | n=200",
  r           = t_r, p = t_p
  ))

results$T2_tomarchio_plus_permuted <- future(run_loglik_scenario(
  x_list      = sim_t2$x_list,
  g           = 2,
  scenario_name = "Tomarchio_plus_permuted",
  subtitle    = "Tomarchio 2x4 | 10 column + 15 permuted outliers | n=200",
  r           = t_r, p = t_p
  ))

results$T3_tomarchio_heavy_col <- future(run_loglik_scenario(
  x_list      = sim_t3$x_list,
  g           = 2,
  scenario_name = "Tomarchio_heavy_column",
  subtitle    = "Tomarchio 2x4 | 20 column-replacement outliers | n=200",
  r           = t_r, p = t_p
  ))

results$V5_extreme_imbalance <- future(run_loglik_scenario(
  x_list      = sim_v5$x_list,
  g           = 3,
  scenario_name = "Viroli_extreme_imbalance",
  subtitle    = "Viroli 3x5 | 3 groups with π=(0.10,0.10,0.80) | 15 permuted outliers",
  r           = v_r, p = v_p
  ))

results$V6_row_spike_outliers <- future(run_loglik_scenario(
  x_list      = sim_v6$x_list,
  g           = 3,
  scenario_name = "Viroli_row_spike",
  subtitle    = "Viroli 3x5 | 3 groups | 15 row-spike outliers (one noisy row) | 0 permuted",
  r           = v_r, p = v_p
  ))

results$V7_column_outliers <- future(run_loglik_scenario(
  x_list      = sim_v7$x_list,
  g           = 3,
  scenario_name = "Viroli_column_outliers",
  subtitle    = "Viroli 3x5 | 3 groups | 15 column-outliers (one noisy column) | 0 permuted",
  r           = v_r, p = v_p
  ))

results$C5_highly_imbalanced <- future(run_loglik_scenario(
  x_list      = sim_c5$x_list,
  g           = 2,
  scenario_name = "Contam_highly_imbalanced",
  subtitle    = "Two-group 4x4 | M1=1,M2=-1 | very imbalanced (10 vs 60) | 5 uniform outliers",
  r           = c_r, p = c_p
  ))

results$C6_strong_separation <- future(run_loglik_scenario(
  x_list      = sim_c6_strong$x_list,
  g           = 2,
  scenario_name = "Contam_strong_separation",
  subtitle    = "Two-group 4x4 | M1=2,M2=-2 | 10 uniform outliers | n=70",
  r           = c_r, p = c_p
  ))

results$C7_weak_separation <- future(run_loglik_scenario(
  x_list      = sim_c7_weak$x_list,
  g           = 2,
  scenario_name = "Contam_weak_separation",
  subtitle    = "Two-group 4x4 | M1=0.25,M2=-0.25 | 10 uniform outliers | n=70",
  r           = c_r, p = c_p
  ))

# ══════════════════════════════════════════════════════════════════════════════
# SCALED SIZE-SENSITIVITY SCENARIOS
# ══════════════════════════════════════════════════════════════════════════════

# 6 grid points: base size → 5 larger sizes
size_grid <- list(
  c(3, 5), c(4, 6), c(5, 7), c(6, 8), c(7, 9), c(8, 10)
)

cat("\n===== Running scaled size-sensitivity grids =====\n")

for (s in seq_along(size_grid)) {
  rg <- size_grid[[s]][1]
  pg <- size_grid[[s]][2]
  n_total <- round(5 * rg * pg / 3)

  cat(sprintf("\n  >>> Size iteration %d: r=%d, p=%d (n=%d)\n", s - 1, rg, pg, n_total))

  set.seed(42 + 700 + s)
  sim_viroli_small <- scaled_viroli_simulation(
    r = rg, p = pg, n = n_total,
    n1 = round(0.3 * n_total), n2 = round(0.4 * n_total),
    n3 = n_total - round(0.3 * n_total) - round(0.4 * n_total),
    n_outliers = 15, signal_strength = 0.5, cov_scale = 1,
    outlier_type = "perm"
  )
  results[[sprintf("V1_small_signal_sz%d", s - 1)]] <- future(run_loglik_scenario(
    x_list = sim_viroli_small$x_list, g = 3,
    scenario_name = sprintf("V1_small_signal_sz%d", s - 1),
    subtitle = sprintf("Scaling %d: Viroli %dx%d | 3 groups | 15 permuted outliers | n=%d", s - 1, rg, pg, n_total),
    r = rg, p = pg, true_k_noise = 15
  ))

  set.seed(42 + 800 + s)
  sim_two_group_scaled <- scaled_two_group_simulation(
    r = rg, p = pg, n1 = round(0.6 * n_total), n2 = round(0.4 * n_total),
    n_outliers = 10, mean_diff = 2
  )
  results[[sprintf("C_base_2group_sz%d", s - 1)]] <- future(run_loglik_scenario(
    x_list = sim_two_group_scaled$x_list, g = 2,
    scenario_name = sprintf("C_base_2group_sz%d", s - 1),
    subtitle = sprintf("Scaling %d: Two-group %dx%d | M1=1,M2=-1 | 10 uniform outliers | n=%d", s - 1, rg, pg, n_total),
    r = rg, p = pg, true_k_noise = 10
  ))

  if (s - 1 >= 1) {
    set.seed(42 + 900 + s)
    sim_tomarchio_scaled <- scaled_tomarchio_simulation(
       r = rg, p = pg, n = n_total, n_outliers = 10
    )
    results[[sprintf("T_base_tomarchio_sz%d", s - 1)]] <- future(run_loglik_scenario(
      x_list = sim_tomarchio_scaled$x_list, g = 2,
      scenario_name = sprintf("T_base_tomarchio_sz%d", s - 1),
      subtitle = sprintf("Scaling %d: Tomarchio %dx%d | 2 groups | 10 col outliers | n=%d", s - 1, rg, pg, n_total),
      r = rg, p = pg, true_k_noise = 10
    ))
  }
}

cat("  Resolving all futures ...\n")
for (nm in names(results)) {
  results[[nm]] <- value(results[[nm]])
}

summary_data <- rbindlist(lapply(seq_along(results), function(i) {
  nm <- names(results)[[i]]
  r  <- results[[nm]]
  data.table(
    idx           = i,
    Scenario      = nm,
    g             = if (grepl("Tomarchio", nm)) 2 else if (grepl("^V", nm)) 3 else 2,
    Selected_k    = r$best_k,
    Max_logLik_k  = ifelse(any(r$plot_df$logLik[complete.cases(r$plot_df$logLik)] > -Inf),
                          r$plot_df$k[complete.cases(r$plot_df$logLik)][
                            which.max(r$plot_df$logLik[complete.cases(r$plot_df$logLik)])
                          ], NA_real_),
    Spearman_k_vs_logLik = r$correlation,
    N_candidates  = nrow(r$plot_df),
    N_valid_logLik= sum(complete.cases(r$plot_df$logLik))
  )
}), fill = TRUE)

setcolorder(summary_data, c("idx", "Scenario", "g", "Selected_k", "Max_logLik_k",
                             "Spearman_k_vs_logLik", "N_candidates", "N_valid_logLik"))

cat("\n===== Cross-scenario summary =====\n")
print(summary_data)

cat(sprintf("\nSaved %d scenario trend plots under %s/\n",
            length(results), normalizePath("loglik/plots")))
