library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)
library(future)
library(future.apply)

plan(multisession, workers = availableCores() - 1)

# ──────────────────────────────────────────────────────────────────────────────
# Shared helper — evaluate one k-grid candidate and return log-likelihood
# ──────────────────────────────────────────────────────────────────────────────
evaluate_k_candidate_with_loglik <- function(idx, x_list, r, p, g,
                                              max_iter, tol,
                                              k_grid, nstart, init,
                                              noise_pi_init, verbose,
                                              outlier_idx) {
  current_k <- k_grid[idx]

  fit_noise <-  Ampharos:::mv_noise_fit_impl(
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

  detected_noise <- which(fit_noise$cluster == 0)
  correct_noise <- length(detected_noise) == length(outlier_idx) &&
                   length(intersect(detected_noise, outlier_idx)) == length(outlier_idx)

  if (length(x_clean) <= g) {
    return(list(
      k            = current_k,
      logLik       = NA_real_,
      ks_statistic = Inf,
      ks_p.value   = NA_real_,
      n_used       = length(x_clean),
      correct_noise = correct_noise
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
      n_used       = length(x_clean),
      correct_noise = correct_noise
    ))
  }

  ks_result <- suppressWarnings(
    tryCatch(
      Ampharos:::mv_noise_ks_score(fit_clean, x_clean),
      error = function(e) list(statistic = Inf, p.value = NA_real_,
                                n_used = length(x_clean))
    )
  )

  list(
    k            = current_k,
    logLik       = tail(fit_clean$logLik, 1),
    ks_statistic = ks_result$statistic,
    ks_p.value   = ks_result$p.value,
    n_used       = ks_result$n_used,
    correct_noise = correct_noise
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

  outlier_idx <- integer(0)
  if (n_outliers > 0) {
    contam <- lapply(seq_len(n_outliers),
                     function(i) matrix(runif(r * p, -3, 3), r, p))
    x_list <- c(x_list, contam)
    outlier_idx <- c(outlier_idx, (length(x_list) - n_outliers + 1):length(x_list))
  }
  if (n_outliers_perm > 0) {
    perm_idx <- sample(seq_along(x_list), n_outliers_perm)
    for (idx in perm_idx) {
      x_list[[idx]] <- matrix(sample(x_list[[idx]]), r, p)
    }
    outlier_idx <- union(outlier_idx, perm_idx)
  }
  if (n_row_spike > 0) {
    spike_idx <- sample(seq_along(x_list), n_row_spike)
    for (idx in spike_idx) {
      row_id <- sample.int(r, 1)
      x_list[[idx]][row_id, ] <- runif(p, -8, 8)
    }
    outlier_idx <- union(outlier_idx, spike_idx)
  }
  if (n_col_out > 0) {
    col_idx <- sample(seq_along(x_list), n_col_out)
    for (idx in col_idx) {
      col_id <- sample.int(p, 1)
      x_list[[idx]][, col_id] <- runif(r, col_out_range[1], col_out_range[2])
    }
    outlier_idx <- union(outlier_idx, col_idx)
  }
  list(x_list = x_list, outlier_idx = outlier_idx)
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
  list(x_list = x_list, outlier_idx = out_idx)
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
  list(x_list = result$x_list, outlier_idx = result$outlier_idx)
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
  list(x_list = x_list, outlier_idx = out_idx)
}
# ──────────────────────────────────────────────────────────────────────────────
run_loglik_scenario <- function(x_list, g, scenario_name,
                                 subtitle, r, p,
                                 true_k_noise = NA_integer_,
                                 outlier_idx = NULL,
                                 nstart = 10, max_iter = 100, tol = 1e-6,
                                 use_kmeans = TRUE, init = "kmeans",
                                 noise_pi_init = 0.05,
                                 save_plots = TRUE, plots_dir = "loglik/plots") {
   for (pkg in c("Ampharos", "ggplot2", "clusterGeneration", "data.table")) {
     if (!requireNamespace(pkg, quietly = TRUE)) stop("Package ", pkg, " required")
   }
       if (!is.null(outlier_idx) && !is.na(true_k_noise) && length(outlier_idx) != true_k_noise) {
        warning("outlier_idx length (", length(outlier_idx), ") != true_k_noise (", true_k_noise, ")")
      }
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
    verbose       = FALSE,
    outlier_idx   = outlier_idx
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

    n_correct <- sum(plot_df$correct_noise, na.rm = TRUE)
    cat(sprintf("  k-grid values with perfect noise recovery: %d / %d\n",
                n_correct, nrow(plot_df)))

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
    plot_df$correct_noise_group <- cumsum(c(TRUE, diff(plot_df$correct_noise) != 0))

    p1 <- ggplot(plot_df, aes(x = k, y = logLik)) +
      geom_vline(xintercept = best_k, color = selected_k_color, linetype = "dashed",
                 linewidth = 0.8) +
      geom_vline(xintercept = k_selection$k_grid[best_ks_idx],
                 color = "steelblue", linetype = "dotdash", linewidth = 0.8) +
      geom_line(aes(color = correct_noise, group = correct_noise_group),
                linewidth = 1) +
      geom_point(data = subset(plot_df, type != "Other"),
                 aes(shape = type), color = "darkorange2", size = 3) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkorange2"),
                         labels = c("Correct noise recovery", "Incorrect recovery"),
                         name = "Noise recovery") +
      scale_x_log10() +
      scale_y_continuous(labels = scales::comma_format()) +
      labs(
        title    = sprintf("Scenario: %s", scenario_name),
        subtitle = subtitle,
        x        = expression(k~("noise height, log scale")),
        y        = expression( Final~log-likelihood ),
        shape    = "",
        colour   = "Noise recovery"
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
    plot_df[, .(k, logLik, ks_statistic, correct_noise, correct_noise_group)],
    id.vars = c("k", "correct_noise", "correct_noise_group"),
    variable.name = "metric", value.name = "value"
  )
  plot_df_long$metric <- factor(plot_df_long$metric,
    levels = c("logLik", "ks_statistic"),
    labels = c("Log-likelihood", "KS statistic")
  )

    p3 <- ggplot(plot_df_long, aes(x = k, y = value)) +
      geom_vline(xintercept = best_k, color = selected_k_color, linetype = "dashed",
                 linewidth = 0.8) +
      geom_line(aes(color = correct_noise, group = interaction(metric, correct_noise_group)),
                linewidth = 1) +
      facet_wrap(~metric, scales = "free_y", ncol = 1) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkorange2"),
                         labels = c("Correct noise recovery", "Incorrect recovery"),
                         name = "Noise recovery") +
      scale_x_log10() +
      labs(
        title    = sprintf("Scenario: %s — estimate_k diagnostics", scenario_name),
        subtitle = if (is.na(true_k_noise)) "Black dashed = selected k; red = correct noise recovery"
                   else sprintf("Green dashed = correct selected k; red = correct noise recovery"),
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

source("loglik/scenarios.R")
