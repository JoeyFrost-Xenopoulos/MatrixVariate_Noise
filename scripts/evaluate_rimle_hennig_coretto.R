#!/usr/bin/env Rscript
# Evaluate RIMLE estimate_k with HENNIG-CORETTO initialization
# Viroli noise setup across matrix dimensions (3,3) through (10,13)

# ---- Load packages ----
if (!requireNamespace("RIMLEMV", quietly = TRUE)) {
  stop("Package 'RIMLEMV' is required. Install it first.")
}
if (!requireNamespace("clusterGeneration", quietly = TRUE)) {
  stop("Package 'clusterGeneration' is required. Install it first.")
}

library(RIMLEMV)
library(clusterGeneration)

# ---- Parameters ----
set.seed(42)
n <- 300
g <- 3
n_noise <- 15
n_reps <- 1
init_method <- "hennig-coretto"
q_initial <- 3

# ---- Dimensions (r, p) ----
dims <- list(
  c(3, 3), c(3, 4), c(3, 5),
  c(4, 3), c(4, 4), c(4, 5),
  c(5, 5), c(5, 6), c(5, 7),
  c(6, 6), c(6, 7), c(6, 8),
  c(7, 7), c(7, 8), c(7, 9),
  c(8, 8), c(8, 9), c(8, 10),
  c(9, 9), c(9, 10), c(9, 11),
  c(10, 10), c(10, 11), c(10, 12), c(10, 13)
)

# ---- Generate Viroli-style data for any (r, p) ----
generate_viroli_data <- function(r, p, n = 300, n_noise = 15, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  props <- c(0.3, 0.4, 0.3)
  n_g <- round(props * n)
  n_g[3] <- n - sum(n_g[1:2])

  M1 <- matrix(0, r, p)
  M1[1, 1] <- 0.5
  if (r >= 2) M1[2, 1] <- 0.5

  M2 <- matrix(0, r, p)

  M3 <- matrix(0, r, p)
  M3[1, 1] <- -0.5
  if (r >= 2) M3[2, 1] <- 0.5

  U1 <- rcorrmatrix(r)
  V1 <- rcorrmatrix(p)
  U2 <- rcorrmatrix(r)
  V2 <- rcorrmatrix(p)
  U3 <- rcorrmatrix(r)
  V3 <- rcorrmatrix(p)

  simulate_group <- function(n, mean_mat, row_cov, col_cov) {
    lapply(seq_len(n), function(i) {
      mean_mat + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov
    })
  }

  x_list <- c(
    simulate_group(n_g[1], M1, U1, V1),
    simulate_group(n_g[2], M2, U2, V2),
    simulate_group(n_g[3], M3, U3, V3)
  )

  outlier_idx <- sample(length(x_list), n_noise)
  for (idx in outlier_idx) {
    x_list[[idx]] <- matrix(sample(x_list[[idx]]), r, p)
  }

  list(x_list = x_list, outlier_idx = outlier_idx)
}

# ---- Run a single estimate_k fit ----
run_single_fit <- function(x_list, init, q = 3) {
  fit <- tryCatch({
    RIMLEMV::rimle_fit(
      x_list = x_list,
      g = g,
      estimate_k = TRUE,
      init = init,
      q = q,
      verbose = TRUE,
      nstart = 10,
      pi_max = 0.5,
      gamma = 1000,
      max_iter = 100,
      tol = 1e-6,
      use_parallel = TRUE,
      n_cores = 8
    )
  }, error = function(e) {
    message("Fit failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(fit)) return(NULL)

  noise_count <- sum(fit$cluster == 0)
  correct_noise <- noise_count == n_noise

  list(
    fit = fit,
    correct_noise = correct_noise,
    noise_count = noise_count,
    selected_k = fit$k_selection$selected_k,
    k_grid = fit$k_selection$k_grid,
    n_used = fit$k_selection$n_used,
    loglik = tail(fit$logLik, 1),
    converged = fit$converged,
    iterations = fit$iterations
  )
}

# ---- Parse final q from verbose output (Hennig-Coretto retries) ----
parse_q_from_verbose <- function(captured, q_initial = 3) {
  pattern <- "Retrying k-grid search with q = ([0-9]+)"
  matches <- regmatches(captured, regexec(pattern, captured))
  qs <- sapply(matches, function(m) if (length(m) > 1) as.integer(m[2]) else NA_integer_)
  qs <- qs[!is.na(qs)]
  if (length(qs) == 0) return(q_initial)
  tail(qs, 1)
}

# ---- Parse n_used values from verbose output ----
parse_n_used_from_verbose <- function(captured) {
  lines <- captured[grepl("Testing k =", captured)]
  n_used <- integer(0)
  
  for (line in lines) {
    if (grepl("n_used =", line)) {
      val <- as.integer(gsub(".*n_used = ([0-9]+).*", "\\1", line))
      n_used <- c(n_used, val)
    } else if (grepl("insufficient retained observations", line)) {
      n_used <- c(n_used, NA_integer_)
    } else {
      n_used <- c(n_used, NA_integer_)
    }
  }
  
  n_used
}

# ---- Main simulation loop ----
results <- list()
grid_results <- list()
row_idx <- 1
grid_row_idx <- 1

for (dim_idx in seq_along(dims)) {
  r <- dims[[dim_idx]][1]
  p <- dims[[dim_idx]][2]
  cat(sprintf("Dimension %d/%d: r=%d, p=%d\n", dim_idx, length(dims), r, p))

  for (rep in seq_len(n_reps)) {
    seed <- 1000 + dim_idx * 100 + rep
    cat(sprintf("  Rep %d/%d (seed=%d)...\n", rep, n_reps, seed))

    data_gen <- generate_viroli_data(r, p, n = n, n_noise = n_noise, seed = seed)
    x_list <- data_gen$x_list

    captured <- capture.output(
      fit_result <- run_single_fit(x_list, init = init_method, q = q_initial),
      type = c("output", "message")
    )

    if (is.null(fit_result)) {
      results[[row_idx]] <- list(
        r = r, p = p, rep = rep, seed = seed,
        init = init_method,
        n_noise = n_noise,
        correct_noise = NA,
        noise_count = NA,
        selected_k = NA,
        loglik = NA,
        converged = NA,
        iterations = NA,
        q_final = ifelse(init_method == "hennig-coretto", parse_q_from_verbose(captured, q_initial), NA)
      )
      row_idx <- row_idx + 1
      next
    }

    q_final <- ifelse(init_method == "hennig-coretto",
                      parse_q_from_verbose(captured, q_initial),
                      NA_integer_)

    results[[row_idx]] <- list(
      r = r, p = p, rep = rep, seed = seed,
      init = init_method,
      n_noise = n_noise,
      correct_noise = fit_result$correct_noise,
      noise_count = fit_result$noise_count,
      selected_k = fit_result$selected_k,
      loglik = fit_result$loglik,
      converged = fit_result$converged,
      iterations = fit_result$iterations,
      q_final = q_final
    )

    parsed_n_used <- parse_n_used_from_verbose(captured)
    n_used_vec <- if (length(parsed_n_used) == length(fit_result$k_grid)) {
      parsed_n_used
    } else {
      fit_result$n_used
    }

    for (k_idx in seq_along(fit_result$k_grid)) {
      grid_results[[grid_row_idx]] <- list(
        r = r, p = p, rep = rep, seed = seed,
        init = init_method,
        n_noise = n_noise,
        k_idx = k_idx,
        k_value = fit_result$k_grid[k_idx],
        n_used = n_used_vec[k_idx],
        selected = fit_result$k_grid[k_idx] == fit_result$selected_k,
        q_final = if (init_method == "hennig-coretto") q_final else NA_integer_
      )
      grid_row_idx <- grid_row_idx + 1
    }

    row_idx <- row_idx + 1
  }
}

# ---- Save results ----
summary_df <- do.call(rbind, lapply(results, as.data.frame))
write.csv(summary_df,
          file = sprintf("results_summary_%s.csv", init_method),
          row.names = FALSE)

grid_df <- do.call(rbind, lapply(grid_results, as.data.frame))
write.csv(grid_df,
          file = sprintf("results_grid_%s.csv", init_method),
          row.names = FALSE)

cat(sprintf("Done! Results saved to results_summary_%s.csv and results_grid_%s.csv\n",
            init_method, init_method))
