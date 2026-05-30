set.seed(123)

n_group1 <- 20
n_group2 <- 20
n_contam_grid <- 2:15
size_grid <- list(
  c(3, 5),
  c(4, 6),
  c(5, 8)
)
row_sd <- 0.5
col_sd <- 0.5
g <- 2
n_grid <- 50
max_iter <- 100
tol <- 1e-06
nstart <- 10
noise_pi_init <- 0.05
verbose <- FALSE

repo_root <- normalizePath(getwd())
script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grep("^--file=", script_args)])
if (length(script_file) > 0) {
  repo_root <- normalizePath(file.path(dirname(script_file[1]), ".."))
} else if (file.exists(file.path(repo_root, "R"))) {
  repo_root <- normalizePath(repo_root)
} else if (file.exists(file.path(dirname(repo_root), "R"))) {
  repo_root <- normalizePath(file.path(repo_root, ".."))
} else {
  stop("Could not determine the repository root containing the R/ directory.")
}
r_dir <- file.path(repo_root, "R")
r_files <- sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))
if (length(r_files) == 0) {
  stop("Could not find package R files under the R/ directory.")
}
invisible(lapply(r_files, source, local = globalenv()))

output_dir <- file.path(repo_root, "results", paste0("noise_k_batch_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

noise_k_grid <- exp(seq(log(1e-01), log(.Machine$double.xmin), length.out = n_grid))

all_curve_rows <- list()
summary_rows <- list()
run_counter <- 0L

for (size_idx in seq_along(size_grid)) {
  size_pair <- size_grid[[size_idx]]
  nrow_matrix <- size_pair[1]
  ncol_matrix <- size_pair[2]
  r <- nrow_matrix
  p <- ncol_matrix

  M1 <- matrix(1, r, p)
  M2 <- matrix(-1, r, p)

  for (n_contam in n_contam_grid) {
    run_counter <- run_counter + 1L
    run_label <- paste0("r", r, "_p", p, "_contam", n_contam)

    x_list <- c(
      lapply(seq_len(n_group1), function(i) {
        mean_matrix <- M1
        row_cov <- diag(row_sd, nrow(mean_matrix))
        col_cov <- diag(col_sd, ncol(mean_matrix))
        mean_matrix + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov
      }),
      lapply(seq_len(n_group2), function(i) {
        mean_matrix <- M2
        row_cov <- diag(row_sd, nrow(mean_matrix))
        col_cov <- diag(col_sd, ncol(mean_matrix))
        mean_matrix + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov
      }),
      lapply(seq_len(n_contam), function(i) matrix(runif(r * p, -15, 15), r, p))
    )

    results <- data.frame(
      noise_k = noise_k_grid,
      ks_statistic = NA_real_,
      converged = NA,
      stringsAsFactors = FALSE
    )

    run_error <- NA_character_
    for (i in seq_along(noise_k_grid)) {
      candidate_k <- noise_k_grid[i]
      fit <- tryCatch(
        matrix_variate_noise_fit(
          x_list = x_list,
          g = g,
          noise_type = "hc",
          noise_k = candidate_k,
          max_iter = max_iter,
          tol = tol,
          nstart = nstart,
          noise_pi_init = noise_pi_init,
          verbose = verbose
        ),
        error = function(e) e
      )

      if (inherits(fit, "error")) {
        run_error <- conditionMessage(fit)
        results$ks_statistic[i] <- NA_real_
        results$converged[i] <- NA
      } else {
        results$ks_statistic[i] <- matrix_noise_ks_score(fit, x_list)$statistic
        results$converged[i] <- isTRUE(fit$converged)
      }
    }

    finite_ks <- results$ks_statistic[is.finite(results$ks_statistic)]
    selected_idx <- NA_integer_
    selected_k <- NA_real_
    selected_ks <- NA_real_
    plot_file <- file.path(output_dir, paste0("ks_curve_", run_label, ".png"))
    curve_file <- file.path(output_dir, paste0("ks_curve_", run_label, ".csv"))

    status <- "ok"
    error_message <- NA_character_

    if (length(finite_ks) == 0) {
      status <- "no_finite_ks"
      error_message <- if (!is.na(run_error)) run_error else "No finite KS statistics were produced for this run."
    } else {
      selected_idx <- which.min(ifelse(is.finite(results$ks_statistic), results$ks_statistic, Inf))
      selected_k <- results$noise_k[selected_idx]
      selected_ks <- results$ks_statistic[selected_idx]

      results$selected <- FALSE
      results$selected[selected_idx] <- TRUE
      results$nrow_matrix <- r
      results$ncol_matrix <- p
      results$n_contam <- n_contam
      results$run_label <- run_label

      grDevices::png(plot_file, width = 900, height = 650, res = 120)
      plot_title <- paste0("HC noise_k vs KS statistic (", r, " x ", p, ", contam = ", n_contam, ")")
      plot(
        results$noise_k,
        results$ks_statistic,
        log = "x",
        type = "b",
        pch = 19,
        col = "steelblue4",
        ylim = range(finite_ks),
        xlab = expression(noise[k]),
        ylab = "KS statistic",
        main = plot_title
      )
      points(selected_k, selected_ks, col = "red", pch = 19, cex = 1.4)
      legend("topright", legend = "Selected k", col = "red", pch = 19, bty = "n")
      grid()
      grDevices::dev.off()

      write.csv(results, curve_file, row.names = FALSE)
      all_curve_rows[[length(all_curve_rows) + 1L]] <- results
    }

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      run_id = run_counter,
      run_label = run_label,
      nrow_matrix = r,
      ncol_matrix = p,
      n_contam = n_contam,
      selected_k = selected_k,
      selected_ks = selected_ks,
      status = status,
      error_message = error_message,
      plot_file = if (status == "ok") plot_file else NA_character_,
      curve_file = if (status == "ok") curve_file else NA_character_,
      stringsAsFactors = FALSE
    )
  }
}

summary_results <- do.call(rbind, summary_rows)
curve_results <- if (length(all_curve_rows) > 0) do.call(rbind, all_curve_rows) else data.frame()

write.csv(summary_results, file.path(output_dir, "noise_k_summary.csv"), row.names = FALSE)
write.csv(curve_results, file.path(output_dir, "noise_k_all_curves.csv"), row.names = FALSE)
saveRDS(summary_results, file.path(output_dir, "noise_k_summary.rds"))
saveRDS(curve_results, file.path(output_dir, "noise_k_all_curves.rds"))

summary_results