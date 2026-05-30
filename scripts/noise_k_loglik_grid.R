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
fine_half_width_log10 <- 1
fine_grid_n <- 31
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
    coarse_plot_file <- file.path(output_dir, paste0("ks_curve_coarse_", run_label, ".png"))
    coarse_curve_file <- file.path(output_dir, paste0("ks_curve_coarse_", run_label, ".csv"))
    refined_plot_file <- file.path(output_dir, paste0("ks_curve_refined_", run_label, ".png"))
    refined_curve_file <- file.path(output_dir, paste0("ks_curve_refined_", run_label, ".csv"))

    status <- "ok"
    error_message <- NA_character_

    if (length(finite_ks) == 0) {
      status <- "no_finite_ks"
      error_message <- if (!is.na(run_error)) run_error else "No finite KS statistics were produced for this run."
    } else {
      selected_idx <- which.min(ifelse(is.finite(results$ks_statistic), results$ks_statistic, Inf))
      selected_k <- results$noise_k[selected_idx]
      selected_ks <- results$ks_statistic[selected_idx]

      fine_lower <- max(.Machine$double.xmin, selected_k / (10^fine_half_width_log10))
      fine_upper <- selected_k * (10^fine_half_width_log10)
      fine_noise_k_grid <- exp(seq(log(fine_lower), log(fine_upper), length.out = fine_grid_n))

      fine_results <- data.frame(
        noise_k = fine_noise_k_grid,
        ks_statistic = NA_real_,
        converged = NA,
        stringsAsFactors = FALSE
      )

      fine_run_error <- NA_character_
      for (i in seq_along(fine_noise_k_grid)) {
        candidate_k <- fine_noise_k_grid[i]
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
          fine_run_error <- conditionMessage(fit)
          fine_results$ks_statistic[i] <- NA_real_
          fine_results$converged[i] <- NA
        } else {
          fine_results$ks_statistic[i] <- matrix_noise_ks_score(fit, x_list)$statistic
          fine_results$converged[i] <- isTRUE(fit$converged)
        }
      }

      fine_finite_ks <- fine_results$ks_statistic[is.finite(fine_results$ks_statistic)]
      if (length(fine_finite_ks) == 0) {
        status <- "no_finite_ks_refined"
        error_message <- if (!is.na(fine_run_error)) fine_run_error else "No finite KS statistics were produced in the refined grid."
      } else {
        fine_selected_idx <- which.min(ifelse(is.finite(fine_results$ks_statistic), fine_results$ks_statistic, Inf))
        fine_selected_k <- fine_results$noise_k[fine_selected_idx]
        fine_selected_ks <- fine_results$ks_statistic[fine_selected_idx]

        fine_results$selected <- FALSE
        fine_results$selected[fine_selected_idx] <- TRUE
        fine_results$nrow_matrix <- r
        fine_results$ncol_matrix <- p
        fine_results$n_contam <- n_contam
        fine_results$run_label <- run_label
        fine_results$grid_type <- "fine"

        selected_k <- fine_selected_k
        selected_ks <- fine_selected_ks

        results$selected <- FALSE
        results$selected[selected_idx] <- TRUE
        results$nrow_matrix <- r
        results$ncol_matrix <- p
        results$n_contam <- n_contam
        results$run_label <- run_label
        results$grid_type <- "coarse"

        grDevices::png(coarse_plot_file, width = 900, height = 650, res = 120)
        coarse_plot_title <- paste0("Coarse HC noise_k vs KS statistic (", r, " x ", p, ", contam = ", n_contam, ")")
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
          main = coarse_plot_title
        )
        points(selected_k, selected_ks, col = "red", pch = 19, cex = 1.4)
        legend("topright", legend = "Selected k", col = "red", pch = 19, bty = "n")
        grid()
        grDevices::dev.off()

        write.csv(results, coarse_curve_file, row.names = FALSE)

        grDevices::png(refined_plot_file, width = 900, height = 650, res = 120)
        plot_title <- paste0("Refined HC noise_k vs KS statistic (", r, " x ", p,", contam = ", n_contam, ")")
        plot(
          fine_results$noise_k,
          fine_results$ks_statistic,
          log = "x",
          type = "b",
          pch = 19,
          col = "steelblue4",
          ylim = range(fine_finite_ks),
          xlab = expression(noise[k]),
          ylab = "KS statistic",
          main = plot_title
        )
        points(fine_selected_k, fine_selected_ks, col = "red", pch = 19, cex = 1.4)
        legend("topright", legend = "Selected k", col = "red", pch = 19, bty = "n")
        grid()
        grDevices::dev.off()

        write.csv(fine_results, refined_curve_file, row.names = FALSE)
        all_curve_rows[[length(all_curve_rows) + 1L]] <- fine_results
      }
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
      coarse_plot_file = if (status == "ok") coarse_plot_file else NA_character_,
      coarse_curve_file = if (status == "ok") coarse_curve_file else NA_character_,
      refined_plot_file = if (status == "ok") refined_plot_file else NA_character_,
      refined_curve_file = if (status == "ok") refined_curve_file else NA_character_,
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