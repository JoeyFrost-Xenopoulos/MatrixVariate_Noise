#' Extensive HC Noise-k Test Scheme
#'
#' This file defines a simulation harness for studying how the automatically
#' selected HC noise height changes with matrix size, sample size, and the
#' amount of contamination.
#'
#' The goal is to make the search for an optimal `k` more empirical: generate
#' many matrix-valued datasets of different shapes, contaminate some of the
#' observations with outliers, fit the matrix-variate noise mixture with
#' `select_noise_k = TRUE`, and then summarize the selected `k` values in a way
#' that makes size-scaling patterns easy to inspect.
#'
#' @param size_grid A list of integer length-2 vectors giving the matrix
#'   dimensions to test. The default covers a broad progression of shapes from
#'   very small to substantially larger matrices.
#' @param n_grid Integer vector of sample sizes to test for each matrix shape.
#' @param contam_grid Integer vector of contamination counts to inject into
#'   each simulated dataset.
#' @param n_sim Integer: number of replicated simulations for each scenario.
#' @param g Integer: number of Gaussian mixture components.
#' @param signal_shift Numeric: separation between the clean Gaussian component
#'   means.
#' @param row_sd Numeric: standard deviation scale for the row-side simulation.
#' @param col_sd Numeric: standard deviation scale for the column-side
#'   simulation.
#' @param contam_low Numeric: lower bound for contaminated matrix entries.
#' @param contam_high Numeric: upper bound for contaminated matrix entries.
#' @param noise_k_grid Numeric vector passed to [matrix_variate_noise_fit()] when
#'   `select_noise_k = TRUE`.
#' @param max_iter Integer: maximum EM iterations.
#' @param tol Numeric: EM convergence tolerance.
#' @param nstart Integer: number of k-means restarts for initialization.
#' @param verbose Logical: print progress messages.
#' @param csv_path Optional path or filename stem for a CSV file to write the
#'   per-simulation results at the end of the run.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A list with three elements:
#'   
#'   * `results`: one row per simulation replicate and scenario.
#'   * `summary`: scenario-level summaries of the selected `k` values.
#'   * `scaling_fit`: a simple regression fit that helps inspect whether the
#'     selected `k` is scaling with dimension, sample size, or contamination.
#'
#' @examples
#' \dontrun{
#' run_optimal_k_test_scheme(
#'   size_grid = default_optimal_k_size_grid(),
#'   n_grid = c(30, 60),
#'   contam_grid = c(0, 5, 10),
#'   n_sim = 5,
#'   csv_path = "output/optimal_k_results.csv",
#'   seed = 1
#' )
#' }
#'
#' @export
run_optimal_k_test_scheme <- function(size_grid = default_optimal_k_size_grid(),
										 n_grid = c(30L, 50L, 100L),
										 contam_grid = c(0L, 5L, 10L, 20L),
										 n_sim = 20L,
										 g = 2L,
										 signal_shift = 6,
										 row_sd = 0.35,
										 col_sd = 0.35,
										 contam_low = -15,
										 contam_high = 15,
										 noise_k_grid = 10^seq(-150, -1, length.out = 25),
										 max_iter = 50L,
										 tol = 1e-06,
										 nstart = 50L,
										 verbose = FALSE,
										 csv_path = NULL,
										 seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	size_grid <- lapply(size_grid, function(item) {
		item <- as.integer(item)
		if (length(item) != 2L || any(!is.finite(item)) || any(item <= 0L)) {
			stop("Each size_grid entry must be a length-2 positive integer vector.")
		}
		item
	})

	n_grid <- as.integer(n_grid)
	contam_grid <- as.integer(contam_grid)
	if (length(n_grid) == 0L || any(n_grid <= 0L)) {
		stop("n_grid must contain at least one positive sample size.")
	}
	if (length(contam_grid) == 0L || any(contam_grid < 0L)) {
		stop("contam_grid must contain at least one non-negative contamination count.")
	}
	if (n_sim <= 0L) {
		stop("n_sim must be positive.")
	}
	if (g < 2L) {
		stop("g must be at least 2.")
	}

	csv_file <- NULL
	if (!is.null(csv_path)) {
		csv_path <- as.character(csv_path)[1]
		if (nzchar(csv_path)) {
			if (dir.exists(csv_path)) {
				csv_file <- file.path(csv_path, "optimal_k_results.csv")
			} else if (grepl("\\.csv$", csv_path, ignore.case = TRUE)) {
				csv_file <- csv_path
			} else {
				csv_file <- paste0(csv_path, ".csv")
			}
			dir.create(dirname(csv_file), recursive = TRUE, showWarnings = FALSE)
		}
	}

	write_optimal_k_results_csv <- function(result_rows, csv_file, verbose = FALSE, completed = FALSE) {
		if (is.null(csv_file)) {
			return(invisible(NULL))
		}
		if (is.null(result_rows) || nrow(result_rows) == 0L) {
			return(invisible(NULL))
		}
		utils::write.csv(result_rows, file = csv_file, row.names = FALSE)
		if (verbose) {
			if (completed) {
				message(sprintf("Wrote optimal-k results to %s", csv_file))
			} else {
				message(sprintf("Wrote partial optimal-k results to %s before stopping due to an error.", csv_file))
			}
		}
		invisible(NULL)
	}

	results <- vector("list", length(size_grid) * length(n_grid) * length(contam_grid) * n_sim)
	result_idx <- 0L

	tryCatch(
		{
			for (size_item in size_grid) {
				r <- size_item[1]
				p <- size_item[2]
				dimension <- r * p

				for (n_obs in n_grid) {
					for (n_contam in contam_grid) {
						n_clean <- n_obs - n_contam
						if (n_clean <= 0L) {
							if (verbose) {
								message(sprintf("Skipping size %dx%d, n = %d, contam = %d because there are no clean observations.", r, p, n_obs, n_contam))
							}
							next
						}

						group_sizes <- rep(floor(n_clean / g), g)
						remainder <- n_clean - sum(group_sizes)
						if (remainder > 0L) {
							group_sizes[seq_len(remainder)] <- group_sizes[seq_len(remainder)] + 1L
						}
						if (any(group_sizes <= 0L)) {
							if (verbose) {
								message(sprintf("Skipping size %dx%d, n = %d, contam = %d because a Gaussian group would be empty.", r, p, n_obs, n_contam))
							}
							next
						}

						for (sim_idx in seq_len(n_sim)) {
							if (verbose && sim_idx %% 5L == 0L) {
								message(sprintf("Testing %dx%d | n = %d | contam = %d | sim %d/%d", r, p, n_obs, n_contam, sim_idx, n_sim))
							}

							dataset <- simulate_optimal_k_dataset(
								r = r,
								p = p,
								group_sizes = group_sizes,
								n_contam = n_contam,
								signal_shift = signal_shift,
								row_sd = row_sd,
								col_sd = col_sd,
								contam_low = contam_low,
								contam_high = contam_high
							)

							fit <- matrix_variate_noise_fit(
								x_list = dataset$x_list,
								g = g,
								noise_type = "hc",
								select_noise_k = TRUE,
								noise_k_grid = noise_k_grid,
								nstart = nstart,
								max_iter = max_iter,
								tol = tol,
								verbose = verbose
							)

							result_idx <- result_idx + 1L
							selected_k <- fit$noise$search$selected_k
							results[[result_idx]] <- data.frame(
								sim = sim_idx,
								r = r,
								p = p,
								dimension = dimension,
								n_obs = n_obs,
								n_clean = n_clean,
								n_contam = n_contam,
								contam_rate = n_contam / n_obs,
								selected_k = selected_k,
								log10_selected_k = log10(selected_k),
								selected_k_per_entry = selected_k^(1 / dimension),
								log10_selected_k_per_entry = log10(selected_k) / dimension,
								ks_statistic = fit$noise$search$ks_statistic,
								ks_p_value = fit$noise$search$ks_p_value,
								noise_rate = mean(fit$cluster == 0L),
								stringsAsFactors = FALSE
							)
						}
					}
				}
			}
			NULL
		},
		error = function(e) {
			partial_results <- if (result_idx > 0L) {
				do.call(rbind, results[seq_len(result_idx)])
			} else {
				NULL
			}
			write_optimal_k_results_csv(partial_results, csv_file, verbose = verbose, completed = FALSE)
			stop(e)
		}
	)

	results <- do.call(rbind, results[seq_len(result_idx)])
	if (is.null(results) || nrow(results) == 0L) {
		stop("No valid simulation scenarios were produced.")
	}

	summary <- summarize_optimal_k_test_scheme(results)
	scaling_fit <- tryCatch(
		stats::lm(log10_selected_k ~ log10(dimension) + log10(n_obs) + contam_rate, data = results),
		error = function(e) NULL
	)

	if (!is.null(csv_file)) {
		write_optimal_k_results_csv(results, csv_file, verbose = verbose, completed = TRUE)
	}

	list(
		results = results,
		summary = summary,
		scaling_fit = scaling_fit
	)
}

#' Default Size Grid for Optimal-k Testing
#'
#' @return A broad list of matrix dimension pairs for the default simulation sweep.
#' @keywords internal
default_optimal_k_size_grid <- function() {
	list(
		c(2, 3),
		c(3, 4),
		c(4, 6),
		c(5, 7),
		c(6, 8),
		c(7, 10),
		c(8, 10),
		c(10, 12),
		c(12, 16),
		c(14, 18),
		c(16, 20),
		c(18, 24),
		c(20, 24),
		c(24, 30),
		c(30, 36),
		c(36, 45)
	)
}

#' Simulate a Matrix Dataset for Optimal-k Testing
#'
#' @param r Number of rows in each matrix.
#' @param p Number of columns in each matrix.
#' @param group_sizes Integer vector of clean observations per Gaussian group.
#' @param n_contam Integer: number of contaminated matrices.
#' @param signal_shift Numeric: separation between clean component means.
#' @param row_sd Numeric: row-side noise scale.
#' @param col_sd Numeric: column-side noise scale.
#' @param contam_low Numeric: lower contamination bound.
#' @param contam_high Numeric: upper contamination bound.
#'
#' @return A list with the simulated matrices and the true labels.
#'
#' @keywords internal
simulate_optimal_k_dataset <- function(r, p, group_sizes, n_contam,
										signal_shift = 6,
										row_sd = 0.35,
										col_sd = 0.35,
										contam_low = -15,
										contam_high = 15) {
	n_groups <- length(group_sizes)
	clean_total <- sum(group_sizes)
	means <- lapply(seq_len(n_groups), function(component) {
		base_level <- (component - 1L) * signal_shift
		matrix(base_level, nrow = r, ncol = p)
	})

	row_cov <- diag(row_sd, r)
	col_cov <- diag(col_sd, p)

	clean_x <- vector("list", clean_total)
	labels <- integer(clean_total + n_contam)
	position <- 0L

	for (component in seq_len(n_groups)) {
		component_size <- group_sizes[component]
		for (i in seq_len(component_size)) {
			position <- position + 1L
			noise <- matrix(rnorm(r * p), r, p)
			clean_x[[position]] <- means[[component]] + row_cov %*% noise %*% col_cov
			labels[position] <- component
		}
	}

	contam_x <- if (n_contam > 0L) {
		lapply(seq_len(n_contam), function(i) {
			matrix(runif(r * p, min = contam_low, max = contam_high), r, p)
		})
	} else {
		list()
	}

	list(
		x_list = c(clean_x, contam_x),
		true_labels = c(labels, rep.int(0L, n_contam)),
		clean_labels = labels,
		contam_idx = if (n_contam > 0L) seq.int(clean_total + 1L, clean_total + n_contam) else integer(0)
	)
}

#' Summarize the Optimal-k Test Scheme
#'
#' @param results Data frame returned by [run_optimal_k_test_scheme()].
#'
#' @return A scenario-level summary data frame.
#'
#' @keywords internal
summarize_optimal_k_test_scheme <- function(results) {
	agg <- stats::aggregate(
		cbind(selected_k, log10_selected_k, selected_k_per_entry, log10_selected_k_per_entry,
		      ks_statistic, noise_rate) ~ r + p + dimension + n_obs + n_clean + n_contam + contam_rate,
		data = results,
		FUN = function(x) c(mean = mean(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
	)

	data.frame(
		r = agg$r,
		p = agg$p,
		dimension = agg$dimension,
		n_obs = agg$n_obs,
		n_clean = agg$n_clean,
		n_contam = agg$n_contam,
		contam_rate = agg$contam_rate,
		selected_k_mean = agg$selected_k[, "mean"],
		selected_k_median = agg$selected_k[, "median"],
		log10_selected_k_mean = agg$log10_selected_k[, "mean"],
		log10_selected_k_median = agg$log10_selected_k[, "median"],
		selected_k_per_entry_mean = agg$selected_k_per_entry[, "mean"],
		selected_k_per_entry_median = agg$selected_k_per_entry[, "median"],
		log10_selected_k_per_entry_mean = agg$log10_selected_k_per_entry[, "mean"],
		log10_selected_k_per_entry_median = agg$log10_selected_k_per_entry[, "median"],
		ks_statistic_mean = agg$ks_statistic[, "mean"],
		ks_statistic_median = agg$ks_statistic[, "median"],
		noise_rate_mean = agg$noise_rate[, "mean"],
		noise_rate_median = agg$noise_rate[, "median"],
		row.names = NULL,
		stringsAsFactors = FALSE
	)
}
