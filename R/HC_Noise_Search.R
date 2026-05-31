#' Select an HC Noise-k Fit
#'
#' @param x_list List of matrices used for fitting.
#' @param g Integer: number of Gaussian mixture components.
#' @param max_iter Integer: maximum EM iterations.
#' @param tol Numeric: EM convergence tolerance.
#' @param nstart Integer: number of k-means restarts.
#' @param noise_k_grid Numeric vector of candidate HC noise heights.
#' @param noise_jitter Numeric: retained for compatibility with the fitter.
#' @param noise_pi_init Numeric: initial noise mixing proportion.
#' @param verbose Logical: print progress messages.
#'
#' @return A fitted model with an added `fit$noise$search` summary.
#' @keywords internal
matrix_noise_hc_select_fit <- function(x_list, g,
									   max_iter = 100,
									   tol = 1e-06,
									   nstart = 10,
									   noise_k_grid = 10^seq(-8, -1, length.out = 50),
									   noise_jitter = 1e-08,
									   noise_pi_init = 0.05,
									   adaptive = TRUE,
									   verbose = FALSE,
									   use_parallel = FALSE,
									   n_cores = NULL) {
	x_list <- matrix_validate_x_list(x_list)

	sanitize_noise_grid <- function(values) {
		values <- as.numeric(values)
		zero_idx <- is.finite(values) & values == 0
		values[zero_idx] <- .Machine$double.xmin
		values <- values[is.finite(values) & values > 0]
		sort(unique(values))
	}

	choose_best_row <- function(results) {
		results <- results[is.finite(results$ks_statistic), , drop = FALSE]
		if (nrow(results) == 0L) {
			return(NULL)
		}
		p_value_rank <- results$ks_p_value
		p_value_rank[!is.finite(p_value_rank)] <- -Inf
		results <- results[order(results$ks_statistic, -p_value_rank, results$noise_k), , drop = FALSE]
		results[1L, , drop = FALSE]
	}

	choose_fallback_row <- function(results) {
		results <- results[is.finite(results$n_used), , drop = FALSE]
		if (nrow(results) == 0L) {
			return(NULL)
		}
		loglik_rank <- results$logLik
		loglik_rank[!is.finite(loglik_rank)] <- -Inf
		results <- results[order(-results$n_used, -loglik_rank, results$noise_k), , drop = FALSE]
		results[1L, , drop = FALSE]
	}

	expand_search_grid <- function(candidate_grid, selected_k) {
		candidate_grid <- sanitize_noise_grid(candidate_grid)
		if (length(candidate_grid) == 0L || !is.finite(selected_k) || selected_k <= 0) {
			return(candidate_grid)
		}

		selected_idx <- match(selected_k, candidate_grid)
		if (is.na(selected_idx)) {
			selected_idx <- which.min(abs(candidate_grid - selected_k))
		}

		lower_candidate <- NA_real_
		upper_candidate <- NA_real_
		if (length(candidate_grid) == 1L) {
			lower_candidate <- max(.Machine$double.xmin, candidate_grid[1L] / 10)
			upper_candidate <- candidate_grid[1L] * 10
		} else {
			if (selected_idx == 1L) {
				lower_ratio <- candidate_grid[2L] / candidate_grid[1L]
				if (!is.finite(lower_ratio) || lower_ratio <= 1) {
					lower_ratio <- 10
				}
				lower_candidate <- max(.Machine$double.xmin, candidate_grid[1L] / lower_ratio)
			}
			if (selected_idx == length(candidate_grid)) {
				upper_ratio <- candidate_grid[selected_idx] / candidate_grid[selected_idx - 1L]
				if (!is.finite(upper_ratio) || upper_ratio <= 1) {
					upper_ratio <- 10
				}
				upper_candidate <- candidate_grid[selected_idx] * upper_ratio
			}
		}

		sanitize_noise_grid(c(candidate_grid, lower_candidate, upper_candidate))
	}

	evaluate_candidate_grid <- function(candidate_grid, round_id) {
		search_results <- list()
		best_k <- NA_real_
		best_statistic <- Inf
		best_p_value <- NA_real_
		fallback_k <- NA_real_
		fallback_n_used <- -Inf
		fallback_loglik <- -Inf

		if (isTRUE(use_parallel) && length(candidate_grid) > 1L) {
			cores <- if (is.null(n_cores)) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
			cores <- max(1L, min(length(candidate_grid), cores))
			cl <- parallel::makeCluster(cores)
			on.exit(parallel::stopCluster(cl), add = TRUE)
			parallel::clusterExport(cl, varlist = c(
				"matrix_validate_x_list", "matrix_log_sum_exp",
				"make_spd", "matrix_mahalanobis", "matrix_variate_log_density",
				"matrix_mixture_kmeans_init", "matrix_variate_noise_fit_impl",
				"matrix_noise_ks_score", "matrix_noise_hc_candidate_eval"
			), envir = environment())
			candidate_results <- parallel::parLapply(
				cl,
				candidate_grid,
				matrix_noise_hc_candidate_eval,
				x_list = x_list,
				g = g,
				max_iter = max_iter,
				tol = tol,
				nstart = nstart,
				noise_jitter = noise_jitter,
				noise_pi_init = noise_pi_init,
				verbose = verbose,
				keep_fit = FALSE
			)
		} else {
			candidate_results <- lapply(candidate_grid, function(candidate_k) {
				if (verbose) message(sprintf("Checking HC noise_k = %.4e", candidate_k))
				matrix_noise_hc_candidate_eval(
					candidate_k = candidate_k,
					x_list = x_list,
					g = g,
					max_iter = max_iter,
					tol = tol,
					nstart = nstart,
					noise_jitter = noise_jitter,
					noise_pi_init = noise_pi_init,
					verbose = verbose,
					keep_fit = FALSE
				)
			})
		}

		for (res in candidate_results) {
			candidate_k <- res$candidate_k
			score <- res$score
			search_results[[length(search_results) + 1L]] <- data.frame(
				round = round_id,
				noise_k = candidate_k,
				ks_statistic = score$statistic,
				ks_p_value = score$p.value,
				n_used = score$n_used,
				logLik = res$logLik,
				converged = res$converged,
				stringsAsFactors = FALSE
			)
			if (is.finite(score$n_used) && (
				score$n_used > fallback_n_used ||
				(score$n_used == fallback_n_used && is.finite(res$logLik) && res$logLik > fallback_loglik)
			)) {
				fallback_k <- candidate_k
				fallback_n_used <- score$n_used
				fallback_loglik <- res$logLik
			}
			if (is.finite(score$statistic) && (
				score$statistic < best_statistic ||
				(identical(score$statistic, best_statistic) && (!is.finite(best_p_value) || score$p.value > best_p_value))
			)) {
				best_k <- candidate_k
				best_statistic <- score$statistic
				best_p_value <- score$p.value
			}
		}

		list(
			results = do.call(rbind, search_results),
			best_k = best_k,
			best_statistic = best_statistic,
			best_p_value = best_p_value,
			fallback_k = fallback_k,
			fallback_n_used = fallback_n_used,
			fallback_loglik = fallback_loglik
		)
	}

	candidate_grid <- matrix_noise_hc_search_grid(noise_k_grid = noise_k_grid, x_list = x_list, adaptive = adaptive)
	candidate_grid <- sanitize_noise_grid(candidate_grid)

	if (length(candidate_grid) == 0) {
		stop("HC noise_k selection failed: candidate grid has no finite positive values after sanitization.")
	}

	search_results <- list()
	best_fit <- NULL
	best_k <- NA_real_
	best_statistic <- Inf
	best_p_value <- NA_real_
	fallback_k <- NA_real_
	fallback_n_used <- -Inf
	fallback_loglik <- -Inf
	adaptive_rounds <- 0L
	adaptive_boundary_hit <- FALSE
	max_rounds <- if (isTRUE(adaptive)) 4L else 1L

	for (round_id in seq_len(max_rounds)) {
		round_eval <- evaluate_candidate_grid(candidate_grid, round_id)
		adaptive_rounds <- round_id
		search_results[[length(search_results) + 1L]] <- round_eval$results

		combined_results <- do.call(rbind, search_results)
		best_row <- choose_best_row(combined_results)
		fallback_row <- choose_fallback_row(combined_results)

		if (!is.null(best_row)) {
			best_k <- best_row$noise_k
			best_statistic <- best_row$ks_statistic
			best_p_value <- best_row$ks_p_value
		}
		if (!is.null(fallback_row)) {
			fallback_k <- fallback_row$noise_k
			fallback_n_used <- fallback_row$n_used
			fallback_loglik <- fallback_row$logLik
		}

		selected_k <- if (is.finite(best_k)) best_k else fallback_k
		if (!isTRUE(adaptive) || round_id >= max_rounds) {
			break
		}
		if (!is.finite(selected_k) || selected_k <= 0) {
			break
		}
		selected_idx <- match(selected_k, candidate_grid)
		if (is.na(selected_idx)) {
			selected_idx <- which.min(abs(candidate_grid - selected_k))
		}
		if (selected_idx %in% c(1L, length(candidate_grid))) {
			expanded_grid <- expand_search_grid(candidate_grid, selected_k)
			if (identical(expanded_grid, candidate_grid)) {
				break
			}
			adaptive_boundary_hit <- TRUE
			candidate_grid <- expanded_grid
			next
		}
		break
	}

	selected_k <- if (is.finite(best_k)) best_k else fallback_k
	if (!is.finite(selected_k) || selected_k <= 0) {
		stop("HC noise_k selection failed: no candidate fit could be retained.")
	}
	selected_result <- matrix_noise_hc_candidate_eval(
		candidate_k = selected_k,
		x_list = x_list,
		g = g,
		max_iter = max_iter,
		tol = tol,
		nstart = nstart,
		noise_jitter = noise_jitter,
		noise_pi_init = noise_pi_init,
		verbose = verbose,
		keep_fit = TRUE
	)
	best_fit <- selected_result$fit
	if (is.null(best_fit)) {
		stop("HC noise_k selection failed: no candidate fit could be retained.")
	}
	best_k <- selected_k
	best_statistic <- selected_result$score$statistic
	best_p_value <- selected_result$score$p.value

	combined_results <- do.call(rbind, search_results)
	best_fit$noise$search <- list(
		enabled = TRUE,
		adaptive = isTRUE(adaptive),
		adaptive_rounds = adaptive_rounds,
		adaptive_boundary_hit = adaptive_boundary_hit,
		criterion = "matrix_ks",
		grid = sort(unique(combined_results$noise_k)),
		results = combined_results,
		selected_k = best_k,
		ks_statistic = best_statistic,
		ks_p_value = best_p_value,
		fallback_used = !is.finite(best_statistic),
		fallback_n_used = if (is.finite(fallback_n_used)) fallback_n_used else NA_real_
	)
	best_fit
}

#' @keywords internal
matrix_noise_ks_score <- function(fit, x_list) {
	keep_idx <- which(fit$cluster > 0)
	if (length(keep_idx) < 2) {
		return(list(statistic = Inf, p.value = NA_real_, n_used = length(keep_idx)))
	}

	distances <- vapply(keep_idx, function(i) {
		component <- fit$cluster[i]
		matrix_mahalanobis(
			x = x_list[[i]],
			mean_matrix = fit$M[[component]],
			row_cov = fit$U[[component]],
			col_cov = fit$V[[component]]
		)
	}, numeric(1))
	distances <- distances[is.finite(distances)]
	if (length(distances) < 2 || length(unique(distances)) < 2) {
		return(list(statistic = Inf, p.value = NA_real_, n_used = length(distances)))
	}

	df <- nrow(x_list[[1]]) * ncol(x_list[[1]])
	test <- tryCatch(
		suppressWarnings(stats::ks.test(distances, "pchisq", df = df)),
		error = function(e) NULL
	)
	if (is.null(test)) {
		return(list(statistic = Inf, p.value = NA_real_, n_used = length(distances)))
	}

	list(
		statistic = unname(test$statistic),
		p.value = unname(test$p.value),
		n_used = length(distances)
	)
}

#' Build an HC Noise Search Grid
#'
#' @param noise_k_grid Candidate HC noise heights supplied by the caller.
#' @param x_list List of matrices used to infer a dimension-aware heuristic.
#' @param adaptive Logical: if `TRUE`, augment the grid with a dimension-aware
#'   heuristic before selection.
#'
#' @return A sorted unique numeric vector of positive candidate noise heights.
#' @keywords internal
matrix_noise_hc_search_grid <- function(noise_k_grid, x_list, adaptive = TRUE) {
	x_list <- matrix_validate_x_list(x_list)

	sanitize_noise_grid <- function(values) {
		values <- as.numeric(values)
		zero_idx <- is.finite(values) & values == 0
		values[zero_idx] <- .Machine$double.xmin
		values <- values[is.finite(values) & values > 0]
		sort(unique(values))
	}

	candidate_grid <- sanitize_noise_grid(noise_k_grid)
	if (length(candidate_grid) == 0) {
		stop("noise_k_grid must contain at least one positive finite value.")
	}
	if (!isTRUE(adaptive)) {
		return(candidate_grid)
	}

	dimension <- nrow(x_list[[1]]) * ncol(x_list[[1]])
	if (is.finite(dimension) && dimension > 0) {
		center_log10 <- -0.75 * dimension
		half_width <- max(6, ceiling(dimension / 2))
		lower_log10 <- max(log10(.Machine$double.xmin), center_log10 - half_width)
		upper_log10 <- center_log10 + half_width
		if (is.finite(lower_log10) && is.finite(upper_log10)) {
			heuristic_grid <- 10^seq(lower_log10, upper_log10, length.out = max(9L, length(candidate_grid)))
			candidate_grid <- sanitize_noise_grid(c(candidate_grid, heuristic_grid))
		}
	}

	sanitize_noise_grid(candidate_grid)
}
