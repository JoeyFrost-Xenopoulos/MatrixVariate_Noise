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
                                       noise_k_grid = 10^seq(-8, -1, length.out = 15),
                                       noise_jitter = 1e-08,
                                       noise_pi_init = 0.05,
                                       verbose = FALSE) {
	sanitize_noise_grid <- function(values) {
		values <- as.numeric(values)
		zero_idx <- is.finite(values) & values == 0
		values[zero_idx] <- .Machine$double.xmin
		values <- values[is.finite(values) & values > 0]
		sort(unique(values))
	}

	candidate_grid <- matrix_noise_hc_search_grid(noise_k_grid = noise_k_grid, x_list = x_list)
	candidate_grid <- sanitize_noise_grid(candidate_grid)
	if (length(candidate_grid) == 0) {
		stop("HC noise_k selection failed: candidate grid has no finite positive values after sanitization.")
	}
	search_results <- list()
	best_fit <- NULL
	best_k <- NA_real_
	best_statistic <- Inf
	best_p_value <- NA_real_
	fallback_fit <- NULL
	fallback_k <- NA_real_
	fallback_n_used <- -Inf
	fallback_loglik <- -Inf
	seen_keys <- character(0)
	max_rounds <- 3L

	for (round_idx in seq_len(max_rounds)) {
		round_grid <- candidate_grid[!(format(candidate_grid, scientific = TRUE, digits = 16) %in% seen_keys)]
		if (length(round_grid) == 0) {
			break
		}

		round_best_k <- NA_real_
		round_best_statistic <- Inf
		round_best_p_value <- NA_real_
		round_fallback_k <- NA_real_
		round_fallback_n_used <- -Inf
		round_fallback_loglik <- -Inf

		for (candidate_k in round_grid) {
			candidate_key <- format(candidate_k, scientific = TRUE, digits = 16)
			seen_keys <- c(seen_keys, candidate_key)
			if (verbose) {
				message(sprintf("Checking HC noise_k = %.4e", candidate_k))
			}
			candidate_fit <- matrix_variate_noise_fit_impl(
				x_list = x_list,
				g = g,
				noise_type = "hc",
				max_iter = max_iter,
				tol = tol,
				nstart = nstart,
				noise_k = candidate_k,
				noise_jitter = noise_jitter,
				noise_pi_init = noise_pi_init,
				verbose = verbose
			)
			score <- matrix_noise_ks_score(candidate_fit, x_list)
			search_results[[length(search_results) + 1L]] <- data.frame(
				round = round_idx,
				noise_k = candidate_k,
				ks_statistic = score$statistic,
				ks_p_value = score$p.value,
				n_used = score$n_used,
				logLik = if (length(candidate_fit$logLik) > 0) tail(candidate_fit$logLik, 1) else NA_real_,
				converged = candidate_fit$converged,
				stringsAsFactors = FALSE
			)

			candidate_loglik <- if (length(candidate_fit$logLik) > 0) tail(candidate_fit$logLik, 1) else NA_real_
			if (is.finite(score$n_used) && (
				score$n_used > fallback_n_used ||
				(score$n_used == fallback_n_used && is.finite(candidate_loglik) && candidate_loglik > fallback_loglik)
			)) {
				fallback_fit <- candidate_fit
				fallback_k <- candidate_k
				fallback_n_used <- score$n_used
				fallback_loglik <- candidate_loglik
			}

			if (is.finite(score$n_used) && (
				score$n_used > round_fallback_n_used ||
				(score$n_used == round_fallback_n_used && is.finite(candidate_loglik) && candidate_loglik > round_fallback_loglik)
			)) {
				round_fallback_k <- candidate_k
				round_fallback_n_used <- score$n_used
				round_fallback_loglik <- candidate_loglik
			}

			if (is.finite(score$statistic) && (
				score$statistic < best_statistic ||
				(identical(score$statistic, best_statistic) && (!is.finite(best_p_value) || score$p.value > best_p_value))
			)) {
				best_fit <- candidate_fit
				best_k <- candidate_k
				best_statistic <- score$statistic
				best_p_value <- score$p.value
			}

			if (is.finite(score$statistic) && (
				score$statistic < round_best_statistic ||
				(identical(score$statistic, round_best_statistic) && (!is.finite(round_best_p_value) || score$p.value > round_best_p_value))
			)) {
				round_best_k <- candidate_k
				round_best_statistic <- score$statistic
				round_best_p_value <- score$p.value
			}
		}

		if (is.null(best_fit)) {
			best_fit <- fallback_fit
			best_k <- fallback_k
			best_statistic <- Inf
			best_p_value <- NA_real_
			if (!is.null(best_fit) && verbose) {
				message(sprintf(
					"HC noise_k selection fell back to the candidate with the most non-noise observations (n_used = %s).",
					format(fallback_n_used, trim = TRUE)
				))
			}
		}

		if (is.null(best_fit)) {
			stop("HC noise_k selection failed: no candidate fit could be retained.")
		}

		pivot_k <- if (is.finite(round_best_k)) round_best_k else round_fallback_k
		if (!is.finite(pivot_k) || pivot_k <= 0) {
			break
		}

		candidate_grid <- sanitize_noise_grid(candidate_grid)
		if (length(candidate_grid) == 0) {
			break
		}

		if (round_idx >= max_rounds) {
			break
		}

		log10_round <- log10(round_grid)
		log10_round <- log10_round[is.finite(log10_round)]
		round_span <- if (length(log10_round) > 0) diff(range(log10_round)) else NA_real_
		if (!is.finite(round_span) || round_span <= 0) {
			round_span <- 1
		}

		pivot_log10 <- log10(pivot_k)
		if (!is.finite(pivot_log10)) {
			break
		}

		# Refine search around the best region and shrink the window each round.
		refine_half_width <- max(round_span / 3, 0.5) / (2^(round_idx - 1))
		lower_log10 <- max(log10(.Machine$double.xmin), pivot_log10 - refine_half_width)
		upper_log10 <- pivot_log10 + refine_half_width
		if (!is.finite(lower_log10) || !is.finite(upper_log10) || lower_log10 >= upper_log10) {
			break
		}

		refine_count <- max(9L, min(21L, length(round_grid) + 2L))
		refined_grid <- 10^seq(lower_log10, upper_log10, length.out = refine_count)
		candidate_grid <- sanitize_noise_grid(c(candidate_grid, refined_grid))

		candidate_grid <- sanitize_noise_grid(candidate_grid)
	}

	best_fit$noise$search <- list(
		enabled = TRUE,
		criterion = "matrix_ks",
		grid = sort(unique(candidate_grid)),
		results = do.call(rbind, search_results),
		selected_k = best_k,
		ks_statistic = best_statistic,
		ks_p_value = best_p_value,
		fallback_used = is.na(best_statistic),
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
#'
#' @return A sorted unique numeric vector of positive candidate noise heights.
#' @keywords internal
matrix_noise_hc_search_grid <- function(noise_k_grid, x_list) {
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
