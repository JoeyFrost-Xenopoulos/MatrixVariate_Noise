#' Score HC Noise Fit with a Matrix KS Test
#'
#' @param fit A fitted noise model.
#' @param x_list List of matrices used for fitting.
#' @return A list with `statistic`, `p.value`, and `n_used`.
#' @keywords internal
matrix_noise_ks_score <- function(fit, x_list) {
	x_list <- matrix_validate_x_list(x_list)
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

#' Evaluate a Single HC Noise Candidate
#'
#' @param candidate_k Candidate noise height.
#' @param x_list List of matrices used for fitting.
#' @param g Number of Gaussian components.
#' @param max_iter Maximum EM iterations.
#' @param tol EM tolerance.
#' @param nstart K-means restarts.
#' @param noise_jitter Jitter passed through to the fitter.
#' @param noise_pi_init Initial noise mixing proportion.
#' @param verbose Logical: print progress messages.
#' @param keep_fit Logical: whether to retain the fitted model.
#' @return A list with the candidate summary and optional fitted model.
#' @keywords internal
matrix_noise_hc_candidate_eval <- function(candidate_k,
											   x_list,
											   g,
											   max_iter,
											   tol,
											   nstart,
											   noise_jitter,
											   noise_pi_init,
											   verbose,
											   keep_fit = TRUE) {
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
	candidate_loglik <- if (length(candidate_fit$logLik) > 0) tail(candidate_fit$logLik, 1) else NA_real_
	list(
		candidate_k = candidate_k,
		score = score,
		logLik = candidate_loglik,
		converged = candidate_fit$converged,
		fit = if (isTRUE(keep_fit)) candidate_fit else NULL
	)
}