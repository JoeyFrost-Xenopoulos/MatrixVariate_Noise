#' Automatic K Selection for RIMLE
#'
#' Selects the optimal noise constant k by evaluating a grid of candidates and
#' scoring each with a Kolmogorov-Smirnov test on Mahalanobis distances of
#' non-noise observations.
#'
#' @param x_list List of matrices.
#' @param g Number of Gaussian components.
#' @param pi_max Maximum noise proportion.
#' @param gamma Eigenratio bound.
#' @param max_iter Maximum ECM iterations.
#' @param tol Convergence tolerance.
#' @param init Initialization method.
#' @param nstart Number of random restarts.
#' @param k_grid Grid of k values to search.
#' @param verbose Logical.
#' @param use_parallel Logical.
#' @param n_cores Number of parallel workers.
#' @return Fitted RIMLE model with k_selection diagnostics.
#' @noRd
rimle_select_k <- function(x_list, g, pi_max = 0.5, gamma = 1000,
			   max_iter = 100, tol = 1e-6,
			   init = c("random", "hennig-coretto", "kmeans"),
			   nstart = 10, k_grid = NULL,
			   verbose = FALSE, use_parallel = FALSE,
			   n_cores = NULL) {
	init <- match.arg(init)
	x_list <- rimle_validate_x_list(x_list)
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	if (is.null(k_grid)) {
		dimension <- r * p
		center_log10 <- -0.75 * dimension
		half_width <- max(6, ceiling(dimension / 2))
		lower_log10 <- max(log10(.Machine$double.xmin), center_log10 - half_width)
		upper_log10 <- center_log10 + half_width
		grid_log10 <- seq(lower_log10, upper_log10, length.out = 30)
		k_grid <- 10^grid_log10
		k_grid <- k_grid[is.finite(k_grid) & k_grid > 0]
		if (length(k_grid) < 2) {
			k_grid <- 10^seq(-16, -1, length.out = 30)
		}
	}

	ks_scores <- rep(Inf, length(k_grid))
	all_ks_results <- vector("list", length(k_grid))

	evaluate_k_candidate <- function(idx) {
		current_k <- k_grid[idx]

		if (verbose) {
			cat("  Testing k =", format(current_k, scientific = TRUE), "... ")
		}

		fit_noise <- tryCatch(
			rimle_fit_impl(x_list = x_list, g = g, k = current_k,
				       pi_max = pi_max, gamma = gamma,
				       max_iter = max_iter, tol = tol,
				       init = init, nstart = nstart,
				       use_parallel = use_parallel, n_cores = n_cores,
				       verbose = FALSE),
			error = function(e) {
				if (verbose) cat("fit failed\n")
				NULL
			}
		)

		if (is.null(fit_noise)) {
			return(list(statistic = Inf, p.value = NA_real_, n_used = 0))
		}

		keep_idx <- which(fit_noise$cluster > 0)
		x_clean <- x_list[keep_idx]

		if (length(x_clean) <= g) {
			if (verbose) cat("insufficient retained observations\n")
			return(list(statistic = Inf, p.value = NA_real_, n_used = length(x_clean)))
		}

		fit_clean <- tryCatch(
			mv_mixture_fit(x_list = x_clean, g = g, max_iter = max_iter, tol = tol, verbose = FALSE),
			error = function(e) {
				if (verbose) cat("refit failed\n")
				NULL
			}
		)

		if (is.null(fit_clean)) {
			return(list(statistic = Inf, p.value = NA_real_, n_used = length(x_clean)))
		}

		ks_result <- tryCatch(
			mv_noise_ks_score(fit_clean, x_clean),
			error = function(e) {
				if (verbose) {
					cat(sprintf("KS scoring failed: %s\n", conditionMessage(e)))
				}
				list(statistic = Inf, p.value = NA_real_, n_used = length(x_clean))
			}
		)

		if (verbose) {
			cat(sprintf("KS = %.4f (n_used = %d)\n", ks_result$statistic, ks_result$n_used))
		}

		ks_result
	}

	grid_results <- lapply(seq_along(k_grid), evaluate_k_candidate)

	for (i in seq_along(k_grid)) {
		ks_scores[i] <- grid_results[[i]]$statistic
		all_ks_results[[i]] <- grid_results[[i]]
	}

	if (all(is.infinite(ks_scores))) {
		stop("All candidate k values failed during KS selection.")
	}

	best_idx <- which.min(ks_scores)
	selected_k <- k_grid[best_idx]

	if (verbose) {
		cat(sprintf("\nSelected optimal k = %e (KS = %.4f)\n",
			    selected_k, ks_scores[best_idx]))
	}

	final_fit <- rimle_fit_impl(x_list = x_list, g = g, k = selected_k,
				    pi_max = pi_max, gamma = gamma,
				    max_iter = max_iter, tol = tol,
				    init = init, nstart = nstart,
				    use_parallel = use_parallel, n_cores = n_cores,
				    verbose = verbose)

	final_fit$k_selection <- list(
		selected_k = selected_k,
		k_grid = k_grid,
		ks_scores = ks_scores,
		ks_pvalues = sapply(all_ks_results, function(x) x$p.value),
		n_used = sapply(all_ks_results, function(x) x$n_used)
	)

	structure(final_fit, class = "rimle_fit")
}
