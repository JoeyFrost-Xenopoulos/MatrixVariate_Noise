#' Core RIMLE Fit Implementation
#'
#' Internal ECM loop for the RIMLE matrix-variate mixture model.
#'
#' @param x_list List of matrices.
#' @param g Number of Gaussian components.
#' @param k Noise constant density.
#' @param pi_max Maximum noise proportion.
#' @param gamma Eigenratio bound.
#' @param max_iter Maximum ECM iterations.
#' @param tol Convergence tolerance.
#' @param init Initialization method.
#' @param nstart Number of random restarts.
#' @param verbose Logical.
#' @return Parameter list without class.
#' @noRd
rimle_fit_impl <- function(x_list, g, k, pi_max = 0.5, gamma = 1000,
			   max_iter = 100, tol = 1e-6,
			   init = c("random", "hennig-coretto", "kmeans"),
			   nstart = 10, verbose = FALSE) {
	init <- match.arg(init)
	x_list <- rimle_validate_x_list(x_list)
	n <- length(x_list)

	best_params <- NULL
	best_loglik <- -Inf

	nstart <- as.integer(nstart)

	for (start in seq_len(nstart)) {
		if (init == "kmeans") {
			params <- rimle_kmeans_init(x_list, g, nstart = 1, verbose = verbose)
			params$pi0 <- max(0.01, 1 - sum(params$pi))
		} else if (init == "random") {
			params <- rimle_random_init(x_list, g)
			params$pi0 <- max(0.01, 1 - sum(params$pi))
		} else {
			params <- rimle_hennig_coretto_init(x_list, g, pi_max = pi_max)
			if (is.null(params$pi0)) params$pi0 <- max(0.01, 1 - sum(params$pi))
		}

		if (params$pi0 + sum(params$pi) > 1) {
			total <- params$pi0 + sum(params$pi)
			params$pi0 <- params$pi0 / total
			params$pi <- params$pi / total
		}

		params$z <- matrix(0, n, g + 1)

		loglik_trace <- numeric(0)
		converged <- FALSE

		for (iteration in seq_len(max_iter)) {
			old_params <- params

			params <- rimle_e_step(x_list, params, g, n, k)
			params <- rimle_cm1_step(x_list, params, g, n, gamma)
			params <- rimle_cm2_step(x_list, params, g, n, k, pi_max)

			z0 <- pmax(params$z[, g + 1], .Machine$double.eps)
			current_loglik <- sum(log(params$pi0 * k) - log(z0))

			if (!is.finite(current_loglik)) current_loglik <- -1e10

			loglik_trace <- c(loglik_trace, current_loglik)

			if (iteration > 1 &&
			    abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) /
			    (abs(loglik_trace[iteration - 1]) + .Machine$double.eps) < tol) {
				converged <- TRUE
				break
			}
		}

		if (verbose && nstart > 1) {
			message(sprintf("Start %d: log-likelihood = %.4f", start, tail(loglik_trace, 1)))
		}

		if (tail(loglik_trace, 1) > best_loglik) {
			best_loglik <- tail(loglik_trace, 1)
			best_params <- params
			best_params$logLik <- loglik_trace
			best_params$iterations <- length(loglik_trace)
			best_params$converged <- converged
			best_params$k <- k
			best_params$pi_max <- pi_max
			best_params$gamma <- gamma
			best_params$g <- g
			best_params$init <- init
		}
	}

	cluster_membership <- max.col(best_params$z, ties.method = "first")
	cluster_membership[cluster_membership == g + 1] <- 0L

	best_params$cluster <- cluster_membership
	best_params$call <- match.call()

	best_params
}

#' Fit RIMLE Matrix-Variate Mixture Model
#'
#' Fits a robust improper mixture log-likelihood estimator (RIMLE) to a list
#' of matrix-valued observations. Supports automatic selection of the noise
#' constant k via a Kolmogorov-Smirnov goodness-of-fit test.
#'
#' @param x_list A non-empty list of numeric matrices, each of the same size.
#' @param g Integer: number of Gaussian mixture components.
#' @param k Numeric: constant noise height (improper uniform density value).
#' @param pi_max Maximum proportion of observations assignable to noise
#'   (default: 0.5).
#' @param gamma Eigenratio bound for covariance regularization (default: 1000).
#' @param max_iter Maximum ECM iterations (default: 100).
#' @param tol Convergence tolerance on the pseudo-log-likelihood trace.
#' @param init Initialization scheme: `"random"`, `"hennig-coretto"`, or `"kmeans"`.
#' @param nstart Number of independent random starts (default: 10).
#' @param estimate_k Logical: if TRUE, automatically select optimal k using
#'   KS goodness-of-fit test.
#' @param k_grid Numeric vector: grid of k values for automatic selection.
#' @param verbose Logical: print iteration progress.
#' @param use_parallel Logical: run k-grid search in parallel.
#' @param n_cores Integer: number of parallel workers (NULL = auto).
#'
#' @return A list of class `"rimle_fit"` containing fitted parameters,
#'   posterior responsibilities, cluster assignments, and diagnostics.
#'
#' @export
rimle_fit <- function(x_list, g, k, pi_max = 0.5, gamma = 1000,
		      max_iter = 100, tol = 1e-6,
		      init = c("random", "hennig-coretto", "kmeans"),
		      nstart = 10, estimate_k = FALSE,
		      k_grid = NULL, verbose = FALSE,
		      use_parallel = FALSE, n_cores = NULL) {
	init <- match.arg(init)
	x_list <- rimle_validate_x_list(x_list)
	n <- length(x_list)

	if (!is.numeric(g) || length(g) != 1 || g < 1) {
		stop("'g' must be a positive integer.")
	}
	g <- as.integer(g)

	if (n < g) {
		stop(sprintf("Number of observations (%d) must be at least as large as number of components (%d).",
			     n, g))
	}

	if (estimate_k) {
		return(rimle_select_k(x_list, g, pi_max = pi_max, gamma = gamma,
				      max_iter = max_iter, tol = tol,
				      init = init, nstart = nstart,
				      k_grid = k_grid, verbose = verbose,
				      use_parallel = use_parallel, n_cores = n_cores))
	}

	if (!is.numeric(k) || length(k) != 1 || k <= 0) {
		stop("'k' must be a positive numeric scalar.")
	}

	if (!is.numeric(pi_max) || length(pi_max) != 1 || pi_max <= 0 || pi_max >= 1) {
		stop("'pi_max' must be in (0, 1).")
	}

	if (!is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
		stop("'gamma' must be positive.")
	}

	params <- rimle_fit_impl(x_list, g, k, pi_max = pi_max, gamma = gamma,
				 max_iter = max_iter, tol = tol,
				 init = init, nstart = nstart,
				 verbose = verbose)

	structure(params, class = "rimle_fit")
}

#' Print Method for RIMLE Fits
#'
#' @param x A fitted `rimle_fit` object.
#' @param ... Additional arguments (unused).
#' @return `x` invisibly.
#' @export
print.rimle_fit <- function(x, ...) {
	if (!inherits(x, "rimle_fit")) {
		stop("'x' must be a fitted 'rimle_fit' model.")
	}

	n <- nrow(x$z)
	r <- nrow(x$M[[1]])
	p <- ncol(x$M[[1]])
	g <- x$g
	sizes <- tabulate(x$cluster[x$cluster > 0], nbins = g)
	noise_count <- sum(x$cluster == 0)

	cat("Robust Improper Mixture Log-Likelihood Estimator (RIMLE)\n")
	cat("========================================================\n\n")
	cat(sprintf("Observations: %d\n", n))
	cat(sprintf("Dimensions:   %d x %d\n", r, p))
	cat(sprintf("Components:   %d Gaussian + 1 noise\n\n", g))

	cat("Mixing Proportions:\n")
	for (i in seq_len(g)) {
		cat(sprintf("  component_%d: %.4f\n", i, x$pi[i]))
	}
	cat(sprintf("  noise:       %.4f\n", x$pi0))

	cat("\nComponent Sizes:\n")
	for (i in seq_len(g)) {
		count <- sizes[i]
		pct <- count / n * 100
		cat(sprintf("  component_%d: %d (%.1f%%)\n", i, count, pct))
	}
	noise_pct <- noise_count / n * 100
	cat(sprintf("  noise:       %d (%.1f%%)\n", noise_count, noise_pct))

	cat("\nConstraints:\n")
	cat(sprintf("  pi_max: %.4f\n", x$pi_max))
	cat(sprintf("  gamma:  %.4f\n", x$gamma))
	cat(sprintf("  k:      %.4e\n", x$k))

	cat("\nConvergence:\n")
	cat(sprintf("  Iterations:     %d\n", x$iterations))
	cat(sprintf("  Converged:      %s\n", if (x$converged) "Yes" else "No"))
	cat(sprintf("  Log-likelihood: %.4f\n", tail(x$logLik, 1)))
	cat(sprintf("  Init:           %s\n", x$init))

	invisible(x)
}

#' Summary Method for RIMLE Fits
#'
#' @param object A fitted `rimle_fit` object.
#' @param ... Additional arguments (unused).
#' @return A list of class `"summary.rimle_fit"`.
#' @export
summary.rimle_fit <- function(object, ...) {
	if (!inherits(object, "rimle_fit")) {
		stop("'object' must be a fitted 'rimle_fit' model.")
	}

	n <- nrow(object$z)
	r <- nrow(object$M[[1]])
	p <- ncol(object$M[[1]])
	g <- object$g
	sizes <- tabulate(object$cluster[object$cluster > 0], nbins = g)
	names(sizes) <- paste0("component_", seq_len(g))
	noise_count <- sum(object$cluster == 0)
	component_sizes <- c(sizes, noise = noise_count)

	n_free <- g * (r * p + r * (r + 1) / 2 + p * (p + 1) / 2 + 1)
	logLik_final <- tail(object$logLik, 1)
	aic <- -2 * logLik_final + 2 * n_free
	bic <- -2 * logLik_final + log(n) * n_free

	mixing_proportions <- c(object$pi, noise = object$pi0)

	res <- list(
		n_obs = n,
		dimensions = c(r = r, p = p),
		n_components = g,
		mixing_proportions = mixing_proportions,
		component_sizes = component_sizes,
		converged = object$converged,
		iterations = object$iterations,
		logLik = logLik_final,
		aic = aic,
		bic = bic,
		n_free_params = n_free,
		pi_max = object$pi_max,
		gamma = object$gamma,
		k = object$k,
		init = object$init
	)

	if (!is.null(object$k_selection)) {
		res$k_selection <- object$k_selection
	}

	structure(res, class = "summary.rimle_fit")
}

#' Print Summary Method for RIMLE Fits
#'
#' @param x A `summary.rimle_fit` object.
#' @param ... Additional arguments (unused).
#' @return `x` invisibly.
#' @export
print.summary.rimle_fit <- function(x, ...) {
	cat("Summary of Robust Improper Mixture Log-Likelihood Estimator (RIMLE)\n")
	cat("==================================================================\n\n")
	cat(sprintf("Observations:      %d\n", x$n_obs))
	cat(sprintf("Dimensions:        %d x %d\n", x$dimensions["r"], x$dimensions["p"]))
	cat(sprintf("Components:        %d Gaussian + 1 noise\n", x$n_components))
	cat(sprintf("Free parameters:   %d\n\n", x$n_free_params))

	cat("Mixing Proportions:\n")
	for (i in seq_along(x$mixing_proportions)) {
		cat(sprintf("  %s: %.4f\n", names(x$mixing_proportions)[i], x$mixing_proportions[i]))
	}

	cat("\nComponent Sizes:\n")
	for (i in seq_along(x$component_sizes)) {
		pct <- x$component_sizes[i] / x$n_obs * 100
		cat(sprintf("  %s: %d (%.1f%%)\n", names(x$component_sizes)[i], x$component_sizes[i], pct))
	}

	cat("\nConvergence:\n")
	cat(sprintf("  Iterations:      %d\n", x$iterations))
	cat(sprintf("  Converged:       %s\n", if (x$converged) "Yes" else "No"))
	cat(sprintf("  Log-likelihood:  %.4f\n", x$logLik))
	cat(sprintf("  AIC:             %.4f\n", x$aic))
	cat(sprintf("  BIC:             %.4f\n", x$bic))

	cat("\nConstraints:\n")
	cat(sprintf("  pi_max: %.4f\n", x$pi_max))
	cat(sprintf("  gamma:  %.4f\n", x$gamma))
	cat(sprintf("  k:      %.4e\n", x$k))
	cat(sprintf("  init:   %s\n", x$init))

	if (!is.null(x$k_selection)) {
		cat(sprintf("\nAutomatic K Selection:\n"))
		cat(sprintf("  Selected k: %.4e\n", x$k_selection$selected_k))
		cat(sprintf("  Grid size:  %d\n", length(x$k_selection$k_grid)))
	}

	invisible(x)
}
