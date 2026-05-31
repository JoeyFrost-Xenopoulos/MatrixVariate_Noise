#' Matrix-Variate Noise Mixture Clustering
#'
#' Fits a matrix-variate Gaussian mixture model with an additional noise
#' component. The noise component can be either Hennig-Coretto style improper
#' constant noise (`hc`) or Banfield-Raftery style bounded uniform noise (`br`).
#'
#' @param x_list A non-empty list of numeric matrices, each of the same size.
#' @param g Integer: number of Gaussian mixture components.
#' @param noise_type Character: `"hc"` for improper constant noise or `"br"`
#'   for convex-hull uniform noise.
#' @param max_iter Integer: maximum EM iterations.
#' @param tol Numeric: convergence tolerance on the log-likelihood trace.
#' @param nstart Integer: number of k-means restarts for initialization.
#' @param noise_k Numeric: constant noise height used when `noise_type = "hc"`.
#' @param select_noise_k Logical: if `TRUE` and `noise_type = "hc"`, grid
#'   search `noise_k_grid` and select the fit with the smallest matrix-KS
#'   statistic on fitted Mahalanobis distances.
#' @param noise_k_grid Numeric vector of candidate HC noise heights to try when
#'   `select_noise_k = TRUE`. The search also augments this grid with a
#'   matrix-dimension-aware log-scale heuristic so it can reach much smaller
#'   densities when the data dimension grows.
#' @param adaptive Logical: if `TRUE`, let the HC selector expand the grid when
#'   the best candidate lands on a boundary.
#' @param noise_jitter Numeric: retained for compatibility; unused by the BR
#'   convex-hull support.
#' @param noise_pi_init Numeric: initial mixing proportion for the noise
#'   component.
#' @param use_parallel Logical: if `TRUE`, evaluate HC candidate noise heights
#'   in parallel during selection.
#' @param n_cores Integer: number of worker processes to use when
#'   `use_parallel = TRUE`.
#' @param verbose Logical: print iteration progress.
#'
#' @return A list containing the fitted mixture parameters, posterior
#'   responsibilities, log-likelihood trace, and a noise summary. The returned
#'   responsibility matrix includes the noise component in the last column, and
#'   `cluster` uses `0` for observations assigned to noise.
#'
#' @details
#' The Gaussian components use the same matrix-variate normal updates as the
#' base fitter in [matrix_variate_mixture_fit()]. The noise component is fixed
#' during the M-step and only its mixing proportion is updated.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' x_list <- list(
#'   matrix(rnorm(12), 3, 4),
#'   matrix(rnorm(12, mean = 2), 3, 4),
#'   matrix(rnorm(12, mean = 8), 3, 4)
#' )
#'
#' fit_hc <- matrix_variate_noise_fit(x_list, g = 2, noise_type = "hc")
#' fit_br <- matrix_variate_noise_fit(x_list, g = 2, noise_type = "br")
#' }
#'
#' @export

matrix_variate_noise_fit <- function(x_list, g,
											 noise_type = c("hc", "br"),
											 max_iter = 1000,
											 tol = 1e-06,
											 nstart = 10,
											 noise_k = 1e-04,
											 select_noise_k = FALSE,
											 noise_k_grid = 10^seq(-8, -1, length.out = 15),
											 adaptive = TRUE,
											 noise_jitter = 1e-08,
											 noise_pi_init = 0.05,
											 use_parallel = FALSE,
											 n_cores = NULL,
											 verbose = FALSE) {
										
	noise_type <- match.arg(noise_type)

	x_list <- matrix_validate_x_list(x_list)

	if (noise_type == "hc" && isTRUE(select_noise_k)) {
		return(matrix_noise_hc_select_fit(
			x_list = x_list,
			g = g,
			max_iter = max_iter,
			tol = tol,
			nstart = nstart,
			noise_k_grid = noise_k_grid,
			adaptive = adaptive,
			noise_jitter = noise_jitter,
			noise_pi_init = noise_pi_init,
			use_parallel = use_parallel,
			n_cores = n_cores,
			verbose = verbose
		))
	}

	if (isTRUE(select_noise_k) && noise_type != "hc") {
		stop("select_noise_k is only supported when noise_type = \"hc\".")
	}

	fit <- matrix_variate_noise_fit_impl(
		x_list = x_list,
		g = g,
		noise_type = noise_type,
		max_iter = max_iter,
		tol = tol,
		nstart = nstart,
		noise_k = noise_k,
		noise_jitter = noise_jitter,
		noise_pi_init = noise_pi_init,
		verbose = verbose
	)

	fit$noise$search <- list(
		enabled = FALSE,
		adaptive = isTRUE(adaptive),
		criterion = "matrix_ks",
		selected_k = if (noise_type == "hc") noise_k else NA_real_
	)
	fit
}

matrix_variate_noise_fit_impl <- function(x_list, g,
											 noise_type = c("hc", "br"),
											 max_iter = 100,
											 tol = 1e-06,
											 nstart = 10,
											 noise_k = 1e-04,
											 noise_jitter = 1e-08,
											 noise_pi_init = 0.05,
											 verbose = FALSE) {
	noise_type <- match.arg(noise_type)

	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	for (x in x_list) {
		if (!is.matrix(x) || nrow(x) != r || ncol(x) != p) {
			stop("All elements of x_list must be matrices with the same dimensions.")
		}
	}

	# k-means init
	params <- matrix_mixture_kmeans_init(x_list, g = g, nstart = nstart)

	# For BR noise compute a convex hull over the vectorized matrices
	noise_support <- NULL
	if (noise_type == "br") {
		noise_support <- matrix_noise_convex_hull_support(x_list, jitter = noise_jitter)
	}

	# Append noise mixing proportion as the last component
	params$pi <- c((1 - noise_pi_init) * params$pi, noise_pi_init)
	names(params$pi) <- c(paste0("component_", seq_len(g)), "noise")

	loglik_trace <- numeric(0)
	responsibilities <- matrix(0, nrow = n, ncol = g + 1)
	colnames(responsibilities) <- c(paste0("component_", seq_len(g)), "noise")

	# Precompute noise log-density vector:
	# HC: constant improper background log(k)
	# BR: uniform within the convex hull (log(1/volume)), -Inf outside
	noise_log_density <- if (noise_type == "hc") {
		rep(log(noise_k), n)
	} else {
		matrix_noise_br_log_density(x_list, noise_support)
	}

	for (iteration in seq_len(max_iter)) {
		log_density <- matrix(NA_real_, nrow = n, ncol = g + 1)

		# E-step
		for (component in seq_len(g)) {
			for (i in seq_len(n)) {
				log_density[i, component] <- log(params$pi[component]) +
					matrix_variate_log_density(
						x = x_list[[i]],
						mean_matrix = params$M[[component]],
						row_cov = params$U[[component]],
						col_cov = params$V[[component]]
					)
			}
		}

		# Noise mixing proportion with noise log-density
		log_density[, g + 1] <- log(params$pi[g + 1]) + noise_log_density

		# Normalize log-densities to posterior responsibilities using log-sum-exp
		for (i in seq_len(n)) {
			row_log_densities <- log_density[i, ]
			normalizer <- matrix_log_sum_exp(row_log_densities)
			responsibilities[i, ] <- exp(row_log_densities - normalizer)
		}

		# Observed-data log-likelihood
		current_loglik <- sum(apply(log_density, 1, matrix_log_sum_exp))
		loglik_trace <- c(loglik_trace, current_loglik)

		if (iteration > 1 && abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
			break
		}

		# M-step: update Gaussian parameters using responsibilities
		component_responsibilities <- responsibilities[, seq_len(g), drop = FALSE]
		component_sizes <- colSums(component_responsibilities)
		noise_size <- sum(responsibilities[, g + 1])
		new_params <- params

		for (component in seq_len(g)) {
			if (component_sizes[component] <= 0) {
				next
			}

			# Effective weights for this component
			weights <- component_responsibilities[, component]
			weights_sum <- component_sizes[component]
			v_for_row <- make_spd(params$V[[component]])

			# Update mean matrix: weighted average
			mean_matrix <- matrix(0, r, p)
			for (i in seq_len(n)) {
				mean_matrix <- mean_matrix + weights[i] * x_list[[i]]
			}
			mean_matrix <- mean_matrix / weights_sum

			# Update row covariance U_g using current V_g (v_for_row)
			row_cov <- matrix(0, r, r)
			for (i in seq_len(n)) {
				centered <- x_list[[i]] - mean_matrix
				row_cov <- row_cov + weights[i] * (centered %*% solve(v_for_row, t(centered)))
			}
			row_cov <- row_cov / (p * weights_sum)
			row_cov <- make_spd(row_cov)
			
			# Identifiability
			row_scale <- r / sum(diag(row_cov))
			row_cov <- make_spd(row_cov * row_scale)

			# Update column covariance V_g using updated U_g
			col_cov <- matrix(0, p, p)
			for (i in seq_len(n)) {
				centered <- x_list[[i]] - mean_matrix
				col_cov <- col_cov + weights[i] * (t(centered) %*% solve(row_cov, centered))
			}
			col_cov <- col_cov / (r * weights_sum)
			col_cov <- make_spd(col_cov)

			# Store updated parameters for this Gaussian component
			new_params$pi[component] <- weights_sum / n
			new_params$M[[component]] <- mean_matrix
			new_params$U[[component]] <- row_cov
			new_params$V[[component]] <- col_cov
		}

		# Update noise mixing proportion
		new_params$pi[g + 1] <- noise_size / n
		new_params$pi <- new_params$pi / sum(new_params$pi)
		params <- new_params

		if (verbose) {
			if (noise_type == "hc") {
				message(sprintf("Iteration %d: log-likelihood = %.4f | noise_k = %.4e", iteration, current_loglik, noise_k))
			} else {
				message(sprintf("Iteration %d: log-likelihood = %.4f | noise_type = %s", iteration, current_loglik, noise_type))
			}
		}
	}

	# Hard assignments pick the component with maximum posterior; map noise -> 0
	cluster_membership <- max.col(responsibilities, ties.method = "first")
	cluster_membership[cluster_membership == g + 1] <- 0L

	list(
		pi = params$pi,
		M = params$M,
		U = params$U,
		V = params$V,
		z = responsibilities,
		cluster = cluster_membership,
		logLik = loglik_trace,
		iterations = length(loglik_trace),
		converged = length(loglik_trace) < max_iter,
		noise = list(
			type = noise_type,
			pi = params$pi[g + 1],
			k = if (noise_type == "hc") noise_k else NA_real_,
			hull = noise_support
		)
	)
}
