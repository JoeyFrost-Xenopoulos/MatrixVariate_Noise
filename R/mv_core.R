#' Validate Matrix List Input
#'
#' @param x_list List of matrices to validate.
#' @return A list of same-sized matrices.
#' @noRd
mv_validate_x_list <- function(x_list) {
	if (!is.list(x_list) || length(x_list) == 0) {
		stop("x_list must be a non-empty list of matrices.")
	}

	if (!is.matrix(x_list[[1]])) {
		stop("First element of x_list is not a matrix.")
	}

	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	if (r == 0 || p == 0) {
		stop("Matrices in x_list must have at least one row and one column.")
	}

	for (idx in seq_along(x_list)) {
		x <- x_list[[idx]]
		if (!is.matrix(x)) {
			stop(sprintf("Element %d of x_list is not a matrix.", idx))
		}
		if (nrow(x) != r || ncol(x) != p) {
			stop(sprintf(
				"Element %d of x_list has dimensions %d x %d, expected %d x %d.",
				idx, nrow(x), ncol(x), r, p
			))
		}
	}

	x_list
}

#' Stable Log-Sum-Exp
#'
#' @param values Numeric vector.
#' @return Numeric scalar.
#' @noRd
mv_log_sum_exp <- function(values) {
	if (!is.numeric(values)) {
		stop("'values' must be a numeric vector.")
	}
	finite_values <- values[is.finite(values)]
	if (length(finite_values) == 0) {
		return(-Inf)
	}
	max_value <- max(finite_values)
	max_value + log(sum(exp(finite_values - max_value)))
}

#' Compute Matrix-Variate Mahalanobis Distance
#'
#' Calculates the Mahalanobis distance between a matrix and a mean matrix
#' under the matrix-variate normal distribution with specified row and column
#' covariance structures.
#'
#' @param x A numeric matrix (r x p): the observation
#' @param mean_matrix A numeric matrix (r x p): the component mean
#' @param row_cov A numeric matrix (r x r): row covariance matrix U
#' @param col_cov A numeric matrix (p x p): column covariance matrix V
#'
#' @return Numeric scalar representing the Mahalanobis distance
#'
#' @details
#'
#' This metric extends the multivariate Mahalanobis distance to account for
#' the matrix structure. The computation uses Cholesky decomposition and
#' forward/backsolve for numerical stability.
#'
#' @noRd
mv_mahalanobis <- function(x, mean_matrix, row_cov, col_cov) {
	if (!is.matrix(x) || !is.numeric(x)) {
		stop("'x' must be a numeric matrix.")
	}
	if (!is.matrix(mean_matrix) || !is.numeric(mean_matrix)) {
		stop("'mean_matrix' must be a numeric matrix.")
	}
	if (!identical(dim(x), dim(mean_matrix))) {
		stop(sprintf(
			"'x' (%d x %d) and 'mean_matrix' (%d x %d) must have the same dimensions.",
			nrow(x), ncol(x), nrow(mean_matrix), ncol(mean_matrix)
		))
	}
	if (!is.matrix(row_cov) || nrow(row_cov) != ncol(row_cov) || nrow(row_cov) != nrow(x)) {
		stop("'row_cov' must be a square matrix with dimension matching nrow(x).")
	}
	if (!is.matrix(col_cov) || nrow(col_cov) != ncol(col_cov) || nrow(col_cov) != ncol(x)) {
		stop("'col_cov' must be a square matrix with dimension matching ncol(x).")
	}
	row_cov <- make_spd(row_cov)
	col_cov <- make_spd(col_cov)
	row_chol <- chol(row_cov)
	col_chol <- chol(col_cov)
	centered <- x - mean_matrix
	row_inv_centered <- backsolve(row_chol, forwardsolve(t(row_chol), centered))
	col_inv <- chol2inv(col_chol)

	sum(row_inv_centered * (centered %*% col_inv))
}

#' Compute Log-Likelihood of Matrix under Matrix-Variate Normal Distribution
#'
#' Evaluates the log-density of a matrix observation under the matrix-variate
#' normal distribution with specified parameters.
#'
#' @param x A numeric matrix (r x p): the observation
#' @param mean_matrix A numeric matrix (r x p): the component mean matrix M
#' @param row_cov A numeric matrix (r x r): row covariance matrix U
#' @param col_cov A numeric matrix (p x p): column covariance matrix V
#'
#' @return Numeric scalar representing the log-density value
#'
#' @details
#'
#' Computation uses Cholesky decomposition for numerical stability and to
#' avoid explicit matrix inversion.
#'
#' @noRd
mv_log_density <- function(x, mean_matrix, row_cov, col_cov) {
	if (!is.matrix(x) || !is.numeric(x)) {
		stop("'x' must be a numeric matrix.")
	}
	if (!is.matrix(mean_matrix) || !is.numeric(mean_matrix)) {
		stop("'mean_matrix' must be a numeric matrix.")
	}
	if (!identical(dim(x), dim(mean_matrix))) {
		stop("'x' and 'mean_matrix' must have the same dimensions.")
	}
	if (!is.matrix(row_cov) || nrow(row_cov) != ncol(row_cov) || nrow(row_cov) != nrow(x)) {
		stop("'row_cov' must be a square matrix with dimension matching nrow(x).")
	}
	if (!is.matrix(col_cov) || nrow(col_cov) != ncol(col_cov) || nrow(col_cov) != ncol(x)) {
		stop("'col_cov' must be a square matrix with dimension matching ncol(x).")
	}
	row_cov <- make_spd(row_cov)
	col_cov <- make_spd(col_cov)
	row_chol <- chol(row_cov)
	col_chol <- chol(col_cov)

	row_logdet <- 2 * sum(log(diag(row_chol)))
	col_logdet <- 2 * sum(log(diag(col_chol)))

	centered <- x - mean_matrix
	row_inv_centered <- backsolve(row_chol, forwardsolve(t(row_chol), centered))
	col_inv <- chol2inv(col_chol)
	trace_form <- sum(row_inv_centered * (centered %*% col_inv))

	r <- nrow(x)
	p <- ncol(x)

	-0.5 * (r * p * log(2 * pi) + p * row_logdet + r * col_logdet + trace_form)
}

#' Compute E-Step Log-Densities for Gaussian Components
#'
#' Evaluates log(pi_g) + log f(X_i | theta_g) for each observation and component.
#'
#' @param x_list List of matrices.
#' @param params Parameter list with pi, M, U, V.
#' @param g Number of Gaussian components.
#' @param n Number of observations.
#' @return A matrix (n x g) of weighted log-densities for Gaussian components.
#' @noRd
mv_e_step_log_density <- function(x_list, params, g, n) {
	log_density <- matrix(NA_real_, nrow = n, ncol = g)
	for (component in seq_len(g)) {
		for (i in seq_len(n)) {
			log_density[i, component] <- log(params$pi[component]) +
				mv_log_density(
					x = x_list[[i]],
					mean_matrix = params$M[[component]],
					row_cov = params$U[[component]],
					col_cov = params$V[[component]]
				)
		}
	}
	log_density
}

#' Normalize Log-Densities to Posterior Responsibilities
#'
#' Applies row-wise log-sum-exp normalization.
#'
#' @param log_density Matrix of log-densities (n x K).
#' @return Matrix of posterior responsibilities (n x K), rows sum to 1.
#' @noRd
mv_normalize_responsibilities <- function(log_density) {
	n <- nrow(log_density)
	responsibilities <- matrix(0, nrow = n, ncol = ncol(log_density))
	for (i in seq_len(n)) {
		normalizer <- mv_log_sum_exp(log_density[i, ])
		if (!is.finite(normalizer)) {
			responsibilities[i, ] <- 1 / ncol(log_density)
		} else {
			responsibilities[i, ] <- exp(log_density[i, ] - normalizer)
			responsibilities[i, ][!is.finite(responsibilities[i, ])] <- 0
		}
	}
	responsibilities
}

#' Observed-Data Log-Likelihood from Log-Densities
#'
#' Sums the row-wise log-sum-exp of a log-density matrix.
#'
#' @param log_density Matrix of log-densities (n x K).
#' @return Numeric scalar: the observed-data log-likelihood.
#' @noRd
mv_loglik <- function(log_density) {
	sum(apply(log_density, 1, mv_log_sum_exp))
}

#' Compute Weighted Mean Matrix (M-Step)
#'
#' @param x_list List of matrices.
#' @param weights Numeric vector of responsibilities for one component.
#' @param weights_sum Sum of weights (effective sample size).
#' @param r Number of rows.
#' @param p Number of columns.
#' @return Weighted mean matrix (r x p).
#' @noRd
mv_weighted_mean <- function(x_list, weights, weights_sum, r, p) {
	n <- length(x_list)
	mean_matrix <- matrix(0, r, p)
	for (i in seq_len(n)) {
		mean_matrix <- mean_matrix + weights[i] * x_list[[i]]
	}
	mean_matrix / weights_sum
}

#' Update Row Covariance (M-Step)
#'
#' Computes the row covariance update: U = (1/(p*n_g)) sum w_i (X_i-M) V^{-1} (X_i-M)^T
#'
#' @param x_list List of matrices.
#' @param mean_matrix Current mean matrix.
#' @param v_inv_target Column covariance used for the update (made SPD internally).
#' @param weights Numeric vector of responsibilities.
#' @param weights_sum Effective sample size.
#' @param r Number of rows.
#' @param p Number of columns.
#' @param scale_trace Logical: if TRUE, enforce tr(U) = r identifiability constraint.
#' @return Updated row covariance matrix (r x r), positive definite.
#' @noRd
mv_update_row_cov <- function(x_list, mean_matrix, v_inv_target, weights,
                                  weights_sum, r, p, scale_trace = TRUE) {
	n <- length(x_list)
	v_spd <- make_spd(v_inv_target)
	row_cov <- matrix(0, r, r)
	for (i in seq_len(n)) {
		centered <- x_list[[i]] - mean_matrix
		row_cov <- row_cov + weights[i] * (centered %*% solve(v_spd, t(centered)))
	}
	row_cov <- row_cov / (p * weights_sum)
	row_cov <- make_spd(row_cov)

	if (scale_trace) {
		row_scale <- r / sum(diag(row_cov))
		row_cov <- make_spd(row_cov * row_scale)
	}
	row_cov
}

#' Update Column Covariance (M-Step)
#'
#' Computes the column covariance update: V = (1/(r*n_g)) sum w_i (X_i-M)^T U^{-1} (X_i-M)
#'
#' @param x_list List of matrices.
#' @param mean_matrix Current mean matrix.
#' @param u_inv_target Row covariance used for the update.
#' @param weights Numeric vector of responsibilities.
#' @param weights_sum Effective sample size.
#' @param r Number of rows.
#' @param p Number of columns.
#' @return Updated column covariance matrix (p x p), positive definite.
#' @noRd
mv_update_col_cov <- function(x_list, mean_matrix, u_inv_target, weights,
                                  weights_sum, r, p) {
	n <- length(x_list)
	col_cov <- matrix(0, p, p)
	for (i in seq_len(n)) {
		centered <- x_list[[i]] - mean_matrix
		col_cov <- col_cov + weights[i] * (t(centered) %*% solve(u_inv_target, centered))
	}
	col_cov <- col_cov / (r * weights_sum)
	make_spd(col_cov)
}

#' Shared EM M-Step
#'
#' Updates the parameters of the first `g` Gaussian components from posterior
#' responsibilities. Shared by `mv_mixture_fit`, `mv_noise_fit_impl`,
#' `mv_short_em_burn_in`, and `mv_mixture_emrefine_init` so the per-component
#' update logic lives in exactly one place.
#'
#' @param params Current parameter list (pi, M, U, V).
#' @param x_list Validated list of matrices.
#' @param responsibilities Posterior responsibility matrix (n x g) for the
#'   Gaussian components only.
#' @param g Number of Gaussian components.
#' @param n Number of observations.
#' @param r Number of rows.
#' @param p Number of columns.
#' @param warn_zero Logical: emit a warning when a component has zero effective
#'   membership (used by full fits and emrefine, not by burn-in).
#' @return A new parameter list with component pi, M, U, V updated. The global
#'   mixing proportions are NOT renormalized here; callers renormalize over the
#'   components they manage (g, or g + 1 for the noise model).
#' @noRd
mv_em_mstep <- function(params, x_list, responsibilities, g, n, r, p,
                          warn_zero = FALSE) {
	component_sizes <- colSums(responsibilities)
	new_params <- params

	for (component in seq_len(g)) {
		if (component_sizes[component] <= 0 || !is.finite(component_sizes[component])) {
			if (warn_zero) {
				warning(
					sprintf(
						"Component %d has zero effective membership; skipping update.",
						component
					),
					call. = FALSE
				)
			}
			next
		}

		weights <- responsibilities[, component]
		weights_sum <- component_sizes[component]

		mean_matrix <- mv_weighted_mean(x_list, weights, weights_sum, r, p)
		row_cov <- mv_update_row_cov(x_list, mean_matrix, params$V[[component]],
		                              weights, weights_sum, r, p)
		col_cov <- mv_update_col_cov(x_list, mean_matrix, row_cov,
		                              weights, weights_sum, r, p)

		new_params$pi[component] <- weights_sum / n
		new_params$M[[component]] <- mean_matrix
		new_params$U[[component]] <- row_cov
		new_params$V[[component]] <- col_cov
	}

	new_params
}

#' Compute Component Parameters from Cluster Assignments
#'
#' Shared logic used by both kmeans and emrefine initialization to compute
#' mean matrices and covariances from initial cluster assignments.
#'
#' @param x_list Validated list of matrices.
#' @param g Number of components.
#' @param cluster_assignments Integer vector of cluster labels (1..g).
#' @param init_method Character label for warnings (e.g. "K-means", "Random").
#' @return A list with pi, M, U, V, cluster.
#' @noRd
mv_compute_init_params <- function(x_list, g, cluster_assignments, init_method = "Initialization") {
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	mixing_proportions <- numeric(g)
	mean_matrices <- vector("list", g)
	row_covariances <- vector("list", g)
	col_covariances <- vector("list", g)

	for (component in seq_len(g)) {
		component_index <- which(cluster_assignments == component)
		if (length(component_index) == 0) {
			warning(sprintf(
				"%s initialization: component %d received no observations; using a deterministic observation as seed.",
				init_method, component
			), call. = FALSE)
			component_index <- ((component - 1L) %% n) + 1L
		}

		component_data <- x_list[component_index]
		mixing_proportions[component] <- length(component_index) / n
		mean_matrices[[component]] <- Reduce(`+`, component_data) / length(component_data)

		row_cov <- matrix(0, r, r)
		col_cov <- matrix(0, p, p)
		for (x in component_data) {
			centered <- x - mean_matrices[[component]]
			row_cov <- row_cov + centered %*% t(centered)
			col_cov <- col_cov + t(centered) %*% centered
		}

		row_cov <- row_cov / (p * length(component_data))
		col_cov <- col_cov / (r * length(component_data))
		row_cov <- make_spd(row_cov)
		col_cov <- make_spd(col_cov)

		row_covariances[[component]] <- row_cov
		col_covariances[[component]] <- col_cov
		row_scale <- r / sum(diag(row_covariances[[component]]))
		row_covariances[[component]] <- row_covariances[[component]] * row_scale
		col_covariances[[component]] <- col_covariances[[component]] / row_scale
		row_covariances[[component]] <- make_spd(row_covariances[[component]])
		col_covariances[[component]] <- make_spd(col_covariances[[component]])
	}

	list(
		pi = mixing_proportions,
		M = mean_matrices,
		U = row_covariances,
		V = col_covariances,
		cluster = cluster_assignments
	)
}

#' Compute Mahalanobis Distances for Non-Noise Observations
#'
#' Computes Mahalanobis distance for each observation assigned to a Gaussian
#' component (cluster > 0), using each observation's assigned component parameters.
#'
#' @param fit Fitted model with cluster, M, U, V.
#' @param x_list List of matrices.
#' @return Numeric vector of finite Mahalanobis distances (sorted).
#' @noRd
mv_component_distances <- function(fit, x_list) {
	keep_idx <- which(fit$cluster > 0)
	if (length(keep_idx) < 2) {
		return(numeric(0))
	}

	distances <- vapply(keep_idx, function(i) {
		comp <- fit$cluster[i]
		mv_mahalanobis(
			x = x_list[[i]],
			mean_matrix = fit$M[[comp]],
			row_cov = fit$U[[comp]],
			col_cov = fit$V[[comp]]
		)
	}, numeric(1))

	distances[is.finite(distances)]
}
