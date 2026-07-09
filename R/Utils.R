#' Validate Matrix List Input
#'
#' @param x_list List of matrices to validate.
#' @return A list of same-sized matrices.
#' @keywords internal
matrix_validate_x_list <- function(x_list) {
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
#' @keywords internal
matrix_log_sum_exp <- function(values) {
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

# --- Shared EM Utilities ---

#' Dispatch Initialization Scheme
#'
#' Selects and runs the appropriate initialization method.
#'
#' @param x_list Validated list of matrices.
#' @param g Number of components.
#' @param init One of "kmeans" or "emrefine".
#' @param nstart Number of k-means restarts (used only for kmeans init).
#' @return Initial parameter list (pi, M, U, V, cluster).
#' @keywords internal
matrix_init_dispatch <- function(x_list, g, init, nstart = 10) {
	if (init == "emrefine") {
		matrix_mixture_emrefine_init(x_list, g = g)
	} else {
		matrix_mixture_kmeans_init(x_list, g = g, nstart = nstart)
	}
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
#' @keywords internal
matrix_e_step_log_density <- function(x_list, params, g, n) {
	log_density <- matrix(NA_real_, nrow = n, ncol = g)
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
	log_density
}

#' Normalize Log-Densities to Posterior Responsibilities
#'
#' Applies row-wise log-sum-exp normalization.
#'
#' @param log_density Matrix of log-densities (n x K).
#' @return Matrix of posterior responsibilities (n x K), rows sum to 1.
#' @keywords internal
matrix_normalize_responsibilities <- function(log_density) {
	n <- nrow(log_density)
	responsibilities <- matrix(0, nrow = n, ncol = ncol(log_density))
	for (i in seq_len(n)) {
		normalizer <- matrix_log_sum_exp(log_density[i, ])
		responsibilities[i, ] <- exp(log_density[i, ] - normalizer)
	}
	responsibilities
}

#' Compute Weighted Mean Matrix (M-Step)
#'
#' @param x_list List of matrices.
#' @param weights Numeric vector of responsibilities for one component.
#' @param weights_sum Sum of weights (effective sample size).
#' @param r Number of rows.
#' @param p Number of columns.
#' @return Weighted mean matrix (r x p).
#' @keywords internal
matrix_weighted_mean <- function(x_list, weights, weights_sum, r, p) {
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
#' @keywords internal
matrix_update_row_cov <- function(x_list, mean_matrix, v_inv_target, weights,
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
#' @keywords internal
matrix_update_col_cov <- function(x_list, mean_matrix, u_inv_target, weights,
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
#' @keywords internal
matrix_compute_init_params <- function(x_list, g, cluster_assignments, init_method = "Initialization") {
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
#' @keywords internal
matrix_component_distances <- function(fit, x_list) {
	keep_idx <- which(fit$cluster > 0)
	if (length(keep_idx) < 2) {
		return(numeric(0))
	}

	distances <- vapply(keep_idx, function(i) {
		comp <- fit$cluster[i]
		matrix_mahalanobis(
			x = x_list[[i]],
			mean_matrix = fit$M[[comp]],
			row_cov = fit$U[[comp]],
			col_cov = fit$V[[comp]]
		)
	}, numeric(1))

	distances[is.finite(distances)]
}