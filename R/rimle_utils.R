#' Validate RIMLE Input
#'
#' @param x_list List of matrices to validate.
#' @return A validated list of same-sized matrices.
#' @noRd
rimle_validate_x_list <- function(x_list) {
	if (!is.list(x_list) || length(x_list) == 0) {
		stop("x_list must be a non-empty list of matrices.")
	}
	if (!is.matrix(x_list[[1]])) {
		stop("First element of x_list is not a matrix.")
	}
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])
	if (r == 0 || p == 0) {
		stop("Matrices must have at least one row and one column.")
	}
	for (idx in seq_along(x_list)) {
		x <- x_list[[idx]]
		if (!is.matrix(x)) {
			stop(sprintf("Element %d of x_list is not a matrix.", idx))
		}
		if (nrow(x) != r || ncol(x) != p) {
			stop(sprintf("Element %d has dimensions %d x %d, expected %d x %d.",
				     idx, nrow(x), ncol(x), r, p))
		}
	}
	x_list
}

#' Stable Log-Sum-Exp
#'
#' @param values Numeric vector.
#' @return Numeric scalar.
#' @noRd
rimle_log_sum_exp <- function(values) {
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

#' Matrix-Variate Normal Log-Density
#'
#' Evaluates the log-density of a matrix observation under the matrix-variate
#' normal distribution with mean matrix M, row covariance U, and column
#' covariance V.
#'
#' @param x A numeric matrix (r x p).
#' @param M A numeric matrix (r x p).
#' @param U A numeric matrix (r x r): row covariance.
#' @param V A numeric matrix (p x p): column covariance.
#' @return Numeric scalar log-density.
#' @noRd
rimle_mv_log_density <- function(x, M, U, V) {
	if (!is.matrix(x) || !is.numeric(x)) stop("'x' must be a numeric matrix.")
	if (!is.matrix(M) || !is.numeric(M)) stop("'M' must be a numeric matrix.")
	if (!identical(dim(x), dim(M))) stop("'x' and 'M' must have the same dimensions.")
	if (!is.matrix(U) || nrow(U) != ncol(U) || nrow(U) != nrow(x)) {
		stop("'U' must be square with dimension matching nrow(x).")
	}
	if (!is.matrix(V) || nrow(V) != ncol(V) || ncol(V) != ncol(x)) {
		stop("'V' must be square with dimension matching ncol(x).")
	}

	U <- make_spd(U)
	V <- make_spd(V)
	U_chol <- chol(U)
	V_chol <- chol(V)

	U_logdet <- 2 * sum(log(diag(U_chol)))
	V_logdet <- 2 * sum(log(diag(V_chol)))

	centered <- x - M
	U_inv_centered <- backsolve(U_chol, forwardsolve(t(U_chol), centered))
	V_inv <- chol2inv(V_chol)
	trace_form <- sum(U_inv_centered * (centered %*% V_inv))

	r <- nrow(x)
	p <- ncol(x)

	-0.5 * (r * p * log(2 * pi) + p * U_logdet + r * V_logdet + trace_form)
}

#' Matrix-Variate Normal Density
#'
#' @param x Observation matrix.
#' @param M Mean matrix.
#' @param U Row covariance.
#' @param V Column covariance.
#' @return Density scalar.
#' @noRd
rimle_mv_density <- function(x, M, U, V) {
	exp(rimle_mv_log_density(x, M, U, V))
}

#' Shrinkage Operator
#'
#' Implements l_gamma(a, m) = min(max(m, a), gamma * m).
#'
#' @param a Eigenvalue.
#' @param m Shrinkage target.
#' @param gamma Constraint bound.
#' @return Shrunk eigenvalue.
#' @noRd
rimle_shrinkage <- function(a, m, gamma) {
	pmax(pmin(pmax(m, a), gamma * m), 1e-10)
}

#' Compute m_* for CM1 Eigenratio Constraint
#'
#' Solves the one-dimensional convex problem for the shrinkage parameter.
#'
#' @param eig_values_list List of eigenvalue vectors for the matrix being updated.
#' @param Zg Effective sample sizes per component.
#' @param gamma Constraint bound.
#' @return Optimal m_*.
#' @noRd
rimle_compute_m_star <- function(eig_values_list, Zg, gamma) {
	g <- length(eig_values_list)

	objective <- function(m) {
		if (m <= 0) return(Inf)
		total <- 0
		for (comp in seq_len(g)) {
			sh <- rimle_shrinkage(eig_values_list[[comp]], m, gamma)
			if (any(sh <= 0)) return(Inf)
			total <- total + Zg[comp] * sum(log(sh) + eig_values_list[[comp]] / sh)
		}
		total
	}

	all_eigs <- unlist(eig_values_list)
	if (length(all_eigs) == 0 || all(!is.finite(all_eigs))) {
		return(1)
	}

	upper <- max(all_eigs, na.rm = TRUE) * 10
	if (!is.finite(upper) || upper <= 0) upper <- 1

	res <- optimize(objective, interval = c(1e-15, upper), maximum = FALSE, tol = 1e-8)
	res$minimum
}

#' Compute omega_* for CM2 Noise Proportion Constraint
#'
#' Finds the root of the constrained mixing proportion equation.
#'
#' @param x_list List of matrices.
#' @param params Current parameters.
#' @param g Number of Gaussian components.
#' @param n Number of observations.
#' @param k Noise constant.
#' @param Z0 Noise effective sample size.
#' @param pi_max Maximum noise proportion.
#' @return Optimal omega_*.
#' @noRd
rimle_compute_omega_star <- function(x_list, params, g, n, k, Z0, pi_max) {
	if (Z0 >= n - 1e-8) return(min(pi_max, 0.99))

	f <- function(omega) {
		if (omega <= 0 || omega >= 1) return(NA_real_)
		total <- 0
		for (i in seq_len(n)) {
			log_terms <- c(log(omega * k))
			for (comp in seq_len(g)) {
				gauss_dens <- rimle_mv_density(x_list[[i]], params$M[[comp]],
							       params$U[[comp]], params$V[[comp]])
				if (gauss_dens > 0) {
					log_terms <- c(log_terms,
						       log((1 - omega) / (n - Z0)) +
						       log(params$Zg[comp]) +
						       log(gauss_dens))
				}
			}
			log_denom <- rimle_log_sum_exp(log_terms)
			if (is.finite(log_denom) && log_denom > -700) {
				total <- total + exp(log(omega * k) - log_denom)
			}
		}
		total - n * pi_max
	}

	tryCatch({
		res <- uniroot(f, interval = c(1e-8, min(pi_max, 1 - 1e-8)),
			       tol = 1e-8, extendInt = "no")
		min(res$root, pi_max)
	}, error = function(e) {
		omegas <- seq(1e-8, pi_max, length.out = 100)
		f_vals <- sapply(omegas, f)
		if (!any(is.finite(f_vals))) return(min(pi_max, 0.5))
		omegas[which.min(abs(f_vals))]
	})
}

#' Matrix-Variate Mahalanobis Distance
#'
#' Computes the Mahalanobis distance between a matrix observation and a mean
#' matrix under the matrix-variate normal with specified row and column
#' covariances. Uses Cholesky decomposition for numerical stability.
#'
#' @param x A numeric matrix (r x p): the observation.
#' @param mean_matrix A numeric matrix (r x p): the component mean.
#' @param row_cov A numeric matrix (r x r): row covariance.
#' @param col_cov A numeric matrix (p x p): column covariance.
#' @return Numeric scalar: squared Mahalanobis distance.
#' @noRd
rimle_mv_mahalanobis <- function(x, mean_matrix, row_cov, col_cov) {
  if (!is.matrix(x) || !is.numeric(x)) stop("'x' must be a numeric matrix.")
  if (!is.matrix(mean_matrix) || !is.numeric(mean_matrix)) stop("'mean_matrix' must be a numeric matrix.")
  if (!identical(dim(x), dim(mean_matrix))) stop("'x' and 'mean_matrix' must have the same dimensions.")
  if (!is.matrix(row_cov) || nrow(row_cov) != ncol(row_cov) || nrow(row_cov) != nrow(x)) {
    stop("'row_cov' must be square with dimension matching nrow(x).")
  }
  if (!is.matrix(col_cov) || nrow(col_cov) != ncol(col_cov) || ncol(col_cov) != ncol(x)) {
    stop("'col_cov' must be square with dimension matching ncol(x).")
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

#' Compute Q1 for Convergence Check
#'
#' @param x_list List of matrices.
#' @param params Current parameters.
#' @param g Number of Gaussian components.
#' @param n Number of observations.
#' @return Q1 value.
#' @noRd
rimle_compute_q1 <- function(x_list, params, g, n) {
	total <- 0
	for (comp in seq_len(g)) {
		for (i in seq_len(n)) {
			total <- total + params$z[i, comp] *
				rimle_mv_log_density(x_list[[i]], params$M[[comp]],
						     params$U[[comp]], params$V[[comp]])
		}
	}
	total
}
