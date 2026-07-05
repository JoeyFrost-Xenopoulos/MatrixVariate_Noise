#' Enforce Positive Definiteness on a Matrix
#'
#' Converts a matrix to symmetric positive definite form using iterative jittering
#' of the diagonal. This is necessary for numerical stability when computing
#' Cholesky decompositions and matrix inverses.
#'
#' @param mat A numeric matrix to be made positive definite
#' @param jitter Initial jitter amount added to diagonal (default: 1e-8)
#' @param max_tries Maximum number of jittering attempts (default: 8)
#'
#' @return A symmetric positive definite matrix
#'
#' @details
#' The function:
#' 1. Symmetrizes the matrix by averaging with its transpose
#' 2. Attempts Cholesky decomposition with increasing jitter amounts
#' 3. Returns the first successful candidate or errors if max_tries exceeded
#'
#' @keywords internal
make_spd <- function(mat, jitter = 1e-8, max_tries = 8) {
	if (!is.matrix(mat) || !is.numeric(mat)) {
		stop("'mat' must be a numeric matrix.")
	}
	if (nrow(mat) != ncol(mat)) {
		stop("'mat' must be a square matrix.")
	}
	mat <- (mat + t(mat)) / 2
	for (k in 0:max_tries) {
		j <- jitter * (10^k)
		candidate <- mat + diag(j, nrow(mat))
		ok <- tryCatch({
			chol(candidate)
			TRUE
		}, error = function(e) FALSE)
		if (ok) return(candidate)
	}
	stop("Could not make covariance matrix positive definite.")
}

#' Compute Matrix-Variate Mahalanobis Distance
#'
#' Calculates the Mahalanobis distance between a matrix and a mean matrix
#' under the matrix-variate normal distribution with specified row and column
#' covariance structures.
#'
#' @param x A numeric matrix (r × p): the observation
#' @param mean_matrix A numeric matrix (r × p): the component mean
#' @param row_cov A numeric matrix (r × r): row covariance matrix U
#' @param col_cov A numeric matrix (p × p): column covariance matrix V
#'
#' @return Numeric scalar representing the Mahalanobis distance
#'
#' @details
#'
#' This metric extends the multivariate Mahalanobis distance to account for
#' the matrix structure. The computation uses Cholesky decomposition and
#' forward/backsolve for numerical stability.
#'
#' @export
matrix_mahalanobis <- function(x, mean_matrix, row_cov, col_cov) {
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
	# U^{-1} and V^{-1}
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
#' @param x A numeric matrix (r × p): the observation
#' @param mean_matrix A numeric matrix (r × p): the component mean matrix M
#' @param row_cov A numeric matrix (r × r): row covariance matrix U
#' @param col_cov A numeric matrix (p × p): column covariance matrix V
#'
#' @return Numeric scalar representing the log-density value
#'
#' @details
#'
#' Computation uses Cholesky decomposition for numerical stability and to
#' avoid explicit matrix inversion.
#'
#' @keywords internal
matrix_variate_log_density <- function(x, mean_matrix, row_cov, col_cov) {
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
	# Cholesky decomposition
	row_cov <- make_spd(row_cov)
	col_cov <- make_spd(col_cov)
	row_chol <- chol(row_cov)
	col_chol <- chol(col_cov)

	# |U| and |V| for the denominator
	row_logdet <- 2 * sum(log(diag(row_chol)))
	col_logdet <- 2 * sum(log(diag(col_chol)))

	# tr(V^{-1} * (X - M)^T * U^{-1} * (X - M))
	centered <- x - mean_matrix
	row_inv_centered <- backsolve(row_chol, forwardsolve(t(row_chol), centered))
	col_inv <- chol2inv(col_chol)
	trace_form <- sum(row_inv_centered * (centered %*% col_inv))

	r <- nrow(x)
	p <- ncol(x)

	# Returns log density value
	-0.5 * (r * p * log(2 * pi) + p * row_logdet + r * col_logdet + trace_form)
}

#' Fit Matrix-Variate Gaussian Mixture Model via EM Algorithm
#'
#' Estimates parameters of a matrix-variate Gaussian mixture model (MGMM)
#' using the Expectation-Maximization algorithm. Performs clustering of
#' matrix-valued observations while accounting for row and column dependencies.
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
#' @param g Integer: number of mixture components
#' @param max_iter Integer: maximum EM iterations (default: 100)
#' @param tol Numeric: convergence tolerance for log-likelihood (default: 1e-6)
#' @param nstart Integer: number of k-means restarts for initialization (default: 10). Ignored unless `init = "kmeans"`.
#' @param init Character: initialization scheme. `"kmeans"` (default), `"random"`, or `"ecme"`.
#' @param verbose Logical: print iteration progress (default: FALSE)
#'
#' @return A list containing:
#' - `pi`: numeric vector of length g with final mixing proportions.
#' - `M`: list of g final component mean matrices.
#' - `U`: list of g final row covariance matrices.
#' - `V`: list of g final column covariance matrices.
#' - `z`: numeric matrix (n × g) of posterior responsibilities.
#' - `cluster`: integer vector of length n with hard cluster assignments.
#' - `logLik`: numeric vector with the log-likelihood trace across iterations.
#' - `iterations`: number of EM iterations performed.
#' - `converged`: logical indicating whether the algorithm converged within `max_iter`.
#'
#' @details
#' The EM algorithm alternates between:
#' **E-step:** Compute posterior responsibilities (soft cluster assignments)
#' **M-step:** Update parameters based on responsibilities
#' 
#' @examples
#' \dontrun{
#' set.seed(123)
#' mean_1 <- matrix(c(1.5, 1.2, 1.0, 1.3, 1.1, 1.4, 1.2, 1.0), nrow=2)
#' mean_2 <- matrix(c(-1.4, -1.0, -1.2, -1.3, -1.1, -1.5, -1.0, -1.2), nrow=2)
#'
#' simulate_matrix_group <- function(n, mean_matrix, row_sd=0.35, col_sd=0.35) {
#'   r <- nrow(mean_matrix); p <- ncol(mean_matrix)
#'   row_cov <- diag(row_sd, r); col_cov <- diag(col_sd, p)
#'   lapply(seq_len(n), function(i) {
#'     noise <- matrix(rnorm(r*p), r, p)
#'     mean_matrix + row_cov %*% noise %*% col_cov
#'   })
#' }
#'
#' x_list <- c(
#'   simulate_matrix_group(15, mean_1),
#'   simulate_matrix_group(15, mean_2)
#' )
#'
#' fit <- matrix_variate_mixture_fit(x_list, g=2, max_iter=50, verbose=TRUE)
#' fit$cluster
#' fit$pi
#' }
#'
#' @export
matrix_variate_mixture_fit <- function(x_list, g, max_iter = 100, tol = 1e-06,
																			 nstart = 10, init = c("kmeans", "random", "ecme", "kmeans++"),
																			 verbose = FALSE) {
	init <- match.arg(init)
	x_list <- matrix_validate_x_list(x_list)
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	if (!is.numeric(g) || length(g) != 1 || g < 1) {
		stop("'g' must be a positive integer specifying the number of mixture components.")
	}
	g <- as.integer(g)

	if (n < g) {
		stop(sprintf(
			"Number of observations (%d) must be at least as large as the number of components (%d).",
			n, g
		))
	}

	params <- matrix_init_dispatch(x_list, g, init, nstart)
	loglik_trace <- numeric(0)
	responsibilities <- matrix(0, n, g)

	# EM loop
	for (iteration in seq_len(max_iter)) {
		# E-step
		log_density <- matrix_e_step_log_density(x_list, params, g, n)
		responsibilities <- matrix_normalize_responsibilities(log_density)

		# Observed data log-likelihood
		current_loglik <- sum(apply(log_density, 1, matrix_log_sum_exp))
		loglik_trace <- c(loglik_trace, current_loglik)

		if (iteration > 1 && abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
			break
		}

		# M-step
		component_sizes <- colSums(responsibilities)
		new_params <- params

		for (component in seq_len(g)) {
			if (component_sizes[component] <= 0) {
				warning(sprintf(
					"Component %d has zero effective membership at iteration %d; skipping update.",
					component, iteration
				), call. = FALSE)
				next
			}

			weights <- responsibilities[, component]
			weights_sum <- component_sizes[component]

			mean_matrix <- matrix_weighted_mean(x_list, weights, weights_sum, r, p)
			row_cov <- matrix_update_row_cov(x_list, mean_matrix, params$V[[component]],
			                                 weights, weights_sum, r, p)
			col_cov <- matrix_update_col_cov(x_list, mean_matrix, row_cov,
			                                 weights, weights_sum, r, p)

			new_params$pi[component] <- weights_sum / n
			new_params$M[[component]] <- mean_matrix
			new_params$U[[component]] <- row_cov
			new_params$V[[component]] <- col_cov
		}

		new_params$pi <- new_params$pi / sum(new_params$pi)
		params <- new_params

		if (verbose) {
			message(sprintf("Iteration %d: log-likelihood = %.4f", iteration, current_loglik))
		}
	}

	cluster_membership <- max.col(responsibilities, ties.method = "first")

	list(
		pi = params$pi,
		M = params$M,
		U = params$U,
		V = params$V,
		z = responsibilities,
		cluster = cluster_membership,
		logLik = loglik_trace,
		iterations = length(loglik_trace),
		converged = length(loglik_trace) < max_iter
	)
}
