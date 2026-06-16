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
#' @keywords internal
matrix_mahalanobis <- function(x, mean_matrix, row_cov, col_cov) {
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
																			 nstart = 10, init = c("kmeans", "random", "ecme"),
																			 verbose = FALSE) {
	init <- match.arg(init)
	x_list <- matrix_validate_x_list(x_list)
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	if (init == "random") {
		params <- matrix_mixture_random_init(x_list, g = g)
	} else if (init == "ecme") {
		params <- matrix_mixture_ecme_init(x_list, g = g)
	} else {
		params <- matrix_mixture_kmeans_init(x_list, g = g, nstart = nstart)
	}
	loglik_trace <- numeric(0)
	responsibilities <- matrix(0, n, g)

	# EM loop
	for (iteration in seq_len(max_iter)) {
		# E-step: Compute responsibilities
		log_density <- matrix(NA_real_, nrow = n, ncol = g)

		for (component in seq_len(g)) {
			for (i in seq_len(n)) {
				# compute P(X_i | component)
				log_density[i, component] <- log(params$pi[component]) +
					matrix_variate_log_density(
						x = x_list[[i]],
						mean_matrix = params$M[[component]],
						row_cov = params$U[[component]],
						col_cov = params$V[[component]]
					)
			}
		}

		# Normalize log-densities to get responsibilities
		for (i in seq_len(n)) {
			row_log_densities <- log_density[i, ]
			normalizer <- matrix_log_sum_exp(row_log_densities)  # numerically stable log-sum-exp
			responsibilities[i, ] <- exp(row_log_densities - normalizer)  # z_hat_ig
		}

		# Compute observed data log-likelihood for convergence check
		current_loglik <- sum(apply(log_density, 1, matrix_log_sum_exp))
		loglik_trace <- c(loglik_trace, current_loglik)

		# Check convergence
		if (iteration > 1 && abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
			break
		}

		# M-step: Update parameters
		component_sizes <- colSums(responsibilities)  # sum of z_hat_ig over all observations i
		new_params <- params

		# Update each component's parameters
		for (component in seq_len(g)) {
			if (component_sizes[component] <= 0) {
				next
			}

			weights <- responsibilities[, component]
			weights_sum <- component_sizes[component]  # effective sample size for this component
			v_for_row <- make_spd(params$V[[component]])

			# Update mean matrix (M-step: M_hat_g)
			mean_matrix <- matrix(0, r, p)
			for (i in seq_len(n)) {
				mean_matrix <- mean_matrix + weights[i] * x_list[[i]]
			}
			mean_matrix <- mean_matrix / weights_sum

			# Update row covariance U_hat_g
			row_cov <- matrix(0, r, r)
			for (i in seq_len(n)) {
				centered <- x_list[[i]] - mean_matrix
				row_cov <- row_cov + weights[i] * (centered %*% solve(v_for_row, t(centered)))
			}
			row_cov <- row_cov / (p * weights_sum)
			row_cov <- make_spd(row_cov)

			# Enforce tr(U) = r
			row_scale <- r / sum(diag(row_cov))
			row_cov <- row_cov * row_scale
			row_cov <- make_spd(row_cov)

			# Update column covariance V_hat_g
			col_cov <- matrix(0, p, p)
			for (i in seq_len(n)) {
				centered <- x_list[[i]] - mean_matrix
				col_cov <- col_cov + weights[i] * (t(centered) %*% solve(row_cov, centered))
			}
			col_cov <- col_cov / (r * weights_sum)
			col_cov <- make_spd(col_cov)

			# Store updated parameters
			new_params$pi[component] <- weights_sum / n  # mixing proportion
			new_params$M[[component]] <- mean_matrix
			new_params$U[[component]] <- row_cov
			new_params$V[[component]] <- col_cov
		}

		# Normalize mixing proportions to sum to 1
		new_params$pi <- new_params$pi / sum(new_params$pi)

		params <- new_params

		if (verbose) {
			message(sprintf("Iteration %d: log-likelihood = %.4f", iteration, current_loglik))
		}
	}

	# Assign each observation to its most likely component
	cluster_membership <- max.col(responsibilities, ties.method = "first")

	# Return fitted model with all parameters and diagnostics
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
