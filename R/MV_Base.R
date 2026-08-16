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
#' @export
#' @noRd
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
#' @param init Character: initialization scheme. `"kmeans"` (default), `"emrefine"`, or `"dbscan"`.
#' @param verbose Logical: print iteration progress (default: FALSE)
#' @param use_parallel Logical: if `TRUE`, run the k-means `nstart` restarts in
#'   parallel via the **future** package. `FALSE` (default) runs sequentially.
#' @param n_cores Integer: number of parallel workers (NULL = auto).
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
#' fit <- mv_mixture_fit(x_list, g=2, max_iter=50, verbose=TRUE)
#' fit$cluster
#' fit$pi
#' }
#'
#' @export
mv_mixture_fit <- function(x_list, g, max_iter = 100, tol = 1e-06,
																			 nstart = 10, init = c("kmeans", "emrefine", "dbscan"),
																			 verbose = FALSE, use_parallel = FALSE,
																			 n_cores = NULL) {
	init <- match.arg(init)
	x_list <- mv_validate_x_list(x_list)
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

	params <- mv_init_dispatch(x_list, g, init, nstart,
	                          use_parallel = use_parallel, n_cores = n_cores)
	loglik_trace <- numeric(0)
	responsibilities <- matrix(0, n, g)
	
	if (verbose) {
	  message(sprintf("Fitting: %d components, max_iter=%d", g, max_iter))
	}

	# EM loop
	for (iteration in seq_len(max_iter)) {
		# E-step
		log_density <- mv_e_step_log_density(x_list, params, g, n)
		responsibilities <- mv_normalize_responsibilities(log_density)

		# Observed data log-likelihood
		current_loglik <- mv_loglik(log_density)
		loglik_trace <- c(loglik_trace, current_loglik)

		if (iteration > 1 && abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
			break
		}

		# M-step
		new_params <- mv_em_mstep(params, x_list, responsibilities, g, n, r, p,
		                          warn_zero = TRUE)
		new_params$pi <- new_params$pi / sum(new_params$pi)
		params <- new_params

		if (verbose) {
			message(sprintf("Iteration %d: log-likelihood = %.4f", iteration, current_loglik))
		}
	}
	
	if (verbose) {
	  message(sprintf("Converged in %d iterations (max_iter=%d).",
	                  length(loglik_trace), max_iter))
	}

	cluster_membership <- max.col(responsibilities, ties.method = "first")

	structure(list(
		pi = params$pi,
		M = params$M,
		U = params$U,
		V = params$V,
		z = responsibilities,
		cluster = cluster_membership,
		logLik = loglik_trace,
		iterations = length(loglik_trace),
		converged = length(loglik_trace) < max_iter
	), class = "mv_mixture_fit")
}

#' Score HC Noise Fit with a Matrix KS Test
#'
#' Supports two scoring modes:
#' - `"onesample"` (default): one-sample KS test comparing Mahalanobis
#'   distances against the theoretical chi-square CDF.
#'
#' @param fit A fitted noise model.
#' @param x_list List of matrices used for fitting.
#'
#' @return A list with `statistic`, `p.value`, `n_used`.
#' @noRd
mv_noise_ks_score <- function(fit, x_list) {
  
  if (is.null(fit$cluster) || is.null(fit$M) ||
      is.null(fit$U) || is.null(fit$V)) {
    stop("'fit' must contain 'cluster', 'M', 'U', and 'V' components.")
  }
  
  x_list <- mv_validate_x_list(x_list)
  
  distances <- mv_component_distances(fit, x_list)
  
  if (length(distances) < 2 || length(unique(distances)) < 2) {
    return(list(
      statistic = Inf,
      p.value = NA_real_,
      n_used = length(distances)
    ))
  }
  
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])
  df <- r * p
  
  test <- tryCatch(
    suppressWarnings(stats::ks.test(distances, "pchisq", df = df)),
    error = function(e) {
      warning(sprintf(
        "One-sample KS test failed (df = %d, n = %d): %s",
        df, length(distances), conditionMessage(e)
      ), call. = FALSE)
      NULL
    }
  )
  
  if (is.null(test)) {
    return(list(
      statistic = Inf,
      p.value = NA_real_,
      n_used = length(distances)
    ))
  }
  
  list(
    statistic = unname(test$statistic),
    p.value = unname(test$p.value),
    n_used = length(distances)
  )
}