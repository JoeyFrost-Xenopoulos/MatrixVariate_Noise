#' Compute Whitening Basis for Matrix Mixture Initialization
#'
#' Computes pooled row and column covariance estimates from the data,
#' scaled to enforce the identifiability constraint tr(U) = r.
#'
#' @param x_list A list of numeric matrices, each of dimension r x p
#'
#' @return A list with elements `mean`, `row_cov`, and `col_cov`.
#' @noRd
mv_init_whitening_basis <- function(x_list) {
	x_list <- mv_validate_x_list(x_list)
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	mean_matrix <- Reduce(`+`, x_list) / n
	row_cov <- matrix(0, r, r)
	col_cov <- matrix(0, p, p)

	for (x in x_list) {
		centered <- x - mean_matrix
		row_cov <- row_cov + centered %*% t(centered)
		col_cov <- col_cov + t(centered) %*% centered
	}

	row_cov <- make_spd(row_cov / (p * n))
	col_cov <- make_spd(col_cov / (r * n))

	row_scale <- r / sum(diag(row_cov))
	row_cov <- make_spd(row_cov * row_scale)
	col_cov <- make_spd(col_cov / row_scale)

	list(mean = mean_matrix, row_cov = row_cov, col_cov = col_cov)
}

#' Whitened Vectorized Matrices
#'
#' Centers and whitens each matrix in x_list using the provided whitening basis,
#' then vectorizes the result row-wise.
#'
#' @param x_list A list of numeric matrices.
#' @param init_basis Whitening basis from `mv_init_whitening_basis`.
#'
#' @return A numeric matrix (n x (r*p)) of whitened, vectorized matrices.
#' @noRd
mv_whitened_vectorized_matrices <- function(x_list, init_basis) {
	row_whitener <- solve(chol(init_basis$row_cov))
	col_whitener <- t(solve(chol(init_basis$col_cov)))

	do.call(rbind, lapply(x_list, function(x) {
		centered <- x - init_basis$mean
		as.vector(row_whitener %*% centered %*% col_whitener)
	}))
}

#' K-Means++ Center Seeding
#'
#' @param x_matrix Numeric matrix of vectorized observations (n x d).
#' @param g Number of centers to select.
#' @param n Number of observations.
#'
#' @return A numeric matrix (g x d) of selected centers.
#' @noRd
mv_kmeanspp_centers <- function(x_matrix, g, n) {
	centers_idx <- integer(g)
	centers_idx[1] <- sample.int(n, 1)
	min_dists <- rep(Inf, n)

	for (component in 2:g) {
		last_center <- x_matrix[centers_idx[component - 1], , drop = FALSE]
		current_dists <- rowSums((x_matrix - matrix(
			last_center,
			nrow = n,
			ncol = ncol(x_matrix),
			byrow = TRUE
		))^2)
		min_dists <- pmin(min_dists, current_dists)

		if (sum(min_dists) <= 0 || !is.finite(sum(min_dists))) {
			centers_idx[component] <- sample.int(n, 1)
		} else {
			probs <- min_dists / sum(min_dists)
			centers_idx[component] <- sample.int(n, 1, prob = probs)
		}
	}

	x_matrix[centers_idx, , drop = FALSE]
}

#' Short EM Burn-In for Initialization
#'
#' Performs a small number of EM iterations to refine initial parameters
#' before the main fitting loop.
#'
#' @param params Initial parameter list.
#' @param x_list List of matrices.
#' @param g Number of components.
#' @param max_iter Maximum iterations (default: 3).
#'
#' @return Refined parameter list with cluster assignments.
#' @noRd
mv_short_em_burn_in <- function(params, x_list, g, max_iter = 3L) {
	x_list <- mv_validate_x_list(x_list)
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	if (!is.numeric(max_iter) || length(max_iter) != 1 || !is.finite(max_iter) || max_iter < 0) {
		stop("'max_iter' must be a non-negative numeric scalar.")
	}
	max_iter <- as.integer(max_iter)

	for (iteration in seq_len(max_iter)) {
		log_density <- mv_e_step_log_density(x_list, params, g, n)
		responsibilities <- mv_normalize_responsibilities(log_density)
		new_params <- mv_em_mstep(params, x_list, responsibilities, g, n, r, p,
		                          warn_zero = FALSE)
		if (sum(new_params$pi) > 0) {
			new_params$pi <- new_params$pi / sum(new_params$pi)
		}
		params <- new_params
	}

	params$cluster <- max.col(mv_normalize_responsibilities(
		mv_e_step_log_density(x_list, params, g, n)
	), ties.method = "first")

	params
}

#' Initialization Log-Likelihood
#'
#' Computes the observed-data log-likelihood for a set of initial parameters.
#'
#' @param params Initial parameter list.
#' @param x_list List of matrices.
#' @param g Number of components.
#'
#' @return Numeric log-likelihood.
#' @noRd
mv_initialization_loglik <- function(params, x_list, g) {
	n <- length(x_list)
	log_density <- mv_e_step_log_density(x_list, params, g, n)
	mv_loglik(log_density)
}

#' K-Means++ Initialization for Matrix Mixture Models
#'
#' @param x_list A list of numeric matrices, each of dimension r x p
#' @param g Integer: number of mixture components
#' @param nstart Integer: number of independent starts (default: 10)
#' @param use_parallel Logical: if TRUE, evaluate the nstart restarts in
#'   parallel (future / multisession workers). FALSE is the sequential
#'   fallback used for debugging.
#' @param n_cores Integer: number of parallel workers (NULL = auto).
#'
#' @return A list containing initial parameters.
#' @noRd
mv_mixture_kmeans_init <- function(x_list, g, nstart = 10,
                                    use_parallel = FALSE, n_cores = NULL) {
	x_list <- mv_validate_x_list(x_list)
	n <- length(x_list)

	if (!is.numeric(nstart) || length(nstart) != 1 || !is.finite(nstart) || nstart < 1) {
		stop("'nstart' must be a positive numeric scalar.")
	}
	nstart <- as.integer(nstart)

	init_basis <- mv_init_whitening_basis(x_list)
	x_matrix <- mv_whitened_vectorized_matrices(x_list, init_basis)

	run_one_restart <- function(restart) {
		centers <- mv_kmeanspp_centers(x_matrix, g, n)
		fit <- tryCatch(
			kmeans(x_matrix, centers = centers, nstart = 1),
			error = function(e) NULL
		)

		if (is.null(fit)) {
			return(list(fit = NULL, score = -Inf))
		}

		candidate <- mv_compute_init_params(x_list, g, fit$cluster, init_method = "K-means")
		candidate <- mv_short_em_burn_in(candidate, x_list, g, max_iter = 3L)
		score <- mv_initialization_loglik(candidate, x_list, g)

		list(fit = candidate, score = score)
	}

	config <- mv_parallel_config(
		use_parallel = use_parallel,
		n_cores = n_cores,
		requested = "restart",
		n_tasks = nstart
	)

	if (config$active) {
		results <- mv_future_lapply(seq_len(nstart), run_one_restart, config)
	} else {
		results <- lapply(seq_len(nstart), run_one_restart)
	}

	best_fit <- NULL
	best_score <- -Inf
	for (res in results) {
		if (is.finite(res$score) && res$score > best_score) {
			best_score <- res$score
			best_fit <- res$fit
		}
	}

	if (is.null(best_fit)) {
		fallback <- kmeans(x_matrix, centers = mv_kmeanspp_centers(x_matrix, g, n), nstart = 1)
		best_fit <- mv_compute_init_params(x_list, g, fallback$cluster, init_method = "K-means")
		best_fit <- mv_short_em_burn_in(best_fit, x_list, g, max_iter = 3L)
	}

	best_fit
}





#' Dispatch Initialization Scheme
#'
#' Selects and runs the appropriate initialization method.
#'
#' @param x_list Validated list of matrices.
#' @param g Number of components.
#' @param nstart Number of k-means restarts.
#' @param use_parallel Logical: enable parallel nstart restarts for kmeans init.
#' @param n_cores Integer: number of parallel workers (NULL = auto).
#' @return Initial parameter list (pi, M, U, V, cluster).
#' @noRd
mv_init_dispatch <- function(x_list, g, nstart = 10,
                             use_parallel = FALSE, n_cores = NULL) {
	mv_mixture_kmeans_init(x_list, g = g, nstart = nstart,
	                       use_parallel = use_parallel, n_cores = n_cores)
}

#' Resolve Parallel Configuration
#'
#' Validates and normalizes the user-facing parallel knobs and returns a list
#' describing whether/when parallel dispatch should happen.
#'
#' @param use_parallel Logical master switch.
#' @param n_cores Integer number of workers, or NULL for auto.
#' @param parallel_strategy Character: "auto", "grid", or "restart".
#' @param requested Character: which parallel layer is asking
#'   ("grid" or "restart"). Ignored unless parallel_strategy != "auto".
#' @param seed Integer seed or NULL.
#' @param n_tasks Integer: number of independent tasks the active parallel
#'   layer will dispatch.
#' @return A list with active (logical), n_cores (integer), strategy
#'   (resolved character), and seed.
#' @noRd
mv_parallel_config <- function(use_parallel, n_cores, parallel_strategy = "auto",
                               requested = NULL, seed = NULL, n_tasks = NULL) {
	stopifnot(
		is.logical(use_parallel), length(use_parallel) == 1, !is.na(use_parallel)
	)

	if (!is.null(n_cores)) {
		if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) {
			stop("'n_cores' must be a positive integer (or NULL for automatic).")
		}
		n_cores <- as.integer(n_cores)
	}

	stopifnot(is.null(seed) || (is.numeric(seed) && length(seed) == 1))
	if (!is.null(seed)) seed <- as.integer(seed)

	if (!is.null(n_tasks)) {
		stopifnot(is.numeric(n_tasks), length(n_tasks) == 1, n_tasks >= 1)
		n_tasks <- as.integer(n_tasks)
	}

	strategy <- match.arg(parallel_strategy, c("auto", "grid", "restart"))

	if (use_parallel && !requireNamespace("future", quietly = TRUE)) {
		stop(
			"'use_parallel = TRUE' requires the 'future' package. ",
			"Install it with install.packages(\"future\") or set use_parallel = FALSE."
		)
	}

	active <- FALSE
	if (use_parallel) {
		if (is.null(requested)) {
			active <- TRUE
		} else if (strategy == "auto") {
			active <- TRUE
		} else {
			active <- (strategy == requested)
		}
	}

	if (active && is.null(n_cores)) {
		n_cores <- mv_optimal_n_cores(n_tasks = n_tasks)
	}

	list(active = active, n_cores = n_cores %||% 1L, strategy = strategy, seed = seed)
}

#' Choose an Optimal Number of Parallel Workers
#'
#' Picks a sensible default worker count for a clustering job.
#'
#' @param n_tasks Integer: number of independent tasks the parallel layer will
#'   dispatch. NULL means "unknown yet".
#' @param max_workers Integer: hard upper bound on workers (default 8).
#' @return A positive integer number of workers to use.
#' @noRd
mv_optimal_n_cores <- function(n_tasks = NULL, max_workers = 8L) {
	if (!requireNamespace("future", quietly = TRUE)) {
		return(1L)
	}
	cores <- tryCatch(
		as.integer(future::availableCores()),
		error = function(e) 1L
	)
	workers <- tryCatch(
		length(future::availableWorkers()),
		error = function(e) cores
	)
	n <- min(cores, workers, max_workers)
	if (!is.null(n_tasks)) {
		n <- min(n, max(1L, as.integer(n_tasks)))
	}
	max(1L, n)
}

#' Run a Future-Lapply Over a List
#'
#' @param X A vector or list to iterate over.
#' @param FUN Function applied to each element of X.
#' @param config List from mv_parallel_config().
#' @param ... Extra arguments forwarded to FUN.
#' @return A list of results.
#' @noRd
mv_future_lapply <- function(X, FUN, config, ...) {
	stopifnot(is.list(config), is.logical(config$active))

	run_task <- function(x, idx, ...) {
		mv_parallel_worker_setup()
		rng_save <- mv_rng_state_save()
		on.exit(mv_rng_state_restore(rng_save), add = TRUE)
		if (!is.null(config$seed)) {
			RNGkind("L'Ecuyer-CMRG")
			set.seed(mv_task_seed(config$seed, idx))
		}
		FUN(x, ...)
	}

	if (!config$active) {
		return(lapply(seq_along(X), function(i) run_task(X[[i]], i, ...)))
	}

	if (!requireNamespace("future", quietly = TRUE) ||
	    !requireNamespace("future.apply", quietly = TRUE)) {
		stop(
			"'use_parallel = TRUE' requires the 'future' and 'future.apply' packages. ",
			"Install them with install.packages(c(\"future\", \"future.apply\")) ",
			"or set use_parallel = FALSE."
		)
	}

	n_workers <- config$n_cores
	n_tasks <- length(X)
	n_workers <- min(n_workers, max(1L, n_tasks))

	prev_plan <- future::plan(
		future::multisession,
		workers = n_workers,
		.init = FALSE
	)
	on.exit({
		ok <- tryCatch(
			future::plan(prev_plan, .init = FALSE),
			error = function(e) FALSE
		)
		if (isFALSE(ok)) {
			warning("Failed to restore previous future::plan(); default plan may have changed.")
		}
	}, add = TRUE)

	future.apply::future_lapply(
		seq_along(X),
		function(i, ...) run_task(X[[i]], i, ...),
		...,
		future.seed = FALSE,
		future.scheduling = FALSE
	)
}