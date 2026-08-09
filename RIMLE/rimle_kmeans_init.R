#' Whitening Basis for RIMLE Initialization
#'
#' Computes pooled row and column covariance estimates from the data,
#' scaled to enforce the identifiability constraint tr(U) = r.
#'
#' @param x_list A list of numeric matrices, each of dimension r x p
#' @return A list with elements `mean`, `row_cov`, and `col_cov`.
#' @noRd
rimle_init_whitening_basis <- function(x_list) {
	x_list <- rimle_validate_x_list(x_list)
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
#' @param init_basis Whitening basis from `rimle_init_whitening_basis`.
#' @return A numeric matrix (n x (r*p)) of whitened, vectorized matrices.
#' @noRd
rimle_whitened_vectorized_matrices <- function(x_list, init_basis) {
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
#' @return A numeric matrix (g x d) of selected centers.
#' @noRd
rimle_kmeanspp_centers <- function(x_matrix, g, n) {
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

#' Compute Initial Parameters from Cluster Assignments
#'
#' @param x_list Validated list of matrices.
#' @param g Number of components.
#' @param cluster_assignments Integer vector of cluster labels (1..g).
#' @return A list with pi, M, U, V.
#' @noRd
rimle_compute_init_params <- function(x_list, g, cluster_assignments) {
	x_list <- rimle_validate_x_list(x_list)
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

		row_scale <- r / sum(diag(row_cov))
		row_cov <- make_spd(row_cov * row_scale)
		col_cov <- make_spd(col_cov / row_scale)

		row_covariances[[component]] <- row_cov
		col_covariances[[component]] <- col_cov
	}

	list(
		pi = mixing_proportions,
		M = mean_matrices,
		U = row_covariances,
		V = col_covariances,
		cluster = cluster_assignments
	)
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
#' @return Refined parameter list with cluster assignments.
#' @noRd
rimle_short_em_burn_in <- function(params, x_list, g, max_iter = 3L) {
	x_list <- rimle_validate_x_list(x_list)
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	if (!is.numeric(max_iter) || length(max_iter) != 1 || !is.finite(max_iter) || max_iter < 0) {
		stop("'max_iter' must be a non-negative numeric scalar.")
	}
	max_iter <- as.integer(max_iter)

	log_density <- matrix(0, n, g)
	for (comp in seq_len(g)) {
		for (i in seq_len(n)) {
			log_density[i, comp] <- log(params$pi[comp]) +
				rimle_mv_log_density(x_list[[i]], params$M[[comp]],
						     params$U[[comp]], params$V[[comp]])
		}
	}
	responsibilities <- matrix(0, n, g)
	for (i in seq_len(n)) {
		normalizer <- rimle_log_sum_exp(log_density[i, ])
		responsibilities[i, ] <- exp(log_density[i, ] - normalizer)
	}

	for (iteration in seq_len(max_iter)) {
		component_sizes <- colSums(responsibilities)
		new_params <- params

		for (component in seq_len(g)) {
			if (component_sizes[component] <= 0) next

			weights <- responsibilities[, component]
			weights_sum <- component_sizes[component]

			mean_matrix <- matrix(0, r, p)
			for (i in seq_len(n)) {
				mean_matrix <- mean_matrix + weights[i] * x_list[[i]]
			}
			mean_matrix <- mean_matrix / weights_sum

			row_cov <- matrix(0, r, r)
			col_cov <- matrix(0, p, p)
			for (i in seq_len(n)) {
				centered <- x_list[[i]] - mean_matrix
				row_cov <- row_cov + weights[i] * (centered %*% t(centered))
				col_cov <- col_cov + weights[i] * (t(centered) %*% centered)
			}
			row_cov <- row_cov / (p * weights_sum)
			col_cov <- col_cov / (r * weights_sum)
			row_cov <- make_spd(row_cov)
			col_cov <- make_spd(col_cov)

			row_scale <- r / sum(diag(row_cov))
			row_cov <- make_spd(row_cov * row_scale)
			col_cov <- make_spd(col_cov / row_scale)

			new_params$pi[component] <- weights_sum / n
			new_params$M[[component]] <- mean_matrix
			new_params$U[[component]] <- row_cov
			new_params$V[[component]] <- col_cov
		}

		params <- new_params

		log_density <- matrix(0, n, g)
		for (comp in seq_len(g)) {
			for (i in seq_len(n)) {
				log_density[i, comp] <- log(params$pi[comp]) +
					rimle_mv_log_density(x_list[[i]], params$M[[comp]],
							     params$U[[comp]], params$V[[comp]])
			}
		}
		for (i in seq_len(n)) {
			normalizer <- rimle_log_sum_exp(log_density[i, ])
			responsibilities[i, ] <- exp(log_density[i, ] - normalizer)
		}
	}

	params$cluster <- max.col(responsibilities, ties.method = "first")
	params
}

#' Initialization Log-Likelihood
#'
#' Computes the observed-data log-likelihood for a set of initial parameters.
#'
#' @param params Initial parameter list.
#' @param x_list List of matrices.
#' @param g Number of components.
#' @return Numeric log-likelihood.
#' @noRd
rimle_initialization_loglik <- function(params, x_list, g) {
	n <- length(x_list)
	log_density <- matrix(0, n, g)
	for (comp in seq_len(g)) {
		for (i in seq_len(n)) {
			log_density[i, comp] <- log(params$pi[comp]) +
				rimle_mv_log_density(x_list[[i]], params$M[[comp]],
						     params$U[[comp]], params$V[[comp]])
		}
	}
	total <- 0
	for (i in seq_len(n)) {
		total <- total + rimle_log_sum_exp(log_density[i, ])
	}
	total
}

#' K-Means++ Initialization for RIMLE
#'
#' @param x_list A list of numeric matrices, each of dimension r x p
#' @param g Integer: number of mixture components
#' @param nstart Integer: number of independent starts (default: 10)
#' @return A list containing initial parameters.
#' @noRd
rimle_kmeans_init <- function(x_list, g, nstart = 10) {
	x_list <- rimle_validate_x_list(x_list)
	n <- length(x_list)

	if (!is.numeric(nstart) || length(nstart) != 1 || !is.finite(nstart) || nstart < 1) {
		stop("'nstart' must be a positive numeric scalar.")
	}
	nstart <- as.integer(nstart)

	init_basis <- rimle_init_whitening_basis(x_list)
	x_matrix <- rimle_whitened_vectorized_matrices(x_list, init_basis)

	run_one_restart <- function(restart) {
		centers <- rimle_kmeanspp_centers(x_matrix, g, n)
		fit <- tryCatch(
			kmeans(x_matrix, centers = centers, nstart = 1),
			error = function(e) NULL
		)

		if (is.null(fit)) {
			return(list(fit = NULL, score = -Inf))
		}

		candidate <- rimle_compute_init_params(x_list, g, fit$cluster)
		candidate <- rimle_short_em_burn_in(candidate, x_list, g, max_iter = 3L)
		score <- rimle_initialization_loglik(candidate, x_list, g)

		list(fit = candidate, score = score)
	}

	results <- lapply(seq_len(nstart), run_one_restart)

	best_fit <- NULL
	best_score <- -Inf
	for (res in results) {
		if (is.finite(res$score) && res$score > best_score) {
			best_score <- res$score
			best_fit <- res$fit
		}
	}

	if (is.null(best_fit)) {
		centers <- rimle_kmeanspp_centers(x_matrix, g, n)
		fit <- kmeans(x_matrix, centers = centers, nstart = 1)
		best_fit <- rimle_compute_init_params(x_list, g, fit$cluster)
		best_fit <- rimle_short_em_burn_in(best_fit, x_list, g, max_iter = 3L)
	}

	best_fit
}
