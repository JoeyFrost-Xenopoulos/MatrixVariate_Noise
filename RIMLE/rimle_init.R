#' Random Initialization for RIMLE
#'
#' Randomly partitions observations into G non-empty groups and computes
#' initial parameters using matrix-normal flip-flop updates.
#'
#' @param x_list List of matrices.
#' @param g Number of components.
#' @return Initial parameter list.
#' @noRd
rimle_random_init <- function(x_list, g) {
	x_list <- rimle_validate_x_list(x_list)
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	cluster_assignments <- sample(seq_len(g), n, replace = TRUE)
	while (length(unique(cluster_assignments)) < g) {
		cluster_assignments <- sample(seq_len(g), n, replace = TRUE)
	}

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

		for (iter in seq_len(20)) {
			col_cov <- make_spd(col_cov)
			row_cov_new <- matrix(0, r, r)
			for (x in component_data) {
				centered <- x - mean_matrices[[component]]
				row_cov_new <- row_cov_new + centered %*% solve(col_cov, t(centered))
			}
			row_cov_new <- row_cov_new / (p * length(component_data))
			row_cov_new <- make_spd(row_cov_new)
			row_scale <- r / sum(diag(row_cov_new))
			row_cov_new <- make_spd(row_cov_new * row_scale)
			col_cov <- col_cov / row_scale
			col_cov <- make_spd(col_cov)
			row_cov <- row_cov_new
		}

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

#' Hennig-Coretto Style Initialization for RIMLE
#'
#' Identifies low-density observations via q-NND and clusters the remaining
#' points using agglomerative hierarchical clustering.
#'
#' @param x_list List of matrices.
#' @param g Number of components.
#' @param pi_max Maximum noise proportion.
#' @param q Nearest neighbor order (default 3).
#' @return Initial parameter list with noise.
#' @noRd
rimle_hennig_coretto_init <- function(x_list, g, pi_max = 0.5, q = 3) {
	x_list <- rimle_validate_x_list(x_list)
	n <- length(x_list)
	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])

	if (!is.numeric(q) || length(q) != 1 || q < 1 || q > n - 1 || q != as.integer(q)) {
		stop(sprintf("'q' must be an integer between 1 and n - 1 (n = %d).", n))
	}
	q <- as.integer(q)

	frob_dists <- matrix(0, n, n)
	for (i in seq_len(n)) {
		for (j in seq_len(n)) {
			frob_dists[i, j] <- norm(x_list[[i]] - x_list[[j]], type = "F")
		}
	}
	diag(frob_dists) <- Inf

	knn_dists <- apply(frob_dists, 1, function(row) sort(row)[q])

	noise_threshold <- quantile(knn_dists, probs = 1 - pi_max, names = FALSE)

	noise_idx <- which(knn_dists > noise_threshold)
	regular_idx <- which(knn_dists <= noise_threshold)

	pi0_init <- max(0.01, length(noise_idx) / n)
	pi0_init <- min(pi0_init, pi_max)

	if (length(regular_idx) < g) {
		warning("Not enough regular observations for HC init; falling back to random init.", call. = FALSE)
		init <- rimle_random_init(x_list, g)
		init$pi0 <- max(0.01, 1 - sum(init$pi))
		return(init)
	}

	regular_data <- x_list[regular_idx]

	if (requireNamespace("mclust", quietly = TRUE)) {
		x_mat <- do.call(rbind, lapply(regular_data, as.vector))
		hc_model <- mclust::hc(modelName = "EEE", data = x_mat)
		cluster_assignments <- mclust::hclass(hc_model, g)
	} else {
		warning("mclust package not available; using base R hclust for initialization.", call. = FALSE)
		x_mat <- do.call(rbind, lapply(regular_data, as.vector))
		d <- dist(x_mat)
		hc <- hclust(d, method = "ward.D2")
		cluster_assignments <- cutree(hc, k = g)
	}

	mixing_proportions <- numeric(g)
	mean_matrices <- vector("list", g)
	row_covariances <- vector("list", g)
	col_covariances <- vector("list", g)

	for (component in seq_len(g)) {
		component_index <- which(cluster_assignments == component)
		if (length(component_index) == 0) {
			component_index <- ((component - 1L) %% length(regular_idx)) + 1L
		}

		component_data <- regular_data[component_index]
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

		for (iter in seq_len(20)) {
			col_cov <- make_spd(col_cov)
			row_cov_new <- matrix(0, r, r)
			for (x in component_data) {
				centered <- x - mean_matrices[[component]]
				row_cov_new <- row_cov_new + centered %*% solve(col_cov, t(centered))
			}
			row_cov_new <- row_cov_new / (p * length(component_data))
			row_cov_new <- make_spd(row_cov_new)
			row_scale <- r / sum(diag(row_cov_new))
			row_cov_new <- make_spd(row_cov_new * row_scale)
			col_cov <- col_cov / row_scale
			col_cov <- make_spd(col_cov)
			row_cov <- row_cov_new
		}

		row_covariances[[component]] <- row_cov
		col_covariances[[component]] <- col_cov
	}

	mixing_proportions <- mixing_proportions / sum(mixing_proportions) * (1 - pi0_init)

	cluster_assignments_full <- integer(n)
	cluster_assignments_full[noise_idx] <- 0L
	cluster_assignments_full[regular_idx] <- cluster_assignments

	list(
		pi0 = pi0_init,
		pi = mixing_proportions,
		M = mean_matrices,
		U = row_covariances,
		V = col_covariances,
		cluster = cluster_assignments_full
	)
}
