#' Check and Fix Eigenratio Constraint at Initialization
#'
#' Repeatedly shrinks eigenvalues until the global eigenratio constraint
#' max_g(λ_max(U_g)*λ_max(V_g)) / min_g(λ_min(U_g)*λ_min(V_g)) <= gamma
#' is satisfied.
#'
#' @param U_list List of row covariance matrices.
#' @param V_list List of column covariance matrices.
#' @param gamma Eigenratio bound.
#' @param max_attempts Maximum shrink attempts.
#' @return List with corrected U_list and V_list.
#' @noRd
rimle_check_and_fix_eigenratio_init <- function(U_list, V_list, gamma,
                                                 max_attempts = 100) {
	g <- length(U_list)
	if (g == 0) return(list(U = U_list, V = V_list))

	for (attempt in seq_len(max_attempts)) {
		all_u_eigs <- unlist(lapply(U_list, function(U) {
			eigen(U, symmetric = TRUE, only.values = TRUE)$values
		}))
		all_v_eigs <- unlist(lapply(V_list, function(V) {
			eigen(V, symmetric = TRUE, only.values = TRUE)$values
		}))

		global_ratio <- max(all_u_eigs * rep(max(all_v_eigs), length(all_u_eigs))) /
			min(all_u_eigs * rep(min(all_v_eigs), length(all_u_eigs)))

		if (global_ratio <= gamma) break

		U_list <- lapply(U_list, function(U) {
			eig <- eigen(U, symmetric = TRUE)
			eig$vectors %*% diag(0.95 * eig$values) %*% t(eig$vectors)
		})
	}

	list(U = U_list, V = V_list)
}

#' Random Initialization for RIMLE
#'
#' Randomly partitions observations into G non-empty groups and computes
#' initial parameters using matrix-normal flip-flop updates.
#'
#' @importFrom mclust hc hclass
#' @param x_list List of matrices.
#' @param g Number of components.
#' @param pi_max Maximum noise proportion.
#' @param gamma Eigenratio bound for constraint enforcement.
#' @return Initial parameter list.
#' @noRd
rimle_random_init <- function(x_list, g, pi_max = 0.5, gamma = 1000) {
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
			V_inv <- solve(col_cov)

			row_cov_new <- matrix(0, r, r)
			for (x in component_data) {
				centered <- x - mean_matrices[[component]]
				row_cov_new <- row_cov_new + centered %*% V_inv %*% t(centered)
			}
			row_cov_new <- row_cov_new / (p * length(component_data))
			row_cov_new <- make_spd(row_cov_new)

			U_inv <- solve(row_cov_new)

			col_cov_new <- matrix(0, p, p)
			for (x in component_data) {
				centered <- x - mean_matrices[[component]]
				col_cov_new <- col_cov_new + t(centered) %*% U_inv %*% centered
			}
			col_cov_new <- col_cov_new / (r * length(component_data))
			col_cov_new <- make_spd(col_cov_new)

			row_scale <- sum(diag(row_cov_new)) / r

			row_cov_new <- make_spd(row_cov_new / row_scale)
			col_cov_new <- make_spd(col_cov_new * row_scale)

			row_cov <- row_cov_new
			col_cov <- col_cov_new
		}

		row_covariances[[component]] <- row_cov
		col_covariances[[component]] <- col_cov
	}

	pi0 <- min(max(0.01, 1 - sum(mixing_proportions)), pi_max)
	remaining <- 1 - pi0
	if (sum(mixing_proportions) > 0) {
		mixing_proportions <- mixing_proportions / sum(mixing_proportions) * remaining
	} else {
		mixing_proportions <- rep(remaining / g, g)
	}

	fixed <- rimle_check_and_fix_eigenratio_init(row_covariances, col_covariances, gamma)

	list(
		pi0 = pi0,
		pi = mixing_proportions,
		M = mean_matrices,
		U = fixed$U,
		V = fixed$V,
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
rimle_hennig_coretto_init <- function(x_list, g, pi_max = 0.5, q = 3, gamma = 1000) {
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
		init <- rimle_random_init(x_list, g, pi_max = pi_max, gamma = gamma)
		return(init)
	}

	regular_data <- x_list[regular_idx]

	if (requireNamespace("mclust", quietly = TRUE)) {
		x_mat <- do.call(rbind, lapply(regular_data, as.vector))
		hc_model <- tryCatch({
			if (!isNamespaceLoaded("mclust")) {
				suppressPackageStartupMessages(library(mclust, character.only = TRUE, quietly = TRUE))
			}
			mclust::hc(modelName = "EEE", data = x_mat)
		}, error = function(e) NULL)
		if (!is.null(hc_model)) {
			cluster_assignments <- mclust::hclass(hc_model, g)
		} else {
			warning("mclust::hc initialization failed; using base R hclust.", call. = FALSE)
			d <- dist(x_mat)
			hc <- hclust(d, method = "ward.D2")
			cluster_assignments <- cutree(hc, k = g)
		}
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
			V_inv <- solve(col_cov)

			row_cov_new <- matrix(0, r, r)
			for (x in component_data) {
				centered <- x - mean_matrices[[component]]
				row_cov_new <- row_cov_new + centered %*% V_inv %*% t(centered)
			}
			row_cov_new <- row_cov_new / (p * length(component_data))
			row_cov_new <- make_spd(row_cov_new)

			U_inv <- solve(row_cov_new)

			col_cov_new <- matrix(0, p, p)
			for (x in component_data) {
				centered <- x - mean_matrices[[component]]
				col_cov_new <- col_cov_new + t(centered) %*% U_inv %*% centered
			}
			col_cov_new <- col_cov_new / (r * length(component_data))
			col_cov_new <- make_spd(col_cov_new)

			row_scale <- sum(diag(row_cov_new)) / r

			row_cov_new <- make_spd(row_cov_new / row_scale)
			col_cov_new <- make_spd(col_cov_new * row_scale)

			row_cov <- row_cov_new
			col_cov <- col_cov_new
		}

		row_covariances[[component]] <- row_cov
		col_covariances[[component]] <- col_cov
	}

	mixing_proportions <- mixing_proportions / sum(mixing_proportions) * (1 - pi0_init)

	cluster_assignments_full <- integer(n)
	cluster_assignments_full[noise_idx] <- 0L
	cluster_assignments_full[regular_idx] <- cluster_assignments

	fixed <- rimle_check_and_fix_eigenratio_init(row_covariances, col_covariances, gamma)

	list(
		pi0 = pi0_init,
		pi = mixing_proportions,
		M = mean_matrices,
		U = fixed$U,
		V = fixed$V,
		cluster = cluster_assignments_full
	)
}

#' K-Means++ Initialization for RIMLE
#'
#' Performs k-means++ clustering in a whitened vectorized space to obtain
#' initial cluster assignments, then derives RIMLE parameters (mean matrices,
#' row/column covariances, mixing proportions) from those assignments. Far-field
#' observations are assigned to the noise component (cluster 0).
#'
#' The matrix-to-vector transformation applies whitening: each matrix is
#' centered by the global mean, then decorrelated and scaled by the pooled
#' row and column covariance Cholesky factors before being vectorized.
#'
#' @param x_list List of matrices.
#' @param g Number of Gaussian components.
#' @param pi_max Maximum noise proportion.
#' @param nstart Number of k-means restarts.
#' @param q Nearest neighbor order for noise threshold.
#' @return Initial parameter list with noise.
#' @noRd
rimle_kmeans_init <- function(x_list, g, pi_max = 0.5, nstart = 10,
                              use_parallel = FALSE, n_cores = NULL, q = 3,
                              gamma = 1000) {
  x_list <- rimle_validate_x_list(x_list)
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  q <- as.integer(q)
  if (q < 1 || q > n - 1) {
    q <- max(1L, min(as.integer(n - 1), as.integer(3)))
  }

  init_basis <- mv_init_whitening_basis(x_list)
  init_basis$row_whitener <- solve(chol(init_basis$row_cov))
  init_basis$col_whitener <- t(solve(chol(init_basis$col_cov)))
  x_matrix <- mv_whitened_vectorized_matrices(x_list, init_basis)

  run_one_restart <- function(restart) {
    centers <- mv_kmeanspp_centers(x_matrix, g, n)
    fit <- tryCatch(
      kmeans(x_matrix, centers = centers, nstart = 1),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(list(cluster = integer(0), noise_idx = integer(0), pi0_init = 0.01,
                  score = Inf))
    }

    candidate <- mv_compute_init_params(x_list, g, fit$cluster, init_method = "K-means")

    whitened_means <- do.call(rbind, lapply(seq_len(g), function(comp) {
      as.vector(init_basis$row_whitener %*% (candidate$M[[comp]] - init_basis$mean) %*% t(init_basis$col_whitener))
    }))

    mahal_dists <- numeric(n)
    for (i in seq_len(n)) {
      min_dist <- Inf
      for (comp in seq_len(g)) {
        d <- sum((x_matrix[i, ] - whitened_means[comp, ])^2)
        if (d < min_dist) min_dist <- d
      }
      mahal_dists[i] <- min_dist
    }

    knn_dists <- numeric(n)
    for (i in seq_len(n)) {
      others <- mahal_dists[-i]
      knn_dists[i] <- sort(others)[min(q, length(others))]
    }

    noise_threshold <- quantile(knn_dists, probs = 1 - pi_max, names = FALSE)
    noise_idx <- which(knn_dists > noise_threshold)
    pi0_init <- max(0.01, length(noise_idx) / n)
    pi0_init <- min(pi0_init, pi_max)

    list(cluster = fit$cluster, noise_idx = noise_idx, pi0_init = pi0_init,
         score = -mean(knn_dists))
  }

  config <- mv_parallel_config(
    use_parallel = use_parallel,
    n_cores = n_cores,
    parallel_strategy = "restart",
    requested = "restart",
    n_tasks = nstart
  )

  if (config$active) {
    results <- mv_future_lapply(seq_len(nstart), run_one_restart, config)
  } else {
    results <- lapply(seq_len(nstart), run_one_restart)
  }

  best_result <- NULL
  best_score <- Inf
  for (res in results) {
    if (length(res$cluster) > 0 && res$score < best_score) {
      best_score <- res$score
      best_result <- res
    }
  }

  if (is.null(best_result)) {
    fallback <- kmeans(x_matrix, centers = mv_kmeanspp_centers(x_matrix, g, n), nstart = 1)
    cluster_assignments <- fallback$cluster
    noise_idx <- integer(0)
    pi0_init <- 0.01
  } else {
    cluster_assignments <- best_result$cluster
    noise_idx <- best_result$noise_idx
    pi0_init <- best_result$pi0_init
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
      V_inv <- solve(col_cov)

      row_cov_new <- matrix(0, r, r)
      for (x in component_data) {
        centered <- x - mean_matrices[[component]]
        row_cov_new <- row_cov_new + centered %*% V_inv %*% t(centered)
      }
      row_cov_new <- row_cov_new / (p * length(component_data))
      row_cov_new <- make_spd(row_cov_new)

      U_inv <- solve(row_cov_new)

      col_cov_new <- matrix(0, p, p)
      for (x in component_data) {
        centered <- x - mean_matrices[[component]]
        col_cov_new <- col_cov_new + t(centered) %*% U_inv %*% centered
      }
      col_cov_new <- col_cov_new / (r * length(component_data))
      col_cov_new <- make_spd(col_cov_new)

      row_scale <- sum(diag(row_cov_new)) / r

      row_cov_new <- make_spd(row_cov_new / row_scale)
      col_cov_new <- make_spd(col_cov_new * row_scale)

      row_cov <- row_cov_new
      col_cov <- col_cov_new
    }

    row_covariances[[component]] <- row_cov
    col_covariances[[component]] <- col_cov
  }

  mixing_proportions <- mixing_proportions / sum(mixing_proportions) * (1 - pi0_init)

  cluster_assignments_full <- integer(n)
  cluster_assignments_full[noise_idx] <- 0L
  cluster_assignments_full[setdiff(seq_len(n), noise_idx)] <- cluster_assignments[setdiff(seq_len(n), noise_idx)]

  fixed <- rimle_check_and_fix_eigenratio_init(row_covariances, col_covariances, gamma)

  list(
    pi0 = pi0_init,
    pi = mixing_proportions,
    M = mean_matrices,
    U = fixed$U,
    V = fixed$V,
    cluster = cluster_assignments_full
  )
}
