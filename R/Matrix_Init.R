#' K-Means++ Initialization for Matrix Mixture Models
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
#' @param g Integer: number of mixture components
#' @param nstart Integer: number of independent starts (default: 10)
#'
#' @return A list containing initial parameters.
#' @keywords internal
matrix_mixture_kmeans_init <- function(x_list, g, nstart = 10) {
  x_list <- matrix_validate_x_list(x_list)
  n <- length(x_list)

  if (!is.numeric(nstart) || length(nstart) != 1 || !is.finite(nstart) || nstart < 1) {
    stop("'nstart' must be a positive numeric scalar.")
  }
  nstart <- as.integer(nstart)

  init_basis <- matrix_init_whitening_basis(x_list)
  x_matrix <- matrix_whitened_vectorized_matrices(x_list, init_basis)

  best_fit <- NULL
  best_score <- -Inf

  for (restart in seq_len(nstart)) {
    centers <- matrix_kmeanspp_centers(x_matrix, g, n)
    fit <- tryCatch(
      kmeans(x_matrix, centers = centers, nstart = 1),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      next
    }

    candidate <- matrix_compute_init_params(x_list, g, fit$cluster, init_method = "K-means")
    candidate <- matrix_short_em_burn_in(candidate, x_list, g, max_iter = 3L)
    score <- matrix_initialization_loglik(candidate, x_list, g)

    if (is.finite(score) && score > best_score) {
      best_score <- score
      best_fit <- candidate
    }
  }

  if (is.null(best_fit)) {
    fallback <- kmeans(x_matrix, centers = matrix_kmeanspp_centers(x_matrix, g, n), nstart = 1)
    best_fit <- matrix_compute_init_params(x_list, g, fallback$cluster, init_method = "K-means")
    best_fit <- matrix_short_em_burn_in(best_fit, x_list, g, max_iter = 3L)
  }

  best_fit
}

matrix_init_whitening_basis <- function(x_list) {
  x_list <- matrix_validate_x_list(x_list)
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

matrix_whitened_vectorized_matrices <- function(x_list, init_basis) {
  row_whitener <- solve(chol(init_basis$row_cov))
  col_whitener <- t(solve(chol(init_basis$col_cov)))

  do.call(rbind, lapply(x_list, function(x) {
    centered <- x - init_basis$mean
    as.vector(row_whitener %*% centered %*% col_whitener)
  }))
}

matrix_initialization_loglik <- function(params, x_list, g) {
  n <- length(x_list)
  log_density <- matrix_e_step_log_density(x_list, params, g, n)
  sum(apply(log_density, 1, matrix_log_sum_exp))
}

matrix_short_em_burn_in <- function(params, x_list, g, max_iter = 3L) {
  x_list <- matrix_validate_x_list(x_list)
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  if (!is.numeric(max_iter) || length(max_iter) != 1 || !is.finite(max_iter) || max_iter < 0) {
    stop("'max_iter' must be a non-negative numeric scalar.")
  }
  max_iter <- as.integer(max_iter)

  for (iteration in seq_len(max_iter)) {
    log_density <- matrix_e_step_log_density(x_list, params, g, n)
    responsibilities <- matrix_normalize_responsibilities(log_density)
    component_sizes <- colSums(responsibilities)
    new_params <- params

    for (component in seq_len(g)) {
      if (component_sizes[component] <= 0) {
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

    if (sum(new_params$pi) > 0) {
      new_params$pi <- new_params$pi / sum(new_params$pi)
    }
    params <- new_params
  }

  params$cluster <- max.col(matrix_normalize_responsibilities(
    matrix_e_step_log_density(x_list, params, g, n)
  ), ties.method = "first")

  params
}

matrix_kmeanspp_centers <- function(x_matrix, g, n) {
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

#' DBSCAN-Based Initialization for Matrix Mixture Models
#'
#' Runs DBSCAN on vectorized, scaled matrices and uses the discovered dense
#' regions to seed a k-means refinement. If DBSCAN finds fewer than `g`
#' clusters, the remaining centers are filled with farthest-point seeds.
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
#' @param g Integer: number of mixture components
#' @param eps Numeric: DBSCAN neighborhood radius. If `NULL`, a heuristic is
#'   used based on the scaled k-nearest-neighbor distances.
#' @param minPts Integer: DBSCAN minimum neighborhood size. If `NULL`, a
#'   dimension-aware heuristic is used.
#'
#' @return A list containing initial parameters.
#' @keywords internal
matrix_mixture_dbscan_init <- function(x_list, g, eps = NULL, minPts = NULL) {
  x_list <- matrix_validate_x_list(x_list)
  n <- length(x_list)

  if (n < 2L) {
    return(matrix_mixture_kmeans_init(x_list, g = g))
  }

  x_matrix <- do.call(rbind, lapply(x_list, function(x) as.vector(x)))
  x_scaled <- scale(x_matrix)
  x_scaled[!is.finite(x_scaled)] <- 0

  if (is.null(minPts)) {
    minPts <- max(3L, min(n, ceiling(log2(n + 1L)) + 1L))
  }
  minPts <- as.integer(max(2L, minPts))
  minPts <- min(minPts, n)

  if (is.null(eps)) {
    eps <- matrix_dbscan_heuristic_eps(x_scaled, minPts)
  }

  cluster_assignments <- matrix_dbscan_cluster_assignments(x_scaled, eps, minPts)
  cluster_table <- sort(table(cluster_assignments[cluster_assignments > 0]), decreasing = TRUE)

  if (length(cluster_table) == 0L) {
    warning(
      "DBSCAN initialization found no dense clusters; falling back to k-means++ initialization.",
      call. = FALSE
    )
    return(matrix_mixture_kmeans_init(x_list, g = g))
  }

  selected_cluster_ids <- as.integer(names(cluster_table))[seq_len(min(g, length(cluster_table)))]
  centers <- do.call(rbind, lapply(selected_cluster_ids, function(component) {
    colMeans(x_matrix[cluster_assignments == component, , drop = FALSE])
  }))

  chosen_indices <- unlist(lapply(selected_cluster_ids, function(component) {
    which(cluster_assignments == component)[1L]
  }), use.names = FALSE)
  chosen_indices <- chosen_indices[is.finite(chosen_indices)]

  if (length(chosen_indices) == 0L) {
    chosen_indices <- sample.int(n, 1L)
  }

  while (length(chosen_indices) < g) {
    current_centers <- x_matrix[chosen_indices, , drop = FALSE]
    current_dists <- rep(Inf, n)

    for (center_row in seq_len(nrow(current_centers))) {
      center <- current_centers[center_row, , drop = FALSE]
      center_dists <- rowSums((x_matrix - matrix(center, nrow = n, ncol = ncol(x_matrix), byrow = TRUE))^2)
      current_dists <- pmin(current_dists, center_dists)
    }

    current_dists[chosen_indices] <- -Inf
    next_index <- which.max(current_dists)

    if (!is.finite(current_dists[next_index])) {
      remaining <- setdiff(seq_len(n), chosen_indices)
      next_index <- sample(remaining, 1L)
    }

    chosen_indices <- c(chosen_indices, next_index)
  }

  if (nrow(centers) < g) {
    centers <- rbind(centers, x_matrix[chosen_indices[(nrow(centers) + 1L):g], , drop = FALSE])
  }

  km <- tryCatch(
    kmeans(x_matrix, centers = centers[seq_len(g), , drop = FALSE], nstart = 1),
    error = function(e) {
      warning(
        sprintf("DBSCAN initialization fallback to k-means++ after k-means failed: %s", conditionMessage(e)),
        call. = FALSE
      )
      matrix_mixture_kmeans_init(x_list, g = g)$cluster
    }
  )

  if (is.integer(km)) {
    return(matrix_compute_init_params(x_list, g, km, init_method = "DBSCAN"))
  }

  matrix_compute_init_params(x_list, g, km$cluster, init_method = "DBSCAN")
}

matrix_dbscan_heuristic_eps <- function(x_matrix, minPts) {
  n <- nrow(x_matrix)
  dist_matrix <- as.matrix(stats::dist(x_matrix))
  sorted_distances <- t(apply(dist_matrix, 1, sort))
  kth_index <- min(minPts, ncol(sorted_distances))
  kth_distances <- sorted_distances[, kth_index]
  finite_distances <- kth_distances[is.finite(kth_distances) & kth_distances > 0]

  if (length(finite_distances) == 0L) {
    upper_distances <- dist_matrix[upper.tri(dist_matrix)]
    finite_distances <- upper_distances[is.finite(upper_distances) & upper_distances > 0]
  }

  if (length(finite_distances) == 0L) {
    return(1)
  }

  stats::median(finite_distances)
}

matrix_dbscan_cluster_assignments <- function(x_matrix, eps, minPts) {
  n <- nrow(x_matrix)
  dist_matrix <- as.matrix(stats::dist(x_matrix))
  visited <- rep(FALSE, n)
  cluster_assignments <- integer(n)
  cluster_id <- 0L

  for (point in seq_len(n)) {
    if (visited[point]) {
      next
    }

    visited[point] <- TRUE
    neighbors <- which(dist_matrix[point, ] <= eps)

    if (length(neighbors) < minPts) {
      next
    }

    cluster_id <- cluster_id + 1L
    cluster_assignments[point] <- cluster_id
    seed_points <- setdiff(neighbors, point)

    while (length(seed_points) > 0L) {
      current_point <- seed_points[1L]
      seed_points <- seed_points[-1L]

      if (!visited[current_point]) {
        visited[current_point] <- TRUE
        current_neighbors <- which(dist_matrix[current_point, ] <= eps)

        if (length(current_neighbors) >= minPts) {
          seed_points <- union(seed_points, setdiff(current_neighbors, current_point))
        }
      }

      if (cluster_assignments[current_point] == 0L) {
        cluster_assignments[current_point] <- cluster_id
      }
    }
  }

  cluster_assignments
}

matrix_mixture_emrefine_init <- function(x_list, g, max_iter = 100) {
  params <- matrix_mixture_kmeans_init(x_list, g = g)
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  responsibilities <- matrix(0, n, g)
  colnames(responsibilities) <- paste0("component_", seq_len(g))

  for (iteration in seq_len(max_iter)) {
    log_density <- matrix_e_step_log_density(x_list, params, g, n)
    responsibilities <- matrix_normalize_responsibilities(log_density)

    component_sizes <- colSums(responsibilities)
    new_params <- params

    for (component in seq_len(g)) {
      if (component_sizes[component] <= 0) {
        warning(sprintf(
          "EM-refine initialization: component %d has zero effective membership at iteration %d; skipping update.",
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
  }

  list(
    pi = params$pi,
    M = params$M,
    U = params$U,
    V = params$V,
    cluster = max.col(responsibilities, ties.method = "first")
  )
}

