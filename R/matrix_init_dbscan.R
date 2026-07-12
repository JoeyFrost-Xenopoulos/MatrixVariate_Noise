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
#' @noRd
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

#' DBSCAN Heuristic Epsilon
#'
#' @param x_matrix Numeric matrix of vectorized observations.
#' @param minPts Minimum neighborhood size.
#'
#' @return Numeric epsilon value.
#' @noRd
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

#' DBSCAN Cluster Assignments
#'
#' @param x_matrix Numeric matrix of vectorized observations.
#' @param eps Neighborhood radius.
#' @param minPts Minimum neighborhood size.
#'
#' @return Integer vector of cluster assignments (0 = noise).
#' @noRd
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
