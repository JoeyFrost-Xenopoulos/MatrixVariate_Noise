#' K-Means Initialization for Matrix Mixture Models
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
#' @param g Integer: number of mixture components
#' @param nstart Integer: number of k-means restarts (default: 10)
#'
#' @return A list containing initial parameters.
#' @keywords internal
matrix_mixture_kmeans_init <- function(x_list, g, nstart = 10) {
  x_list <- matrix_validate_x_list(x_list)

  # vectorize and run kmeans for init
  x_matrix <- do.call(rbind, lapply(x_list, function(x) as.vector(x)))
  km <- kmeans(x_matrix, centers = g, nstart = nstart)

  matrix_compute_init_params(x_list, g, km$cluster, init_method = "K-means")
}

#' K-Means++ Initialization for Matrix Mixture Models
#'
#' Implements the k-means++ seeding algorithm (Arthur & Vassilvitskii, 2007)
#' on vectorized matrix observations, then runs standard k-means from those
#' seeds. The D^2 weighting produces better-spread initial centers than
#' uniform random seeding.
#'
#' @param x_list A list of numeric matrices, each of dimension r x p
#' @param g Integer: number of mixture components
#' @param nstart Integer: number of k-means restarts from the pp-seeded
#'   centers (default: 10)
#'
#' @return A list containing initial parameters.
#' @keywords internal
matrix_mixture_kmeanspp_init <- function(x_list, g, nstart = 10) {
  x_list <- matrix_validate_x_list(x_list)
  n <- length(x_list)

  x_matrix <- do.call(rbind, lapply(x_list, function(x) as.vector(x)))

  # k-means++ center seeding (D^2 weighting)
  centers_idx <- integer(g)
  centers_idx[1] <- sample.int(n, 1)

  for (k in 2:g) {
    # Squared distance from each point to its nearest chosen center
    dists <- apply(x_matrix[centers_idx[1:(k - 1)], , drop = FALSE], 1, function(c) {
      rowSums((x_matrix - matrix(c, nrow = n, ncol = ncol(x_matrix), byrow = TRUE))^2)
    })
    if (is.matrix(dists)) {
      min_dists <- apply(dists, 1, min)
    } else {
      min_dists <- dists
    }
    # Sample proportional to D^2
    probs <- min_dists / sum(min_dists)
    centers_idx[k] <- sample.int(n, 1, prob = probs)
  }

  centers <- x_matrix[centers_idx, , drop = FALSE]

  # Run k-means from the pp-seeded centers
  km <- kmeans(x_matrix, centers = centers, nstart = nstart)

  matrix_compute_init_params(x_list, g, km$cluster, init_method = "K-means++")
}

matrix_mixture_random_init <- function(x_list, g) {
  x_list <- matrix_validate_x_list(x_list)
  n <- length(x_list)

  z <- sample(seq_len(g), size = n, replace = TRUE)

  matrix_compute_init_params(x_list, g, z, init_method = "Random")
}

matrix_mixture_ecme_init <- function(x_list, g, max_iter = 5) {
  params <- matrix_mixture_random_init(x_list, g = g)
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
          "ECME initialization: component %d has zero effective membership at iteration %d; skipping update.",
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

