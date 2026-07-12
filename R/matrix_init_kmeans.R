#' K-Means++ Initialization for Matrix Mixture Models
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
#' @param g Integer: number of mixture components
#' @param nstart Integer: number of independent starts (default: 10)
#'
#' @return A list containing initial parameters.
#' @noRd
mv_mixture_kmeans_init <- function(x_list, g, nstart = 10) {
  x_list <- mv_validate_x_list(x_list)
  n <- length(x_list)

  if (!is.numeric(nstart) || length(nstart) != 1 || !is.finite(nstart) || nstart < 1) {
    stop("'nstart' must be a positive numeric scalar.")
  }
  nstart <- as.integer(nstart)

  init_basis <- mv_init_whitening_basis(x_list)
  x_matrix <- mv_whitened_vectorized_matrices(x_list, init_basis)

  best_fit <- NULL
  best_score <- -Inf

  for (restart in seq_len(nstart)) {
    centers <- mv_kmeanspp_centers(x_matrix, g, n)
    fit <- tryCatch(
      kmeans(x_matrix, centers = centers, nstart = 1),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      next
    }

    candidate <- mv_compute_init_params(x_list, g, fit$cluster, init_method = "K-means")
    candidate <- mv_short_em_burn_in(candidate, x_list, g, max_iter = 3L)
    score <- mv_initialization_loglik(candidate, x_list, g)

    if (is.finite(score) && score > best_score) {
      best_score <- score
      best_fit <- candidate
    }
  }

  if (is.null(best_fit)) {
    fallback <- kmeans(x_matrix, centers = mv_kmeanspp_centers(x_matrix, g, n), nstart = 1)
    best_fit <- mv_compute_init_params(x_list, g, fallback$cluster, init_method = "K-means")
    best_fit <- mv_short_em_burn_in(best_fit, x_list, g, max_iter = 3L)
  }

  best_fit
}

#' K-Means++ Center Seeding
#'
#' @param x_matrix Numeric matrix of vectorized observations (n × d).
#' @param g Number of centers to select.
#' @param n Number of observations.
#'
#' @return A numeric matrix (g × d) of selected centers.
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
