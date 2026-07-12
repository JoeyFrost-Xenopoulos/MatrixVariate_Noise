#' EM-Refine Initialization for Matrix Mixture Models
#'
#' Starts from k-means++ initial parameters and refines them by running
#' a fixed number of EM iterations.
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
#' @param g Integer: number of mixture components
#' @param max_iter Integer: maximum EM iterations (default: 100)
#'
#' @return A list containing initial parameters.
#' @noRd
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
