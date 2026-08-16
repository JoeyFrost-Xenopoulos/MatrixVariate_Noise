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
mv_mixture_emrefine_init <- function(x_list, g, max_iter = 100) {
  params <- mv_mixture_kmeans_init(x_list, g = g)
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  responsibilities <- matrix(0, n, g)
  colnames(responsibilities) <- paste0("component_", seq_len(g))

  for (iteration in seq_len(max_iter)) {
    log_density <- mv_e_step_log_density(x_list, params, g, n)
    responsibilities <- mv_normalize_responsibilities(log_density)

    new_params <- mv_em_mstep(params, x_list, responsibilities, g, n, r, p,
                              warn_zero = TRUE)
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
