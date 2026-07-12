#' Compute Whitening Basis for Matrix Mixture Initialization
#'
#' Computes pooled row and column covariance estimates from the data,
#' scaled to enforce the identifiability constraint tr(U) = r.
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
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
#' @return A numeric matrix (n × (r*p)) of whitened, vectorized matrices.
#' @noRd
mv_whitened_vectorized_matrices <- function(x_list, init_basis) {
  row_whitener <- solve(chol(init_basis$row_cov))
  col_whitener <- t(solve(chol(init_basis$col_cov)))

  do.call(rbind, lapply(x_list, function(x) {
    centered <- x - init_basis$mean
    as.vector(row_whitener %*% centered %*% col_whitener)
  }))
}
