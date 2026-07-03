#' Score HC Noise Fit with a Matrix KS Test
#'
#' Supports two scoring modes:
#' - `"onesample"` (default): one-sample KS test comparing Mahalanobis
#'   distances against the theoretical $\chi^2(r p)$ CDF.
#'
#' @param fit A fitted noise model.
#' @param x_list List of matrices used for fitting.
#'
#' @return A list with `statistic`, `p.value`, `n_used`.
#' @keywords internal
matrix_noise_ks_score <- function(fit, x_list) {
  
  if (is.null(fit$cluster) || is.null(fit$M) ||
      is.null(fit$U) || is.null(fit$V)) {
    stop("'fit' must contain 'cluster', 'M', 'U', and 'V' components.")
  }
  
  x_list <- matrix_validate_x_list(x_list)
  
  distances <- matrix_component_distances(fit, x_list)
  
  if (length(distances) < 2 || length(unique(distances)) < 2) {
    return(list(
      statistic = Inf,
      p.value = NA_real_,
      n_used = length(distances)
    ))
  }
  
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])
  df <- r * p
  
  test <- tryCatch(
    suppressWarnings(stats::ks.test(distances, "pchisq", df = df)),
    error = function(e) {
      warning(sprintf(
        "One-sample KS test failed (df = %d, n = %d): %s",
        df, length(distances), conditionMessage(e)
      ), call. = FALSE)
      NULL
    }
  )
  
  if (is.null(test)) {
    return(list(
      statistic = Inf,
      p.value = NA_real_,
      n_used = length(distances)
    ))
  }
  
  list(
    statistic = unname(test$statistic),
    p.value = unname(test$p.value),
    n_used = length(distances)
  )
}

#' Generate Dimension-Aware Heuristic Grid for HC Noise
#'
#' Creates a grid of candidate noise_k values based on matrix dimensions.
#' The heuristic centers the grid around 10^(-0.75 * dimension) where
#' dimension = rows * cols.
#'
#' @param x_list List of matrices used for fitting.
#' @param n_points Integer: number of points in the grid.
#' @return Numeric vector of candidate noise_k values.
#' @keywords internal
matrix_noise_hc_heuristic_grid <- function(x_list, n_points = 30) {
  x_list <- matrix_validate_x_list(x_list)
  
  dimension <- nrow(x_list[[1]]) * ncol(x_list[[1]])
  
  if (!is.finite(dimension) || dimension <= 0) {
    # Fallback to default grid
    return(10^seq(-16, -1, length.out = n_points))
  }
  
  # Center at -0.75 * dimension (empirical heuristic)
  center_log10 <- -0.75 * dimension
  
  # Width adapts to dimension: larger dimension needs wider search
  half_width <- max(6, ceiling(dimension / 2))
  
  # Ensure we don't go below machine precision
  lower_log10 <- max(log10(.Machine$double.xmin), center_log10 - half_width)
  upper_log10 <- center_log10 + half_width
  
  # Generate grid on log10 scale
  grid_log10 <- seq(lower_log10, upper_log10, length.out = n_points)
  grid <- 10^grid_log10
  
  # Remove any inf or NaN values
  grid <- grid[is.finite(grid) & grid > 0]
  
  # Ensure we have at least some points
  if (length(grid) < 2) {
    grid <- 10^seq(-16, -1, length.out = n_points)
  }
  
  sort(unique(grid))
}