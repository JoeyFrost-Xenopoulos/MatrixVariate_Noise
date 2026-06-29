#' Score HC Noise Fit with a Matrix KS Test
#'
#' Supports two scoring modes:
#' - `"onesample"` (default): one-sample KS test comparing Mahalanobis
#'   distances against the theoretical $\chi^2(r p)$ CDF.
#' - `"twosample"`: two-sample KS test comparing Mahalanobis distance
#'   distributions between the two largest non-noise mixture components.
#'   Under a correct model, both groups should follow the same
#'   $\chi^2(r p)$ distribution, so the KS statistic should be small.
#'
#' @param fit A fitted noise model.
#' @param x_list List of matrices used for fitting.
#' @param ks_type Character: `"onesample"` or `"twosample"`.
#'
#' @return A list with `statistic`, `p.value`, `n_used`, and (for
#'   `"twosample"`) `n_group1`, `n_group2`, and `ks_type`.
#' @keywords internal
matrix_noise_ks_score <- function(fit, x_list,
                                   ks_type = c("onesample", "twosample")) {
  ks_type <- match.arg(ks_type)

  if (is.null(fit$cluster) || is.null(fit$M) || is.null(fit$U) || is.null(fit$V)) {
    stop("'fit' must contain 'cluster', 'M', 'U', and 'V' components.")
  }
  x_list <- matrix_validate_x_list(x_list)

  distances <- matrix_component_distances(fit, x_list)
  if (length(distances) < 2 || length(unique(distances)) < 2) {
    return(list(statistic = Inf, p.value = NA_real_, n_used = length(distances),
                ks_type = ks_type))
  }

  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])
  df <- r * p

  if (ks_type == "onesample") {
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
      return(list(statistic = Inf, p.value = NA_real_, n_used = length(distances),
                  ks_type = ks_type))
    }

    list(
      statistic = unname(test$statistic),
      p.value = unname(test$p.value),
      n_used = length(distances),
      ks_type = ks_type
    )
  } else {
    keep_idx <- which(fit$cluster > 0)
    clusters <- fit$cluster[keep_idx]

    component_counts <- table(clusters)
    valid_components <- as.integer(names(component_counts)[component_counts >= 2])

    if (length(valid_components) >= 2) {
      sorted_components <- valid_components[order(
        component_counts[as.character(valid_components)], decreasing = TRUE)]
      comp1 <- sorted_components[1]
      comp2 <- sorted_components[2]
      group1 <- distances[clusters == comp1]
      group2 <- distances[clusters == comp2]
    } else {
      set.seed(42)
      half <- floor(length(distances) / 2)
      split_idx <- sample(seq_along(distances), size = half)
      group1 <- distances[split_idx]
      group2 <- distances[-split_idx]
    }

    if (length(group1) < 2 || length(group2) < 2) {
      return(list(statistic = Inf, p.value = NA_real_, n_used = length(distances),
                  ks_type = ks_type))
    }

    test <- tryCatch(
      suppressWarnings(stats::ks.test(group1, group2)),
      error = function(e) {
        warning(sprintf(
          "Two-sample KS test failed: %s",
          conditionMessage(e)
        ), call. = FALSE)
        NULL
      }
    )

    if (is.null(test)) {
      return(list(statistic = Inf, p.value = NA_real_, n_used = length(distances),
                  ks_type = ks_type))
    }

    list(
      statistic = unname(test$statistic),
      p.value = unname(test$p.value),
      n_used = length(distances),
      ks_type = ks_type,
      n_group1 = length(group1),
      n_group2 = length(group2)
    )
  }
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