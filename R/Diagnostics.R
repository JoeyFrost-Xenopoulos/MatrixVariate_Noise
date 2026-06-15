#' Matrix Diagnostics
#'
#' @export
matrix_noise_ecdf_vs_cdf_plot <- function(fit, x_list,
                                          main = "ECDF vs Theoretical CDF (Chi-square)",
                                          xlab = "Mahalanobis distance",
                                          ylab = "CDF",
                                          plot_type = "base") {
  
  plot_type <- match.arg(plot_type)
  x_list <- matrix_validate_x_list(x_list)
  
  keep_idx <- which(fit$cluster > 0)
  if (length(keep_idx) < 2) {
    stop("Not enough non-noise points to compute ECDF.")
  }
  
  # Recompute distances exactly as in KS scoring
  distances <- vapply(keep_idx, function(i) {
    comp <- fit$cluster[i]
    
    matrix_mahalanobis(
      x = x_list[[i]],
      mean_matrix = fit$M[[comp]],
      row_cov = fit$U[[comp]],
      col_cov = fit$V[[comp]]
    )
  }, numeric(1))
  
  distances <- distances[is.finite(distances)]
  distances <- sort(distances)
  
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])
  df <- r * p
  
  ecdf_fn <- stats::ecdf(distances)
  
  x_grid <- unique(stats::quantile(distances, probs = seq(0, 1, length.out = 200)))
  ecdf_vals <- ecdf_fn(x_grid)
  cdf_vals <- stats::pchisq(x_grid, df = df)
  
  if (plot_type == "base") {
    plot(x_grid, ecdf_vals, type = "l", lwd = 2,
         col = "black", ylim = c(0, 1),
         xlab = xlab, ylab = ylab, main = main)
    lines(x_grid, cdf_vals, col = "red", lwd = 2, lty = 2)
    legend("bottomright",
           legend = c("Empirical CDF", "Theoretical χ² CDF"),
           col = c("black", "red"),
           lty = c(1, 2),
           lwd = 2,
           bty = "n")
  }
}