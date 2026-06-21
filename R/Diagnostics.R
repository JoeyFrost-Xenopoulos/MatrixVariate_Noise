#' Diagnose whether k-means WCSS is better/worse than true classes
#'
#' @param data Either a list of matrices (all same dimensions) or a numeric matrix
#'   where each row is an observation (already vectorised).
#' @param true_labels An integer vector of true class labels (1..K). Can be NULL.
#' @param kmeans_labels An integer vector of cluster assignments from k-means (1..g).
#'   If NULL, k-means is run internally with `nstart` restarts.
#' @param g Integer: number of clusters for k-means (only needed if `kmeans_labels` is NULL).
#' @param nstart Integer: number of restarts for k-means (default 10).
#' @param verbose Logical: if TRUE, prints the diagnostic message.
#'
#' @return Invisibly a list with:
#'   \item{A}{WCSS of the k-means solution}
#'   \item{B}{WCSS of the true classes (NA if not provided)}
#'   \item{recommendation}{A character string describing the outcome}
#'
#' @details
#' WCSS is the total within-cluster sum of squares (sum of squared Euclidean
#' distances of points to their cluster centroid).
#'
#' - If `B < A`, the true classes have lower WCSS → a better clustering
#'   exists. The k-means solution is likely stuck in a local optimum, so
#'   a smarter initialisation (e.g. k-means++ or increasing `nstart`) may help.
#' - If `A < B`, k-means found a tighter Euclidean clustering than the true
#'   labels. This usually means the true structure is not well captured by
#'   Euclidean distances in the current feature space, and no initialisation
#'   trick will fix that.
#'
#' @export
diagnose_kmeans_wcss <- function(data,
                                 true_labels = NULL,
                                 kmeans_labels = NULL,
                                 g = NULL,
                                 nstart = 10,
                                 verbose = TRUE) {
  
  # Convert list of matrices to vectorised data matrix if needed
  if (is.list(data)) {
    if (!all(sapply(data, is.matrix))) {
      stop("If 'data' is a list, all elements must be matrices")
    }
    dims <- lapply(data, dim)
    if (length(unique(dims)) > 1) {
      stop("All matrices must have the same dimensions")
    }
    data <- do.call(rbind, lapply(data, as.vector))
  }
  
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  n <- nrow(data)
  
  # Obtain k-means clustering if not provided
  if (is.null(kmeans_labels)) {
    if (is.null(g)) stop("Please provide either 'kmeans_labels' or 'g'")
    km <- kmeans(data, centers = g, nstart = nstart)
    kmeans_labels <- km$cluster
    A <- km$tot.withinss
  } else {
    if (length(kmeans_labels) != n) {
      stop("'kmeans_labels' length must match number of rows in 'data'")
    }
    # Compute WCSS manually for the given labels
    A <- 0
    for (cl in unique(kmeans_labels)) {
      idx <- which(kmeans_labels == cl)
      if (length(idx) <= 1) next
      centroid <- colMeans(data[idx, , drop = FALSE])
      A <- A + sum(scale(data[idx, , drop = FALSE], center = centroid, scale = FALSE)^2)
    }
  }
  
  # Compute WCSS for true classes if available
  B <- NA
  if (!is.null(true_labels)) {
    if (length(true_labels) != n) {
      stop("'true_labels' length must match number of rows in 'data'")
    }
    B <- 0
    for (cl in unique(true_labels)) {
      idx <- which(true_labels == cl)
      if (length(idx) <= 1) next
      centroid <- colMeans(data[idx, , drop = FALSE])
      B <- B + sum(scale(data[idx, , drop = FALSE], center = centroid, scale = FALSE)^2)
    }
  }
  
  # Build recommendation message
  if (is.na(B)) {
    msg <- "True class WCSS (B) not provided. Cannot compare."
    recommendation <- "unknown"
  } else if (B < A) {
    msg <- sprintf(
      "WCSS k-means (A) = %.2f > WCSS true classes (B) = %.2f.\n→ True classes have a tighter Euclidean clustering. k-means is stuck in a local optimum. You CAN likely improve by using k-means++ or increasing 'nstart'.",
      A, B
    )
    recommendation <- "try_kmeanspp"
  } else {
    msg <- sprintf(
      "WCSS k-means (A) = %.2f <= WCSS true classes (B) = %.2f.\n→ k-means found a tighter (but not necessarily more meaningful) clustering. The true structure is probably not Euclidean–separable in this vector space. Better initialisation will NOT help; reconsider features or use a different model.",
      A, B
    )
    recommendation <- "euclidean_insufficient"
  }
  
  if (verbose) cat(msg, "\n")
  
  invisible(list(
    A = A,
    B = B,
    recommendation = recommendation,
    message = msg
  ))
}

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