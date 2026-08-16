#' Print and Summary Methods for Fitted Model Objects
#'
#' @description
#' `print` methods display a compact overview of a fitted model.
#' `summary` methods return (and print) a detailed breakdown including
#' component sizes, mixing proportions, convergence status, noise summary,
#' information criteria, and more.
#'
#' @name print_summary
#' @rdname print_summary
#'
#' @param x A fitted model object.
#' @param object A fitted model object.
#' @param ... Additional arguments (currently unused).
#'
#' @return For `print` methods, `x` is returned invisibly.
#'   For `summary` methods, a list with class `"summary.*"` is returned.
NULL

#' @rdname print_summary
#' @export
print.mv_mixture_fit <- function(x, ...) {
  if (!inherits(x, "mv_mixture_fit")) {
    stop("'x' must be a fitted 'mv_mixture_fit' model.")
  }

  n <- nrow(x$z)
  r <- nrow(x$M[[1]])
  p <- ncol(x$M[[1]])
  g <- length(x$pi)
  sizes <- tabulate(x$cluster, nbins = g)

  cat("Matrix-Variate Gaussian Mixture Model\n")
  cat("======================================\n\n")
  cat(sprintf("Observations: %d\n", n))
  cat(sprintf("Dimensions:   %d x %d\n", r, p))
  cat(sprintf("Components:   %d\n\n", g))

  cat("Mixing Proportions:\n")
  for (i in seq_len(g)) {
    cat(sprintf("  component_%d: %.4f\n", i, x$pi[i]))
  }

  cat("\nComponent Sizes:\n")
  for (i in seq_len(g)) {
    count <- sizes[i]
    pct <- count / n * 100
    cat(sprintf("  component_%d: %d (%.1f%%)\n", i, count, pct))
  }

  cat("\nConvergence:\n")
  cat(sprintf("  Iterations:     %d\n", x$iterations))
  cat(sprintf("  Converged:      %s\n", if (x$converged) "Yes" else "No"))
  cat(sprintf("  Log-likelihood: %.4f\n", tail(x$logLik, 1)))

  invisible(x)
}

#' @rdname print_summary
#' @export
summary.mv_mixture_fit <- function(object, ...) {
  if (!inherits(object, "mv_mixture_fit")) {
    stop("'object' must be a fitted 'mv_mixture_fit' model.")
  }

  n <- nrow(object$z)
  r <- nrow(object$M[[1]])
  p <- ncol(object$M[[1]])
  g <- length(object$pi)
  sizes <- tabulate(object$cluster, nbins = g)
  names(sizes) <- paste0("component_", seq_len(g))

  n_free <- g * (r * p + r * (r + 1) / 2 + p * (p + 1) / 2 + 1) - 1
  logLik_final <- tail(object$logLik, 1)
  aic <- -2 * logLik_final + 2 * n_free
  bic <- -2 * logLik_final + log(n) * n_free

  mixing_proportions <- object$pi
  names(mixing_proportions) <- paste0("component_", seq_along(mixing_proportions))

  structure(list(
    n_obs = n,
    dimensions = c(r = r, p = p),
    n_components = g,
    mixing_proportions = mixing_proportions,
    component_sizes = sizes,
    converged = object$converged,
    iterations = object$iterations,
    logLik = logLik_final,
    aic = aic,
    bic = bic,
    n_free_params = n_free
  ), class = "summary.mv_mixture_fit")
}

#' @rdname print_summary
#' @export
print.summary.mv_mixture_fit <- function(x, ...) {
  cat("Summary of Matrix-Variate Gaussian Mixture Model\n")
  cat("=================================================\n\n")
  cat(sprintf("Observations:      %d\n", x$n_obs))
  cat(sprintf("Dimensions:        %d x %d\n", x$dimensions["r"], x$dimensions["p"]))
  cat(sprintf("Components:        %d\n", x$n_components))
  cat(sprintf("Free parameters:   %d\n\n", x$n_free_params))

  cat("Mixing Proportions:\n")
  for (i in seq_along(x$mixing_proportions)) {
    cat(sprintf("  %s: %.4f\n", names(x$mixing_proportions)[i], x$mixing_proportions[i]))
  }

  cat("\nComponent Sizes:\n")
  for (i in seq_along(x$component_sizes)) {
    pct <- x$component_sizes[i] / x$n_obs * 100
    cat(sprintf("  %s: %d (%.1f%%)\n", names(x$component_sizes)[i], x$component_sizes[i], pct))
  }

  cat("\nConvergence:\n")
  cat(sprintf("  Iterations:      %d\n", x$iterations))
  cat(sprintf("  Converged:       %s\n", if (x$converged) "Yes" else "No"))
  cat(sprintf("  Log-likelihood:  %.4f\n", x$logLik))
  cat(sprintf("  AIC:             %.4f\n", x$aic))
  cat(sprintf("  BIC:             %.4f\n", x$bic))

  invisible(x)
}

#' @rdname print_summary
#' @export
print.mv_noise_fit <- function(x, ...) {
  if (!inherits(x, "mv_noise_fit")) {
    stop("'x' must be a fitted 'mv_noise_fit' model.")
  }

  n <- nrow(x$z)
  r <- nrow(x$M[[1]])
  p <- ncol(x$M[[1]])
  g <- length(x$pi) - 1
  sizes <- tabulate(x$cluster[x$cluster > 0], nbins = g)
  noise_count <- sum(x$cluster == 0)

  cat("Matrix-Variate Gaussian Mixture Model with Noise\n")
  cat("=================================================\n\n")
  cat(sprintf("Observations: %d\n", n))
  cat(sprintf("Dimensions:   %d x %d\n", r, p))
  cat(sprintf("Components:   %d Gaussian + 1 noise\n\n", g))

  cat("Mixing Proportions:\n")
  for (i in seq_len(g)) {
    cat(sprintf("  component_%d: %.4f\n", i, x$pi[i]))
  }
  cat(sprintf("  noise:       %.4f\n", x$pi[g + 1]))

  cat("\nComponent Sizes:\n")
  for (i in seq_len(g)) {
    count <- sizes[i]
    pct <- count / n * 100
    cat(sprintf("  component_%d: %d (%.1f%%)\n", i, count, pct))
  }
  noise_pct <- noise_count / n * 100
  cat(sprintf("  noise:       %d (%.1f%%)\n", noise_count, noise_pct))

  cat("\nNoise:\n")
  cat(sprintf("  Type: %s\n", x$noise$type))
  if (x$noise$type == "hc") {
    cat(sprintf("  k:    %.4e\n", x$noise$k))
  }

  cat("\nConvergence:\n")
  cat(sprintf("  Iterations:     %d\n", x$iterations))
  cat(sprintf("  Converged:      %s\n", if (x$converged) "Yes" else "No"))
  cat(sprintf("  Log-likelihood: %.4f\n", tail(x$logLik, 1)))

  if (!is.null(x$k_selection)) {
    cat(sprintf("\nK Selection (automatic):\n"))
    cat(sprintf("  Selected k: %.4e\n", x$k_selection$selected_k))
    cat(sprintf("  Grid size:  %d\n", length(x$k_selection$k_grid)))
  }

  invisible(x)
}

#' @rdname print_summary
#' @export
summary.mv_noise_fit <- function(object, ...) {
  if (!inherits(object, "mv_noise_fit")) {
    stop("'object' must be a fitted 'mv_noise_fit' model.")
  }

  n <- nrow(object$z)
  r <- nrow(object$M[[1]])
  p <- ncol(object$M[[1]])
  g <- length(object$pi) - 1
  sizes <- tabulate(object$cluster[object$cluster > 0], nbins = g)
  names(sizes) <- paste0("component_", seq_len(g))
  noise_count <- sum(object$cluster == 0)
  component_sizes <- c(sizes, noise = noise_count)

  n_free <- g * (r * p + r * (r + 1) / 2 + p * (p + 1) / 2 + 1)
  logLik_final <- tail(object$logLik, 1)
  aic <- -2 * logLik_final + 2 * n_free
  bic <- -2 * logLik_final + log(n) * n_free

  res <- list(
    n_obs = n,
    dimensions = c(r = r, p = p),
    n_components = g,
    mixing_proportions = object$pi,
    component_sizes = component_sizes,
    converged = object$converged,
    iterations = object$iterations,
    logLik = logLik_final,
    aic = aic,
    bic = bic,
    n_free_params = n_free,
    noise_type = object$noise$type,
    noise_k = object$noise$k
  )

  if (!is.null(object$k_selection)) {
    res$k_selection <- object$k_selection
  }

  structure(res, class = "summary.mv_noise_fit")
}

#' @rdname print_summary
#' @export
print.summary.mv_noise_fit <- function(x, ...) {
  cat("Summary of Matrix-Variate Gaussian Mixture Model with Noise\n")
  cat("============================================================\n\n")
  cat(sprintf("Observations:      %d\n", x$n_obs))
  cat(sprintf("Dimensions:        %d x %d\n", x$dimensions["r"], x$dimensions["p"]))
  cat(sprintf("Components:        %d Gaussian + 1 noise\n", x$n_components))
  cat(sprintf("Free parameters:   %d\n\n", x$n_free_params))

  cat("Mixing Proportions:\n")
  for (i in seq_along(x$mixing_proportions)) {
    cat(sprintf("  %s: %.4f\n", names(x$mixing_proportions)[i], x$mixing_proportions[i]))
  }

  cat("\nComponent Sizes:\n")
  for (i in seq_along(x$component_sizes)) {
    pct <- x$component_sizes[i] / x$n_obs * 100
    cat(sprintf("  %s: %d (%.1f%%)\n", names(x$component_sizes)[i], x$component_sizes[i], pct))
  }

  cat("\nNoise:\n")
  cat(sprintf("  Type: %s\n", x$noise_type))
  if (x$noise_type == "hc") {
    cat(sprintf("  k:    %.4e\n", x$noise_k))
  }

  cat("\nConvergence:\n")
  cat(sprintf("  Iterations:      %d\n", x$iterations))
  cat(sprintf("  Converged:       %s\n", if (x$converged) "Yes" else "No"))
  cat(sprintf("  Log-likelihood:  %.4f\n", x$logLik))
  cat(sprintf("  AIC:             %.4f\n", x$aic))
  cat(sprintf("  BIC:             %.4f\n", x$bic))

  if (!is.null(x$k_selection)) {
    cat(sprintf("\nAutomatic K Selection:\n"))
    cat(sprintf("  Selected k: %.4e\n", x$k_selection$selected_k))
    cat(sprintf("  Grid size:  %d\n", length(x$k_selection$k_grid)))
  }

  invisible(x)
}
