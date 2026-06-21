#' Matrix-Variate Gaussian Mixture with MM Algorithm
#'
#' Implements the MM (minorization-maximization) algorithm for matrix-variate
#' Gaussian mixture models with optional noise component. The MM surrogate ensures
#' monotonic log-likelihood increase at each iteration.
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
#' @param g Integer: number of mixture components
#' @param method Character: estimation framework. `"em"` for standard EM or
#'   `"mm"` for MM minorization-maximization.
#' @param noise_type Character: `"hc"` for improper constant noise or `"br"`
#'   for convex-hull uniform noise. Default `NULL` (no noise).
#' @param max_iter Integer: maximum iterations.
#' @param tol Numeric: convergence tolerance.
#' @param nstart Integer: number of k-means restarts for initialization.
#' @param noise_k Numeric: noise height for HC noise.
#' @param init Character: initialization scheme (`"kmeans"`, `"random"`, `"ecme"`).
#' @param verbose Logical: print iteration progress.
#'
#' @return A list containing fitted parameters, responsibilities, and diagnostics.
#' @export
matrix_mm_fit <- function(x_list,
                          g,
                          method = c("em", "mm"),
                          noise_type = c("hc", "br"),
                          max_iter = 100,
                          tol = 1e-06,
                          nstart = 10,
                          noise_k = 1e-04,
                          noise_pi_init = 0.05,
                          init = c("kmeans", "random", "ecme", "kmeans++"),
                          verbose = FALSE) {
  method <- match.arg(method)
  init <- match.arg(init)
  x_list <- matrix_validate_x_list(x_list)
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  has_noise <- !is.null(noise_type)
  if (has_noise) {
    noise_type <- match.arg(noise_type)
    n_components <- g + 1
  } else {
    n_components <- g
  }

  # Initialize
  if (init == "random") {
    params <- matrix_mixture_random_init(x_list, g = g)
  } else if (init == "ecme") {
    params <- matrix_mixture_ecme_init(x_list, g = g)
  } else {
    params <- matrix_mixture_kmeans_init(x_list, g = g, nstart = nstart)
  }

  if (has_noise) {
    params$pi <- c((1 - noise_pi_init) * params$pi, noise_pi_init)
    names(params$pi) <- c(paste0("component_", seq_len(g)), "noise")
  } else {
    names(params$pi) <- paste0("component_", seq_len(g))
    params$M <- params$M[seq_len(g)]
    params$U <- params$U[seq_len(g)]
    params$V <- params$V[seq_len(g)]
  }

  loglik_trace <- numeric(0)
  responsibilities <- matrix(0, nrow = n, ncol = n_components)
  colnames(responsibilities) <- if (has_noise) {
    c(paste0("component_", seq_len(g)), "noise")
  } else {
    paste0("component_", seq_len(g))
  }

  noise_log_density <- if (has_noise && noise_type == "hc") {
    rep(log(noise_k), n)
  }

  noise_support <- NULL
  if (has_noise && noise_type == "br") {
    noise_support <- matrix_noise_convex_hull_support(x_list, jitter = 1e-08)
    noise_log_density <- matrix_noise_br_log_density(x_list, noise_support)
  }

  is_mm <- (method == "mm")

  for (iteration in seq_len(max_iter)) {
    prev_params <- params

    # E-step: compute posteriors using CURRENT parameters
    log_density <- matrix(NA_real_, nrow = n, ncol = n_components)

    for (component in seq_len(g)) {
      for (i in seq_len(n)) {
        log_density[i, component] <- log(params$pi[component]) +
          matrix_variate_log_density(
            x = x_list[[i]],
            mean_matrix = params$M[[component]],
            row_cov = params$U[[component]],
            col_cov = params$V[[component]]
          )
      }
    }

    if (has_noise) {
      log_density[, n_components] <- log(params$pi[n_components]) + noise_log_density
    }

    for (i in seq_len(n)) {
      row_log_densities <- log_density[i, ]
      normalizer <- matrix_log_sum_exp(row_log_densities)
      responsibilities[i, ] <- exp(row_log_densities - normalizer)
    }

    current_loglik <- sum(apply(log_density, 1, matrix_log_sum_exp))
    loglik_trace <- c(loglik_trace, current_loglik)

    if (iteration > 1 &&
        abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
      break
    }

    # M-step
    component_sizes <- colSums(responsibilities)
    new_params <- params

    for (component in seq_len(g)) {
      if (component_sizes[component] <= 0) {
        next
      }

      weights <- responsibilities[, component]
      weights_sum <- component_sizes[component]

      # Mean update (same for EM and MM)
      new_M <- matrix(0, r, p)
      for (i in seq_len(n)) {
        new_M <- new_M + weights[i] * x_list[[i]]
      }
      new_M <- new_M / weights_sum

      if (is_mm) {
        # MM: use CURRENT (previous iteration) U and V for both updates
        U_curr <- make_spd(prev_params$U[[component]])
        V_curr <- make_spd(prev_params$V[[component]])

        # U update using CURRENT V
        row_cov <- matrix(0, r, r)
        for (i in seq_len(n)) {
          centered <- x_list[[i]] - new_M
          row_cov <- row_cov + weights[i] * (centered %*% solve(V_curr, t(centered)))
        }
        row_cov <- row_cov / (p * weights_sum)
        row_cov <- make_spd(row_cov)

        # V update using CURRENT U
        col_cov <- matrix(0, p, p)
        for (i in seq_len(n)) {
          centered <- x_list[[i]] - new_M
          col_cov <- col_cov + weights[i] * (t(centered) %*% solve(U_curr, centered))
        }
        col_cov <- col_cov / (r * weights_sum)
        col_cov <- make_spd(col_cov)

        new_params$U[[component]] <- row_cov
        new_params$V[[component]] <- col_cov
      } else {
        # EM: U uses CURRENT V, V uses NEW U (sequential)
        v_for_row <- make_spd(params$V[[component]])

        row_cov <- matrix(0, r, r)
        for (i in seq_len(n)) {
          centered <- x_list[[i]] - new_M
          row_cov <- row_cov + weights[i] * (centered %*% solve(v_for_row, t(centered)))
        }
        row_cov <- row_cov / (p * weights_sum)
        row_cov <- make_spd(row_cov)

        row_scale <- r / sum(diag(row_cov))
        row_cov <- make_spd(row_cov * row_scale)

        col_cov <- matrix(0, p, p)
        for (i in seq_len(n)) {
          centered <- x_list[[i]] - new_M
          col_cov <- col_cov + weights[i] * (t(centered) %*% solve(row_cov, centered))
        }
        col_cov <- col_cov / (r * weights_sum)
        col_cov <- make_spd(col_cov)

        new_params$U[[component]] <- row_cov
        new_params$V[[component]] <- col_cov
      }

      new_params$M[[component]] <- new_M
    }

    new_params$pi <- component_sizes / n
    params <- new_params

    if (verbose) {
      message(sprintf("[%s] Iteration %d: log-likelihood = %.4f",
                      ifelse(is_mm, "MM", "EM"),
                      iteration, current_loglik))
    }
  }

  cluster_membership <- max.col(responsibilities, ties.method = "first")
  if (has_noise) {
    cluster_membership[cluster_membership == n_components] <- 0L
  }

  result <- list(
    pi = params$pi,
    M = params$M,
    U = params$U,
    V = params$V,
    z = responsibilities,
    cluster = cluster_membership,
    logLik = loglik_trace,
    iterations = length(loglik_trace),
    converged = length(loglik_trace) < max_iter,
    method = method
  )

  if (has_noise) {
    result$noise <- list(
      type = noise_type,
      pi = params$pi[n_components],
      k = if (noise_type == "hc") noise_k else NA_real_,
      hull = noise_support
    )
  }

  result
}
