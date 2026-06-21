#' Matrix-Variate Noise Mixture Clustering with Automatic K Selection
#'
#' Fits a matrix-variate Gaussian mixture model with a basic noise component.
#' The noise component can be either Hennig-Coretto style improper constant
#' noise (`hc`) or Banfield-Raftery style bounded uniform noise (`br`).
#'
#' @param x_list A non-empty list of numeric matrices, each of the same size.
#' @param g Integer: number of Gaussian mixture components.
#' @param noise_type Character: `"hc"` for improper constant noise or `"br"`
#'   for convex-hull uniform noise.
#' @param max_iter Integer: maximum EM iterations.
#' @param tol Numeric: convergence tolerance on the log-likelihood trace.
#' @param init Character: initialization scheme. `"kmeans"` for k-means based
#'   initialization (default), `"random"` for random cluster assignments, or
#'   `"ecme"` for ECME-style pre-EM initialization.
#' @param noise_k Numeric: constant noise height used when `noise_type = "hc"`.
#'   If `estimate_k = TRUE`, this is ignored.
#' @param estimate_k Logical: if TRUE, automatically select optimal noise_k
#'   using KS goodness-of-fit test.
#' @param k_grid Numeric vector: grid of k values to search over when
#'   `estimate_k = TRUE`. If NULL, automatically generates dimension-aware grid.
#' @param adaptive_grid Logical: if TRUE and k_grid is NULL, generate
#'   dimension-aware heuristic grid based on matrix dimensions.
#' @param noise_pi_init Numeric: initial mixing proportion for the noise
#'   component.
#' @param verbose Logical: print iteration progress.
#'
#' @return A list containing the fitted mixture parameters, posterior
#'   responsibilities, log-likelihood trace, and a noise summary. If
#'   `estimate_k = TRUE`, also includes `k_grid` and `ks_scores`.
#'
#' @export
matrix_variate_noise_fit <- function(x_list,
                                      g,
                                      noise_type = c("hc", "br"),
                                      max_iter = 1000,
                                      tol = 1e-06,
                                      nstart = 100,
                                      noise_k = 1e-04,
                                      estimate_k = FALSE,
                                      k_grid = NULL,
                                      adaptive_grid = TRUE,
                                      noise_pi_init = 0.05,
                                      init = c("kmeans", "random", "ecme"),
                                      verbose = FALSE) {
  noise_type <- match.arg(noise_type)
  init <- match.arg(init)
  x_list <- matrix_validate_x_list(x_list)
  
  # Automatic k selection for HC noise
  if (noise_type == "hc" && estimate_k) {
    if (verbose)
      cat("Selecting optimal k using KS test...\n")
    
    # Generate dimension-aware grid if not provided
    if (is.null(k_grid)) {
      if (adaptive_grid) {
        k_grid <- matrix_noise_hc_heuristic_grid(x_list)
        
        if (verbose) {
          cat(
            sprintf(
              "Using grid search: [%e, %e] with %d candidates\n",
              min(k_grid),
              max(k_grid),
              length(k_grid)
            )
          )
        }
      } else {
        k_grid <- 10^seq(-16, -1, length.out = 30)
        
        if (verbose) {
          cat("Using default fixed grid\n")
        }
      }
    }
    
    ks_scores <- rep(Inf, length(k_grid))
    all_ks_results <- vector("list", length(k_grid))
    
    
    for (i in seq_along(k_grid)) {
      current_k <- k_grid[i]
      
      if (verbose) {
        cat("  Testing k =",
            format(current_k, scientific = TRUE), "... ")
      }
      
      # Step 1: HC fit with candidate k
      fit_noise <- matrix_variate_noise_fit_impl(
        x_list = x_list,
        g = g,
        noise_type = "hc",
        max_iter = max_iter,
        tol = tol,
        nstart = nstart,
        noise_k = current_k,
        noise_jitter = NULL,
        noise_pi_init = noise_pi_init,
        init = init,
        verbose = FALSE
      )
      
      # Step 2: Remove observations classified as noise
      keep_idx <- fit_noise$cluster != 0
      x_clean <- x_list[keep_idx]
      
      # Not enough observations to refit
      if (length(x_clean) <= g) {
        ks_result <- list(
          statistic = Inf,
          p.value = NA_real_,
          n_used = length(x_clean)
        )
        
        ks_scores[i] <- Inf
        all_ks_results[[i]] <- ks_result
        
        if (verbose) {
          cat("insufficient retained observations\n")
        }
        
        next
      }
      
      # Step 3: Refit Gaussian mixture on cleaned subset
      fit_clean <- tryCatch(
        matrix_variate_mixture_fit(
          x_list = x_clean,
          g = g,
          max_iter = max_iter,
          tol = tol,
          verbose = FALSE
        ),
        error = function(e)
          NULL
      )
      
      if (is.null(fit_clean)) {
        ks_result <- list(
          statistic = Inf,
          p.value = NA_real_,
          n_used = length(x_clean)
        )
        
        ks_scores[i] <- Inf
        all_ks_results[[i]] <- ks_result
        
        if (verbose) {
          cat("refit failed\n")
        }
        
        next
      }
      
      # Step 4: KS goodness-of-fit score
      ks_result <- tryCatch(
        matrix_noise_ks_score(fit_clean, x_clean),
        error = function(e) {
          list(
            statistic = Inf,
            p.value = NA_real_,
            n_used = length(x_clean)
          )
        }
      )
      
      ks_scores[i] <- ks_result$statistic
      all_ks_results[[i]] <- ks_result
      
      if (verbose) {
        cat(sprintf(
          "KS = %.4f (n_used = %d)\n",
          ks_result$statistic,
          ks_result$n_used
        ))
      }
    }
    
    # Select optimal k
    if (all(is.infinite(ks_scores))) {
      stop("All candidate k values failed during KS selection.")
    }
    
    best_idx <- which.min(ks_scores)
    selected_k <- k_grid[best_idx]
    
    if (verbose) {
      cat(
        sprintf(
          "\nSelected optimal k = %e (KS statistic = %.4f, p-value = %.4f)\n",
          selected_k,
          ks_scores[best_idx],
          all_ks_results[[best_idx]]$p.value
        )
      )
    }
    
    # Final HC fit on FULL dataset using selected k
    final_fit <- matrix_variate_noise_fit_impl(
      x_list = x_list,
      g = g,
      noise_type = "hc",
      max_iter = max_iter,
      tol = tol,
      nstart = nstart,
      noise_k = selected_k,
      noise_jitter = NULL,
      noise_pi_init = noise_pi_init,
      init = init,
      verbose = FALSE
    )
    
    # Attach selection diagnostics
    final_fit$k_selection <- list(
      selected_k = selected_k,
      k_grid = k_grid,
      ks_scores = ks_scores,
      ks_pvalues = sapply(all_ks_results, function(x)
        x$p.value),
      n_used = sapply(all_ks_results, function(x)
        x$n_used),
      adaptive_grid = adaptive_grid
    )
    
    return(final_fit)
  }
  
  # Standard fitting (no automatic selection)
  matrix_variate_noise_fit_impl(
    x_list = x_list,
    g = g,
    noise_type = noise_type,
    max_iter = max_iter,
    tol = tol,
    nstart = nstart,
    noise_k = noise_k,
    noise_jitter = NULL,
    noise_pi_init = noise_pi_init,
    init = init,
    verbose = verbose
  )
}

matrix_variate_noise_fit_impl <- function(x_list,
                                          g,
                                          noise_type = c("hc", "br"),
                                          max_iter = 100,
                                          tol = 1e-06,
                                          nstart = 10,
                                          noise_k = 1e-04,
                                          noise_jitter = 1e-08,
                                          noise_pi_init = 0.05,
                                          init = c("kmeans", "random", "ecme"),
                                          verbose = FALSE) {
  noise_type <- match.arg(noise_type)
  init <- match.arg(init)
  
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])
  
  for (x in x_list) {
    if (!is.matrix(x) || nrow(x) != r || ncol(x) != p) {
      stop("All elements of x_list must be matrices with the same dimensions.")
    }
  }
  
  if (init == "random") {
    params <- matrix_mixture_random_init(x_list, g = g)
  } else if (init == "ecme") {
    params <- matrix_mixture_ecme_init(x_list, g = g)
  } else {
    params <- matrix_mixture_kmeans_init(x_list, g = g, nstart = nstart)
  }
  
  # For BR noise compute a convex hull over the vectorized matrices
  noise_support <- NULL
  if (noise_type == "br") {
    noise_support <- matrix_noise_convex_hull_support(x_list, jitter = noise_jitter)
  }
  
  # Append noise mixing proportion as the last component
  params$pi <- c((1 - noise_pi_init) * params$pi, noise_pi_init)
  names(params$pi) <- c(paste0("component_", seq_len(g)), "noise")
  
  loglik_trace <- numeric(0)
  responsibilities <- matrix(0, nrow = n, ncol = g + 1)
  colnames(responsibilities) <- c(paste0("component_", seq_len(g)), "noise")
  
  # Precompute noise log-density vector:
  # HC: constant improper background log(k)
  # BR: uniform within the convex hull (log(1/volume)), -Inf outside
  noise_log_density <- if (noise_type == "hc") {
    rep(log(noise_k), n)
  } else {
    matrix_noise_br_log_density(x_list, noise_support)
  }
  
  for (iteration in seq_len(max_iter)) {
    log_density <- matrix(NA_real_, nrow = n, ncol = g + 1)
    
    # E-step
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
    
    # Noise mixing proportion with noise log-density
    log_density[, g + 1] <- log(params$pi[g + 1]) + noise_log_density
    
    # Normalize log-densities to posterior responsibilities using log-sum-exp
    for (i in seq_len(n)) {
      row_log_densities <- log_density[i, ]
      normalizer <- matrix_log_sum_exp(row_log_densities)
      responsibilities[i, ] <- exp(row_log_densities - normalizer)
    }
    
    # Observed-data log-likelihood
    current_loglik <- sum(apply(log_density, 1, matrix_log_sum_exp))
    loglik_trace <- c(loglik_trace, current_loglik)
    
    if (iteration > 1 &&
        abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
      break
    }
    
    # M-step: update Gaussian parameters using responsibilities
    component_responsibilities <- responsibilities[, seq_len(g), drop = FALSE]
    component_sizes <- colSums(component_responsibilities)
    noise_size <- sum(responsibilities[, g + 1])
    new_params <- params
    
    for (component in seq_len(g)) {
      if (component_sizes[component] <= 0) {
        next
      }
      
      # Effective weights for this component
      weights <- component_responsibilities[, component]
      weights_sum <- component_sizes[component]
      v_for_row <- make_spd(params$V[[component]])
      
      # Update mean matrix: weighted average
      mean_matrix <- matrix(0, r, p)
      for (i in seq_len(n)) {
        mean_matrix <- mean_matrix + weights[i] * x_list[[i]]
      }
      mean_matrix <- mean_matrix / weights_sum
      
      # Update row covariance U_g using current V_g (v_for_row)
      row_cov <- matrix(0, r, r)
      for (i in seq_len(n)) {
        centered <- x_list[[i]] - mean_matrix
        row_cov <- row_cov + weights[i] * (centered %*% solve(v_for_row, t(centered)))
      }
      row_cov <- row_cov / (p * weights_sum)
      row_cov <- make_spd(row_cov)
      
      # Identifiability
      row_scale <- r / sum(diag(row_cov))
      row_cov <- make_spd(row_cov * row_scale)
      
      # Update column covariance V_g using updated U_g
      col_cov <- matrix(0, p, p)
      for (i in seq_len(n)) {
        centered <- x_list[[i]] - mean_matrix
        col_cov <- col_cov + weights[i] * (t(centered) %*% solve(row_cov, centered))
      }
      col_cov <- col_cov / (r * weights_sum)
      col_cov <- make_spd(col_cov)
      
      # Store updated parameters for this Gaussian component
      new_params$pi[component] <- weights_sum / n
      new_params$M[[component]] <- mean_matrix
      new_params$U[[component]] <- row_cov
      new_params$V[[component]] <- col_cov
    }
    
    # Update noise mixing proportion
    new_params$pi[g + 1] <- noise_size / n
    new_params$pi <- new_params$pi / sum(new_params$pi)
    params <- new_params
    
    if (verbose) {
      if (noise_type == "hc") {
        message(
          sprintf(
            "Iteration %d: log-likelihood = %.4f | noise_k = %.4e",
            iteration,
            current_loglik,
            noise_k
          )
        )
      } else {
        message(
          sprintf(
            "Iteration %d: log-likelihood = %.4f | noise_type = %s",
            iteration,
            current_loglik,
            noise_type
          )
        )
      }
    }
  }
  
  # Hard assignments pick the component with maximum posterior; map noise -> 0
  cluster_membership <- max.col(responsibilities, ties.method = "first")
  cluster_membership[cluster_membership == g + 1] <- 0L
  
  list(
    pi = params$pi,
    M = params$M,
    U = params$U,
    V = params$V,
    z = responsibilities,
    cluster = cluster_membership,
    logLik = loglik_trace,
    iterations = length(loglik_trace),
    converged = length(loglik_trace) < max_iter,
    noise = list(
      type = noise_type,
      pi = params$pi[g + 1],
      k = if (noise_type == "hc")
        noise_k
      else
        NA_real_,
      hull = noise_support
    )
  )
}
