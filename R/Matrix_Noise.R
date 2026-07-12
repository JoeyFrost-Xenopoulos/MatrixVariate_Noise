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
#' @param init Character: initialization scheme. `"kmeans"` for the
#'   k-means++ seeded initialization (default), `"emrefine"` for EM-refine
#'   pre-initialization, or `"dbscan"` for DBSCAN-seeded initialization.
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
mv_noise_fit <- function(x_list,
                                      g,
                                      noise_type = c("hc", "br"),
                                      max_iter = 100,
                                      tol = 1e-06,
                                      nstart = 100,
                                      noise_k = 1e-04,
                                      estimate_k = FALSE,
                                      k_grid = NULL,
                                      adaptive_grid = TRUE,
                                      noise_pi_init = 0.05,
                                      init = c("kmeans", "emrefine", "dbscan"),
                                      verbose = FALSE) {
  noise_type <- match.arg(noise_type)
  init <- match.arg(init)
  x_list <- mv_validate_x_list(x_list)

  if (!is.numeric(g) || length(g) != 1 || g < 1) {
    stop("'g' must be a positive integer specifying the number of mixture components.")
  }
  g <- as.integer(g)

  n <- length(x_list)
  if (n < g) {
    stop(sprintf(
      "Number of observations (%d) must be at least as large as the number of components (%d).",
      n, g
    ))
  }

  if (noise_type == "hc" && !estimate_k) {
    if (!is.numeric(noise_k) || length(noise_k) != 1 || noise_k <= 0) {
      stop("'noise_k' must be a positive numeric scalar.")
    }
  }

  # Automatic k selection for HC noise
  if (noise_type == "hc" && estimate_k) {
    if (verbose)
      cat("Selecting optimal k using %s KS test...\n")
    
    # Generate dimension-aware grid if not provided
    if (is.null(k_grid)) {
      if (adaptive_grid) {
        k_grid <- mv_noise_hc_heuristic_grid(x_list)
        
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
      fit_noise <- mv_noise_fit_impl(
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
        mv_mixture_fit(
          x_list = x_clean,
          g = g,
          max_iter = max_iter,
          tol = tol,
          verbose = FALSE
        ),
        error = function(e) {
          warning(sprintf(
            "Gaussian mixture refit failed for k = %e: %s",
            current_k, conditionMessage(e)
          ), call. = FALSE)
          NULL
        }
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
        mv_noise_ks_score(fit_clean, x_clean),
        error = function(e) {
          warning(sprintf(
            "KS scoring failed for k = %e: %s",
            current_k, conditionMessage(e)
          ), call. = FALSE)
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
    final_fit <- mv_noise_fit_impl(
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
  mv_noise_fit_impl(
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

#' Core Matrix-Variate Noise Mixture Fit
#'
#' Internal implementation of the EM loop for the matrix-variate mixture with a
#' noise component. Called by `mv_noise_fit()`.
#'
#' @noRd
mv_noise_fit_impl <- function(x_list,
                                          g,
                                          noise_type = c("hc", "br"),
                                          max_iter = 100,
                                          tol = 1e-06,
                                          nstart = 10,
                                          noise_k = 1e-04,
                                          noise_jitter = 1e-08,
                                          noise_pi_init = 0.05,
                                          init = c("kmeans", "emrefine", "dbscan"),
                                          verbose = FALSE) {
  noise_type <- match.arg(noise_type)
  init <- match.arg(init)
  x_list <- mv_validate_x_list(x_list)
  
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])
  
  params <- mv_init_dispatch(x_list, g, init, nstart)
  
  # For BR noise compute a convex hull over the vectorized matrices
  noise_support <- NULL
  if (noise_type == "br") {
    noise_support <- mv_noise_convex_hull_support(x_list, jitter = noise_jitter)
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
    mv_noise_br_log_density(x_list, noise_support)
  }
  
  for (iteration in seq_len(max_iter)) {
    # E-step: Gaussian components
    log_density_gauss <- mv_e_step_log_density(x_list, params, g, n)
    log_density <- cbind(log_density_gauss, log(params$pi[g + 1]) + noise_log_density)
    
    responsibilities <- mv_normalize_responsibilities(log_density)
    
    # Observed-data log-likelihood
    current_loglik <- sum(apply(log_density, 1, mv_log_sum_exp))
    loglik_trace <- c(loglik_trace, current_loglik)
    
    if (iteration > 1 &&
        abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
      break
    }
    
    # M-step
    component_responsibilities <- responsibilities[, seq_len(g), drop = FALSE]
    component_sizes <- colSums(component_responsibilities)
    noise_size <- sum(responsibilities[, g + 1])
    new_params <- params
    
    for (component in seq_len(g)) {
      if (component_sizes[component] <= 0) {
        warning(sprintf(
          "Component %d has zero effective membership at iteration %d; skipping update.",
          component, iteration
        ), call. = FALSE)
        next
      }
      
      weights <- component_responsibilities[, component]
      weights_sum <- component_sizes[component]
      
      mean_matrix <- mv_weighted_mean(x_list, weights, weights_sum, r, p)
      row_cov <- mv_update_row_cov(x_list, mean_matrix, params$V[[component]],
                                       weights, weights_sum, r, p)
      col_cov <- mv_update_col_cov(x_list, mean_matrix, row_cov,
                                       weights, weights_sum, r, p)
      
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
        message(sprintf("Iteration %d: log-likelihood = %.4f | noise_k = %.4e",
                        iteration, current_loglik, noise_k))
      } else {
        message(sprintf("Iteration %d: log-likelihood = %.4f | noise_type = %s",
                        iteration, current_loglik, noise_type))
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
