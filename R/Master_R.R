#' #' K-Means Initialization for Matrix Mixture Models
#' #'
#' #' @param x_list A list of numeric matrices, each of dimension r × p
#' #' @param g Integer: number of mixture components
#' #' @param nstart Integer: number of k-means restarts (default: 10)
#' #'
#' #' @return A list containing initial parameters.
#' #' @keywords internal
#' matrix_mixture_kmeans_init <- function(x_list, g, nstart = 10) {
#'   x_list <- matrix_validate_x_list(x_list)
#'   
#'   n <- length(x_list)
#'   r <- nrow(x_list[[1]])
#'   p <- ncol(x_list[[1]])
#'   
#'   # vectorize and run kmeans for init
#'   x_matrix <- do.call(rbind, lapply(x_list, function(x) as.vector(x)))
#'   km <- kmeans(x_matrix, centers = g, nstart = nstart)
#'   z <- km$cluster
#'   
#'   mixing_proportions <- numeric(g)
#'   mean_matrices <- vector("list", g)
#'   row_covariances <- vector("list", g)
#'   col_covariances <- vector("list", g)
#'   
#'   # For each component, compute sample mean and covariances from k-means clusters
#'   for (component in seq_len(g)) {
#'     component_index <- which(z == component)
#'     if (length(component_index) == 0) {
#'       component_index <- sample.int(n, 1)
#'     }
#'     
#'     component_data <- x_list[component_index]
#'     mixing_proportions[component] <- length(component_index) / n
#'     mean_matrices[[component]] <- Reduce(`+`, component_data) / length(component_data)
#'     
#'     row_cov <- matrix(0, r, r)
#'     col_cov <- matrix(0, p, p)
#'     for (x in component_data) {
#'       centered <- x - mean_matrices[[component]]
#'       row_cov <- row_cov + centered %*% t(centered)
#'       col_cov <- col_cov + t(centered) %*% centered
#'     }
#'     
#'     row_cov <- row_cov / (p * length(component_data))
#'     col_cov <- col_cov / (r * length(component_data))
#'     row_cov <- make_spd(row_cov)
#'     col_cov <- make_spd(col_cov)
#'     
#'     row_covariances[[component]] <- row_cov
#'     col_covariances[[component]] <- col_cov
#'     row_scale <- r / sum(diag(row_covariances[[component]]))
#'     row_covariances[[component]] <- row_covariances[[component]] * row_scale
#'     col_covariances[[component]] <- col_covariances[[component]] / row_scale
#'     row_covariances[[component]] <- make_spd(row_covariances[[component]])
#'     col_covariances[[component]] <- make_spd(col_covariances[[component]])
#'   }
#'   
#'   list(
#'     pi = mixing_proportions,
#'     M = mean_matrices,
#'     U = row_covariances,
#'     V = col_covariances,
#'     cluster = z
#'   )
#' }
#' 
#' #' Matrix-Variate Noise Mixture Clustering with Automatic K Selection
#' #'
#' #' Fits a matrix-variate Gaussian mixture model with a basic noise component.
#' #' The noise component can be either Hennig-Coretto style improper constant
#' #' noise (`hc`) or Banfield-Raftery style bounded uniform noise (`br`).
#' #'
#' #' @param x_list A non-empty list of numeric matrices, each of the same size.
#' #' @param g Integer: number of Gaussian mixture components.
#' #' @param noise_type Character: `"hc"` for improper constant noise or `"br"`
#' #'   for convex-hull uniform noise.
#' #' @param max_iter Integer: maximum EM iterations.
#' #' @param tol Numeric: convergence tolerance on the log-likelihood trace.
#' #' @param nstart Integer: number of k-means restarts for initialization.
#' #' @param noise_k Numeric: constant noise height used when `noise_type = "hc"`.
#' #'   If `estimate_k = TRUE`, this is ignored.
#' #' @param estimate_k Logical: if TRUE, automatically select optimal noise_k
#' #'   using KS goodness-of-fit test.
#' #' @param k_grid Numeric vector: grid of k values to search over when
#' #'   `estimate_k = TRUE`. If NULL, automatically generates dimension-aware grid.
#' #' @param adaptive_grid Logical: if TRUE and k_grid is NULL, generate
#' #'   dimension-aware heuristic grid based on matrix dimensions.
#' #' @param noise_pi_init Numeric: initial mixing proportion for the noise
#' #'   component.
#' #' @param verbose Logical: print iteration progress.
#' #'
#' #' @return A list containing the fitted mixture parameters, posterior
#' #'   responsibilities, log-likelihood trace, and a noise summary. If
#' #'   `estimate_k = TRUE`, also includes `k_grid` and `ks_scores`.
#' #'
#' #' @export
#' matrix_variate_noise_fit <- function(x_list,
#'                                      g,
#'                                      noise_type = c("hc", "br"),
#'                                      max_iter = 1000,
#'                                      tol = 1e-06,
#'                                      nstart = 100,
#'                                      noise_k = 1e-04,
#'                                      estimate_k = FALSE,
#'                                      k_grid = NULL,
#'                                      adaptive_grid = TRUE,
#'                                      noise_pi_init = 0.05,
#'                                      verbose = FALSE) {
#'   noise_type <- match.arg(noise_type)
#'   x_list <- matrix_validate_x_list(x_list)
#'   
#'   # Automatic k selection for HC noise
#'   if (noise_type == "hc" && estimate_k) {
#'     if (verbose)
#'       cat("Selecting optimal k using KS test...\n")
#'     
#'     # Generate dimension-aware grid if not provided
#'     if (is.null(k_grid)) {
#'       if (adaptive_grid) {
#'         k_grid <- matrix_noise_hc_heuristic_grid(x_list)
#'         
#'         if (verbose) {
#'           cat(
#'             sprintf(
#'               "Using grid search: [%e, %e] with %d candidates\n",
#'               min(k_grid),
#'               max(k_grid),
#'               length(k_grid)
#'             )
#'           )
#'         }
#'       } else {
#'         k_grid <- 10^seq(-16, -1, length.out = 30)
#'         
#'         if (verbose) {
#'           cat("Using default fixed grid\n")
#'         }
#'       }
#'     }
#'     
#'     ks_scores <- rep(Inf, length(k_grid))
#'     all_ks_results <- vector("list", length(k_grid))
#'     
#'     
#'     for (i in seq_along(k_grid)) {
#'       current_k <- k_grid[i]
#'       
#'       if (verbose) {
#'         cat("  Testing k =",
#'             format(current_k, scientific = TRUE), "... ")
#'       }
#'       
#'       # Step 1: HC fit with candidate k
#'       fit_noise <- matrix_variate_noise_fit_impl(
#'         x_list = x_list,
#'         g = g,
#'         noise_type = "hc",
#'         max_iter = max_iter,
#'         tol = tol,
#'         nstart = nstart,
#'         noise_k = current_k,
#'         noise_jitter = NULL,
#'         noise_pi_init = noise_pi_init,
#'         verbose = FALSE
#'       )
#'       
#'       # Step 2: Remove observations classified as noise
#'       keep_idx <- fit_noise$cluster != 0
#'       x_clean <- x_list[keep_idx]
#'       
#'       # Not enough observations to refit
#'       if (length(x_clean) <= g) {
#'         ks_result <- list(
#'           statistic = Inf,
#'           p.value = NA_real_,
#'           n_used = length(x_clean)
#'         )
#'         
#'         ks_scores[i] <- Inf
#'         all_ks_results[[i]] <- ks_result
#'         
#'         if (verbose) {
#'           cat("insufficient retained observations\n")
#'         }
#'         
#'         next
#'       }
#'       
#'       # Step 3: Refit Gaussian mixture on cleaned subset
#'       fit_clean <- tryCatch(
#'         matrix_variate_mixture_fit(
#'           x_list = x_clean,
#'           g = g,
#'           max_iter = max_iter,
#'           tol = tol,
#'           verbose = FALSE
#'         ),
#'         error = function(e)
#'           NULL
#'       )
#'       
#'       if (is.null(fit_clean)) {
#'         ks_result <- list(
#'           statistic = Inf,
#'           p.value = NA_real_,
#'           n_used = length(x_clean)
#'         )
#'         
#'         ks_scores[i] <- Inf
#'         all_ks_results[[i]] <- ks_result
#'         
#'         if (verbose) {
#'           cat("refit failed\n")
#'         }
#'         
#'         next
#'       }
#'       
#'       # Step 4: KS goodness-of-fit score
#'       ks_result <- tryCatch(
#'         matrix_noise_ks_score(fit_clean, x_clean),
#'         error = function(e) {
#'           list(
#'             statistic = Inf,
#'             p.value = NA_real_,
#'             n_used = length(x_clean)
#'           )
#'         }
#'       )
#'       
#'       ks_scores[i] <- ks_result$statistic
#'       all_ks_results[[i]] <- ks_result
#'       
#'       if (verbose) {
#'         cat(sprintf(
#'           "KS = %.4f (n_used = %d)\n",
#'           ks_result$statistic,
#'           ks_result$n_used
#'         ))
#'       }
#'     }
#'     
#'     # Select optimal k
#'     if (all(is.infinite(ks_scores))) {
#'       stop("All candidate k values failed during KS selection.")
#'     }
#'     
#'     best_idx <- which.min(ks_scores)
#'     selected_k <- k_grid[best_idx]
#'     
#'     if (verbose) {
#'       cat(
#'         sprintf(
#'           "\nSelected optimal k = %e (KS statistic = %.4f, p-value = %.4f)\n",
#'           selected_k,
#'           ks_scores[best_idx],
#'           all_ks_results[[best_idx]]$p.value
#'         )
#'       )
#'     }
#'     
#'     # Final HC fit on FULL dataset using selected k
#'     final_fit <- matrix_variate_noise_fit_impl(
#'       x_list = x_list,
#'       g = g,
#'       noise_type = "hc",
#'       max_iter = max_iter,
#'       tol = tol,
#'       nstart = nstart,
#'       noise_k = selected_k,
#'       noise_jitter = NULL,
#'       noise_pi_init = noise_pi_init,
#'       verbose = FALSE
#'     )
#'     
#'     # Attach selection diagnostics
#'     final_fit$k_selection <- list(
#'       selected_k = selected_k,
#'       k_grid = k_grid,
#'       ks_scores = ks_scores,
#'       ks_pvalues = sapply(all_ks_results, function(x)
#'         x$p.value),
#'       n_used = sapply(all_ks_results, function(x)
#'         x$n_used),
#'       adaptive_grid = adaptive_grid
#'     )
#'     
#'     return(final_fit)
#'   }
#'   
#'   # Standard fitting (no automatic selection)
#'   matrix_variate_noise_fit_impl(
#'     x_list = x_list,
#'     g = g,
#'     noise_type = noise_type,
#'     max_iter = max_iter,
#'     tol = tol,
#'     nstart = nstart,
#'     noise_k = noise_k,
#'     noise_jitter = NULL,
#'     noise_pi_init = noise_pi_init,
#'     verbose = verbose
#'   )
#' }
#' 
#' matrix_variate_noise_fit_impl <- function(x_list,
#'                                           g,
#'                                           noise_type = c("hc", "br"),
#'                                           max_iter = 100,
#'                                           tol = 1e-06,
#'                                           nstart = 10,
#'                                           noise_k = 1e-04,
#'                                           noise_jitter = 1e-08,
#'                                           noise_pi_init = 0.05,
#'                                           verbose = FALSE) {
#'   noise_type <- match.arg(noise_type)
#'   
#'   n <- length(x_list)
#'   r <- nrow(x_list[[1]])
#'   p <- ncol(x_list[[1]])
#'   
#'   for (x in x_list) {
#'     if (!is.matrix(x) || nrow(x) != r || ncol(x) != p) {
#'       stop("All elements of x_list must be matrices with the same dimensions.")
#'     }
#'   }
#'   
#'   # k-means init
#'   params <- matrix_mixture_kmeans_init(x_list, g = g, nstart = nstart)
#'   
#'   # For BR noise compute a convex hull over the vectorized matrices
#'   noise_support <- NULL
#'   if (noise_type == "br") {
#'     noise_support <- matrix_noise_convex_hull_support(x_list, jitter = noise_jitter)
#'   }
#'   
#'   # Append noise mixing proportion as the last component
#'   params$pi <- c((1 - noise_pi_init) * params$pi, noise_pi_init)
#'   names(params$pi) <- c(paste0("component_", seq_len(g)), "noise")
#'   
#'   loglik_trace <- numeric(0)
#'   responsibilities <- matrix(0, nrow = n, ncol = g + 1)
#'   colnames(responsibilities) <- c(paste0("component_", seq_len(g)), "noise")
#'   
#'   # Precompute noise log-density vector:
#'   # HC: constant improper background log(k)
#'   # BR: uniform within the convex hull (log(1/volume)), -Inf outside
#'   noise_log_density <- if (noise_type == "hc") {
#'     rep(log(noise_k), n)
#'   } else {
#'     matrix_noise_br_log_density(x_list, noise_support)
#'   }
#'   
#'   for (iteration in seq_len(max_iter)) {
#'     log_density <- matrix(NA_real_, nrow = n, ncol = g + 1)
#'     
#'     # E-step
#'     for (component in seq_len(g)) {
#'       for (i in seq_len(n)) {
#'         log_density[i, component] <- log(params$pi[component]) +
#'           matrix_variate_log_density(
#'             x = x_list[[i]],
#'             mean_matrix = params$M[[component]],
#'             row_cov = params$U[[component]],
#'             col_cov = params$V[[component]]
#'           )
#'       }
#'     }
#'     
#'     # Noise mixing proportion with noise log-density
#'     log_density[, g + 1] <- log(params$pi[g + 1]) + noise_log_density
#'     
#'     # Normalize log-densities to posterior responsibilities using log-sum-exp
#'     for (i in seq_len(n)) {
#'       row_log_densities <- log_density[i, ]
#'       normalizer <- matrix_log_sum_exp(row_log_densities)
#'       responsibilities[i, ] <- exp(row_log_densities - normalizer)
#'     }
#'     
#'     # Observed-data log-likelihood
#'     current_loglik <- sum(apply(log_density, 1, matrix_log_sum_exp))
#'     loglik_trace <- c(loglik_trace, current_loglik)
#'     
#'     if (iteration > 1 &&
#'         abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
#'       break
#'     }
#'     
#'     # M-step: update Gaussian parameters using responsibilities
#'     component_responsibilities <- responsibilities[, seq_len(g), drop = FALSE]
#'     component_sizes <- colSums(component_responsibilities)
#'     noise_size <- sum(responsibilities[, g + 1])
#'     new_params <- params
#'     
#'     for (component in seq_len(g)) {
#'       if (component_sizes[component] <= 0) {
#'         next
#'       }
#'       
#'       # Effective weights for this component
#'       weights <- component_responsibilities[, component]
#'       weights_sum <- component_sizes[component]
#'       v_for_row <- make_spd(params$V[[component]])
#'       
#'       # Update mean matrix: weighted average
#'       mean_matrix <- matrix(0, r, p)
#'       for (i in seq_len(n)) {
#'         mean_matrix <- mean_matrix + weights[i] * x_list[[i]]
#'       }
#'       mean_matrix <- mean_matrix / weights_sum
#'       
#'       # Update row covariance U_g using current V_g (v_for_row)
#'       row_cov <- matrix(0, r, r)
#'       for (i in seq_len(n)) {
#'         centered <- x_list[[i]] - mean_matrix
#'         row_cov <- row_cov + weights[i] * (centered %*% solve(v_for_row, t(centered)))
#'       }
#'       row_cov <- row_cov / (p * weights_sum)
#'       row_cov <- make_spd(row_cov)
#'       
#'       # Identifiability
#'       row_scale <- r / sum(diag(row_cov))
#'       row_cov <- make_spd(row_cov * row_scale)
#'       
#'       # Update column covariance V_g using updated U_g
#'       col_cov <- matrix(0, p, p)
#'       for (i in seq_len(n)) {
#'         centered <- x_list[[i]] - mean_matrix
#'         col_cov <- col_cov + weights[i] * (t(centered) %*% solve(row_cov, centered))
#'       }
#'       col_cov <- col_cov / (r * weights_sum)
#'       col_cov <- make_spd(col_cov)
#'       
#'       # Store updated parameters for this Gaussian component
#'       new_params$pi[component] <- weights_sum / n
#'       new_params$M[[component]] <- mean_matrix
#'       new_params$U[[component]] <- row_cov
#'       new_params$V[[component]] <- col_cov
#'     }
#'     
#'     # Update noise mixing proportion
#'     new_params$pi[g + 1] <- noise_size / n
#'     new_params$pi <- new_params$pi / sum(new_params$pi)
#'     params <- new_params
#'     
#'     if (verbose) {
#'       if (noise_type == "hc") {
#'         message(
#'           sprintf(
#'             "Iteration %d: log-likelihood = %.4f | noise_k = %.4e",
#'             iteration,
#'             current_loglik,
#'             noise_k
#'           )
#'         )
#'       } else {
#'         message(
#'           sprintf(
#'             "Iteration %d: log-likelihood = %.4f | noise_type = %s",
#'             iteration,
#'             current_loglik,
#'             noise_type
#'           )
#'         )
#'       }
#'     }
#'   }
#'   
#'   # Hard assignments pick the component with maximum posterior; map noise -> 0
#'   cluster_membership <- max.col(responsibilities, ties.method = "first")
#'   cluster_membership[cluster_membership == g + 1] <- 0L
#'   
#'   list(
#'     pi = params$pi,
#'     M = params$M,
#'     U = params$U,
#'     V = params$V,
#'     z = responsibilities,
#'     cluster = cluster_membership,
#'     logLik = loglik_trace,
#'     iterations = length(loglik_trace),
#'     converged = length(loglik_trace) < max_iter,
#'     noise = list(
#'       type = noise_type,
#'       pi = params$pi[g + 1],
#'       k = if (noise_type == "hc")
#'         noise_k
#'       else
#'         NA_real_,
#'       hull = noise_support
#'     )
#'   )
#' }
#' 
#' #' Enforce Positive Definiteness on a Matrix
#' #'
#' #' Converts a matrix to symmetric positive definite form using iterative jittering
#' #' of the diagonal. This is necessary for numerical stability when computing
#' #' Cholesky decompositions and matrix inverses.
#' #'
#' #' @param mat A numeric matrix to be made positive definite
#' #' @param jitter Initial jitter amount added to diagonal (default: 1e-8)
#' #' @param max_tries Maximum number of jittering attempts (default: 8)
#' #'
#' #' @return A symmetric positive definite matrix
#' #'
#' #' @details
#' #' The function:
#' #' 1. Symmetrizes the matrix by averaging with its transpose
#' #' 2. Attempts Cholesky decomposition with increasing jitter amounts
#' #' 3. Returns the first successful candidate or errors if max_tries exceeded
#' #'
#' #' @keywords internal
#' make_spd <- function(mat, jitter = 1e-8, max_tries = 8) {
#'   mat <- (mat + t(mat)) / 2
#'   for (k in 0:max_tries) {
#'     j <- jitter * (10^k)
#'     candidate <- mat + diag(j, nrow(mat))
#'     ok <- tryCatch({
#'       chol(candidate)
#'       TRUE
#'     }, error = function(e) FALSE)
#'     if (ok) return(candidate)
#'   }
#'   stop("Could not make covariance matrix positive definite.")
#' }
#' 
#' #' Compute Matrix-Variate Mahalanobis Distance
#' #'
#' #' Calculates the Mahalanobis distance between a matrix and a mean matrix
#' #' under the matrix-variate normal distribution with specified row and column
#' #' covariance structures.
#' #'
#' #' @param x A numeric matrix (r × p): the observation
#' #' @param mean_matrix A numeric matrix (r × p): the component mean
#' #' @param row_cov A numeric matrix (r × r): row covariance matrix U
#' #' @param col_cov A numeric matrix (p × p): column covariance matrix V
#' #'
#' #' @return Numeric scalar representing the Mahalanobis distance
#' #'
#' #' @details
#' #'
#' #' This metric extends the multivariate Mahalanobis distance to account for
#' #' the matrix structure. The computation uses Cholesky decomposition and
#' #' forward/backsolve for numerical stability.
#' #'
#' #' @keywords internal
#' matrix_mahalanobis <- function(x, mean_matrix, row_cov, col_cov) {
#'   # U^{-1} and V^{-1}
#'   row_cov <- make_spd(row_cov)
#'   col_cov <- make_spd(col_cov)
#'   row_chol <- chol(row_cov)
#'   col_chol <- chol(col_cov)
#'   centered <- x - mean_matrix
#'   row_inv_centered <- backsolve(row_chol, forwardsolve(t(row_chol), centered))
#'   col_inv <- chol2inv(col_chol)
#'   
#'   sum(row_inv_centered * (centered %*% col_inv))
#' }
#' 
#' #' Compute Log-Likelihood of Matrix under Matrix-Variate Normal Distribution
#' #'
#' #' Evaluates the log-density of a matrix observation under the matrix-variate
#' #' normal distribution with specified parameters.
#' #'
#' #' @param x A numeric matrix (r × p): the observation
#' #' @param mean_matrix A numeric matrix (r × p): the component mean matrix M
#' #' @param row_cov A numeric matrix (r × r): row covariance matrix U
#' #' @param col_cov A numeric matrix (p × p): column covariance matrix V
#' #'
#' #' @return Numeric scalar representing the log-density value
#' #'
#' #' @details
#' #'
#' #' Computation uses Cholesky decomposition for numerical stability and to
#' #' avoid explicit matrix inversion.
#' #'
#' #' @keywords internal
#' matrix_variate_log_density <- function(x, mean_matrix, row_cov, col_cov) {
#'   # Cholesky decomposition
#'   row_cov <- make_spd(row_cov)
#'   col_cov <- make_spd(col_cov)
#'   row_chol <- chol(row_cov)
#'   col_chol <- chol(col_cov)
#'   
#'   # |U| and |V| for the denominator
#'   row_logdet <- 2 * sum(log(diag(row_chol)))
#'   col_logdet <- 2 * sum(log(diag(col_chol)))
#'   
#'   # tr(V^{-1} * (X - M)^T * U^{-1} * (X - M))
#'   centered <- x - mean_matrix
#'   row_inv_centered <- backsolve(row_chol, forwardsolve(t(row_chol), centered))
#'   col_inv <- chol2inv(col_chol)
#'   trace_form <- sum(row_inv_centered * (centered %*% col_inv))
#'   
#'   r <- nrow(x)
#'   p <- ncol(x)
#'   
#'   # Returns log density value
#'   -0.5 * (r * p * log(2 * pi) + p * row_logdet + r * col_logdet + trace_form)
#' }
#' 
#' #' Fit Matrix-Variate Gaussian Mixture Model via EM Algorithm
#' #'
#' #' Estimates parameters of a matrix-variate Gaussian mixture model (MGMM)
#' #' using the Expectation-Maximization algorithm. Performs clustering of
#' #' matrix-valued observations while accounting for row and column dependencies.
#' #'
#' #' @param x_list A list of numeric matrices, each of dimension r × p
#' #' @param g Integer: number of mixture components
#' #' @param max_iter Integer: maximum EM iterations (default: 100)
#' #' @param tol Numeric: convergence tolerance for log-likelihood (default: 1e-6)
#' #' @param nstart Integer: number of k-means restarts for initialization (default: 10)
#' #' @param verbose Logical: print iteration progress (default: FALSE)
#' #'
#' #' @return A list containing:
#' #' - `pi`: numeric vector of length g with final mixing proportions.
#' #' - `M`: list of g final component mean matrices.
#' #' - `U`: list of g final row covariance matrices.
#' #' - `V`: list of g final column covariance matrices.
#' #' - `z`: numeric matrix (n × g) of posterior responsibilities.
#' #' - `cluster`: integer vector of length n with hard cluster assignments.
#' #' - `logLik`: numeric vector with the log-likelihood trace across iterations.
#' #' - `iterations`: number of EM iterations performed.
#' #' - `converged`: logical indicating whether the algorithm converged within `max_iter`.
#' #'
#' #' @details
#' #' The EM algorithm alternates between:
#' #' **E-step:** Compute posterior responsibilities (soft cluster assignments)
#' #' **M-step:** Update parameters based on responsibilities
#' #' 
#' #' @examples
#' #' \dontrun{
#' #' set.seed(123)
#' #' mean_1 <- matrix(c(1.5, 1.2, 1.0, 1.3, 1.1, 1.4, 1.2, 1.0), nrow=2)
#' #' mean_2 <- matrix(c(-1.4, -1.0, -1.2, -1.3, -1.1, -1.5, -1.0, -1.2), nrow=2)
#' #'
#' #' simulate_matrix_group <- function(n, mean_matrix, row_sd=0.35, col_sd=0.35) {
#' #'   r <- nrow(mean_matrix); p <- ncol(mean_matrix)
#' #'   row_cov <- diag(row_sd, r); col_cov <- diag(col_sd, p)
#' #'   lapply(seq_len(n), function(i) {
#' #'     noise <- matrix(rnorm(r*p), r, p)
#' #'     mean_matrix + row_cov %*% noise %*% col_cov
#' #'   })
#' #' }
#' #'
#' #' x_list <- c(
#' #'   simulate_matrix_group(15, mean_1),
#' #'   simulate_matrix_group(15, mean_2)
#' #' )
#' #'
#' #' fit <- matrix_variate_mixture_fit(x_list, g=2, max_iter=50, verbose=TRUE)
#' #' fit$cluster
#' #' fit$pi
#' #' }
#' #'
#' #' @export
#' matrix_variate_mixture_fit <- function(x_list, g, max_iter = 100, tol = 1e-06,
#'                                        nstart = 10, verbose = FALSE) {
#'   x_list <- matrix_validate_x_list(x_list)
#'   n <- length(x_list)
#'   r <- nrow(x_list[[1]])
#'   p <- ncol(x_list[[1]])
#'   
#'   # Initialize parameters using k-means
#'   params <- matrix_mixture_kmeans_init(x_list, g = g, nstart = nstart)
#'   loglik_trace <- numeric(0)
#'   responsibilities <- matrix(0, n, g)  # Will hold posterior probabilities P(z_ig = 1 | X_i)
#'   
#'   # EM loop
#'   for (iteration in seq_len(max_iter)) {
#'     # E-step: Compute responsibilities
#'     log_density <- matrix(NA_real_, nrow = n, ncol = g)
#'     
#'     for (component in seq_len(g)) {
#'       for (i in seq_len(n)) {
#'         # compute P(X_i | component)
#'         log_density[i, component] <- log(params$pi[component]) +
#'           matrix_variate_log_density(
#'             x = x_list[[i]],
#'             mean_matrix = params$M[[component]],
#'             row_cov = params$U[[component]],
#'             col_cov = params$V[[component]]
#'           )
#'       }
#'     }
#'     
#'     # Normalize log-densities to get responsibilities
#'     for (i in seq_len(n)) {
#'       row_log_densities <- log_density[i, ]
#'       normalizer <- matrix_log_sum_exp(row_log_densities)  # numerically stable log-sum-exp
#'       responsibilities[i, ] <- exp(row_log_densities - normalizer)  # z_hat_ig
#'     }
#'     
#'     # Compute observed data log-likelihood for convergence check
#'     current_loglik <- sum(apply(log_density, 1, matrix_log_sum_exp))
#'     loglik_trace <- c(loglik_trace, current_loglik)
#'     
#'     # Check convergence
#'     if (iteration > 1 && abs(loglik_trace[iteration] - loglik_trace[iteration - 1]) < tol) {
#'       break
#'     }
#'     
#'     # M-step: Update parameters
#'     component_sizes <- colSums(responsibilities)  # sum of z_hat_ig over all observations i
#'     new_params <- params
#'     
#'     # Update each component's parameters
#'     for (component in seq_len(g)) {
#'       if (component_sizes[component] <= 0) {
#'         next
#'       }
#'       
#'       weights <- responsibilities[, component]
#'       weights_sum <- component_sizes[component]  # effective sample size for this component
#'       v_for_row <- make_spd(params$V[[component]])
#'       
#'       # Update mean matrix (M-step: M_hat_g)
#'       mean_matrix <- matrix(0, r, p)
#'       for (i in seq_len(n)) {
#'         mean_matrix <- mean_matrix + weights[i] * x_list[[i]]
#'       }
#'       mean_matrix <- mean_matrix / weights_sum
#'       
#'       # Update row covariance U_hat_g
#'       row_cov <- matrix(0, r, r)
#'       for (i in seq_len(n)) {
#'         centered <- x_list[[i]] - mean_matrix
#'         row_cov <- row_cov + weights[i] * (centered %*% solve(v_for_row, t(centered)))
#'       }
#'       row_cov <- row_cov / (p * weights_sum)
#'       row_cov <- make_spd(row_cov)
#'       
#'       # Enforce tr(U) = r
#'       row_scale <- r / sum(diag(row_cov))
#'       row_cov <- row_cov * row_scale
#'       row_cov <- make_spd(row_cov)
#'       
#'       # Update column covariance V_hat_g
#'       col_cov <- matrix(0, p, p)
#'       for (i in seq_len(n)) {
#'         centered <- x_list[[i]] - mean_matrix
#'         col_cov <- col_cov + weights[i] * (t(centered) %*% solve(row_cov, centered))
#'       }
#'       col_cov <- col_cov / (r * weights_sum)
#'       col_cov <- make_spd(col_cov)
#'       
#'       # Store updated parameters
#'       new_params$pi[component] <- weights_sum / n  # mixing proportion
#'       new_params$M[[component]] <- mean_matrix
#'       new_params$U[[component]] <- row_cov
#'       new_params$V[[component]] <- col_cov
#'     }
#'     
#'     # Normalize mixing proportions to sum to 1
#'     new_params$pi <- new_params$pi / sum(new_params$pi)
#'     
#'     params <- new_params
#'     
#'     if (verbose) {
#'       message(sprintf("Iteration %d: log-likelihood = %.4f", iteration, current_loglik))
#'     }
#'   }
#'   
#'   # Assign each observation to its most likely component
#'   cluster_membership <- max.col(responsibilities, ties.method = "first")
#'   
#'   # Return fitted model with all parameters and diagnostics
#'   list(
#'     pi = params$pi,
#'     M = params$M,
#'     U = params$U,
#'     V = params$V,
#'     z = responsibilities,
#'     cluster = cluster_membership,
#'     logLik = loglik_trace,
#'     iterations = length(loglik_trace),
#'     converged = length(loglik_trace) < max_iter
#'   )
#' }
#' 
#' #' Score HC Noise Fit with a Matrix KS Test (Chi-Square)
#' #'
#' #' @param fit A fitted noise model.
#' #' @param x_list List of matrices used for fitting.
#' #' @return A list with `statistic`, `p.value`, and `n_used`.
#' #' @keywords internal
#' matrix_noise_ks_score <- function(fit, x_list) {
#'   x_list <- matrix_validate_x_list(x_list)
#'   keep_idx <- which(fit$cluster > 0)
#'   if (length(keep_idx) < 2) {
#'     return(list(statistic = Inf, p.value = NA_real_, n_used = length(keep_idx)))
#'   }
#'   
#'   distances <- vapply(keep_idx, function(i) {
#'     component <- fit$cluster[i]
#'     matrix_mahalanobis(
#'       x = x_list[[i]],
#'       mean_matrix = fit$M[[component]],
#'       row_cov = fit$U[[component]],
#'       col_cov = fit$V[[component]]
#'     )
#'   }, numeric(1))
#'   distances <- distances[is.finite(distances)]
#'   if (length(distances) < 2 || length(unique(distances)) < 2) {
#'     return(list(statistic = Inf, p.value = NA_real_, n_used = length(distances)))
#'   }
#'   
#'   # Get dimensions
#'   r <- nrow(x_list[[1]])
#'   p <- ncol(x_list[[1]])
#'   df <- r * p
#'   
#'   # One-sample KS test against Chi-squared distribution
#'   test <- tryCatch(
#'     suppressWarnings(stats::ks.test(distances, "pchisq", df = df)),
#'     error = function(e) NULL
#'   )
#'   
#'   if (is.null(test)) {
#'     return(list(statistic = Inf, p.value = NA_real_, n_used = length(distances)))
#'   }
#'   
#'   list(
#'     statistic = unname(test$statistic),
#'     p.value = unname(test$p.value),
#'     n_used = length(distances)
#'   )
#' }
#' 
#' #' Generate Dimension-Aware Heuristic Grid for HC Noise
#' #'
#' #' Creates a grid of candidate noise_k values based on matrix dimensions.
#' #' The heuristic centers the grid around 10^(-0.75 * dimension) where
#' #' dimension = rows * cols.
#' #'
#' #' @param x_list List of matrices used for fitting.
#' #' @param n_points Integer: number of points in the grid.
#' #' @return Numeric vector of candidate noise_k values.
#' #' @keywords internal
#' matrix_noise_hc_heuristic_grid <- function(x_list, n_points = 30) {
#'   x_list <- matrix_validate_x_list(x_list)
#'   
#'   dimension <- nrow(x_list[[1]]) * ncol(x_list[[1]])
#'   
#'   if (!is.finite(dimension) || dimension <= 0) {
#'     # Fallback to default grid
#'     return(10^seq(-16, -1, length.out = n_points))
#'   }
#'   
#'   # Center at -0.75 * dimension (empirical heuristic)
#'   center_log10 <- -0.75 * dimension
#'   
#'   # Width adapts to dimension: larger dimension needs wider search
#'   half_width <- max(6, ceiling(dimension / 2))
#'   
#'   # Ensure we don't go below machine precision
#'   lower_log10 <- max(log10(.Machine$double.xmin), center_log10 - half_width)
#'   upper_log10 <- center_log10 + half_width
#'   
#'   # Generate grid on log10 scale
#'   grid_log10 <- seq(lower_log10, upper_log10, length.out = n_points)
#'   grid <- 10^grid_log10
#'   
#'   # Remove any inf or NaN values
#'   grid <- grid[is.finite(grid) & grid > 0]
#'   
#'   # Ensure we have at least some points
#'   if (length(grid) < 2) {
#'     grid <- 10^seq(-16, -1, length.out = n_points)
#'   }
#'   
#'   sort(unique(grid))
#' }