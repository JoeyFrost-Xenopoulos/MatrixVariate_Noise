#' K-Means Initialization for Matrix Mixture Models
#'
#' @param x_list A list of numeric matrices, each of dimension r × p
#' @param g Integer: number of mixture components
#' @param nstart Integer: number of k-means restarts (default: 10)
#'
#' @return A list containing initial parameters.
#' @keywords internal
matrix_mixture_kmeans_init <- function(x_list, g, nstart = 10) {
  x_list <- matrix_validate_x_list(x_list)

  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  # vectorize and run kmeans for init
  x_matrix <- do.call(rbind, lapply(x_list, function(x) as.vector(x)))
  km <- kmeans(x_matrix, centers = g, nstart = nstart)
  z <- km$cluster

  mixing_proportions <- numeric(g)
  mean_matrices <- vector("list", g)
  row_covariances <- vector("list", g)
  col_covariances <- vector("list", g)

  # For each component, compute sample mean and covariances from k-means clusters
  for (component in seq_len(g)) {
    component_index <- which(z == component)
    if (length(component_index) == 0) {
      component_index <- sample.int(n, 1)
    }

    component_data <- x_list[component_index]
    mixing_proportions[component] <- length(component_index) / n
    mean_matrices[[component]] <- Reduce(`+`, component_data) / length(component_data)

    row_cov <- matrix(0, r, r)
    col_cov <- matrix(0, p, p)
    for (x in component_data) {
      centered <- x - mean_matrices[[component]]
      row_cov <- row_cov + centered %*% t(centered)
      col_cov <- col_cov + t(centered) %*% centered
    }

    row_cov <- row_cov / (p * length(component_data))
    col_cov <- col_cov / (r * length(component_data))
    row_cov <- make_spd(row_cov)
    col_cov <- make_spd(col_cov)

    row_covariances[[component]] <- row_cov
    col_covariances[[component]] <- col_cov
    row_scale <- r / sum(diag(row_covariances[[component]]))
    row_covariances[[component]] <- row_covariances[[component]] * row_scale
    col_covariances[[component]] <- col_covariances[[component]] / row_scale
    row_covariances[[component]] <- make_spd(row_covariances[[component]])
    col_covariances[[component]] <- make_spd(col_covariances[[component]])
  }

  list(
    pi = mixing_proportions,
    M = mean_matrices,
    U = row_covariances,
    V = col_covariances,
    cluster = z
  )
}

matrix_mixture_random_init <- function(x_list, g) {
  x_list <- matrix_validate_x_list(x_list)
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  z <- sample(seq_len(g), size = n, replace = TRUE)
  mixing_proportions <- numeric(g)
  mean_matrices <- vector("list", g)
  row_covariances <- vector("list", g)
  col_covariances <- vector("list", g)

  for (component in seq_len(g)) {
    component_index <- which(z == component)
    if (length(component_index) == 0) {
      component_index <- sample.int(n, 1)
    }

    component_data <- x_list[component_index]
    mixing_proportions[component] <- length(component_index) / n
    mean_matrices[[component]] <- Reduce(`+`, component_data) / length(component_data)

    row_cov <- matrix(0, r, r)
    col_cov <- matrix(0, p, p)
    for (x in component_data) {
      centered <- x - mean_matrices[[component]]
      row_cov <- row_cov + centered %*% t(centered)
      col_cov <- col_cov + t(centered) %*% centered
    }

    row_cov <- row_cov / (p * length(component_data))
    col_cov <- col_cov / (r * length(component_data))
    row_cov <- make_spd(row_cov)
    col_cov <- make_spd(col_cov)

    row_covariances[[component]] <- row_cov
    col_covariances[[component]] <- col_cov
    row_scale <- r / sum(diag(row_covariances[[component]]))
    row_covariances[[component]] <- row_covariances[[component]] * row_scale
    col_covariances[[component]] <- col_covariances[[component]] / row_scale
    row_covariances[[component]] <- make_spd(row_covariances[[component]])
    col_covariances[[component]] <- make_spd(col_covariances[[component]])
  }

  list(
    pi = mixing_proportions,
    M = mean_matrices,
    U = row_covariances,
    V = col_covariances,
    cluster = z
  )
}

matrix_mixture_ecme_init <- function(x_list, g, max_iter = 5) {
  params <- matrix_mixture_random_init(x_list, g = g)
  n <- length(x_list)
  r <- nrow(x_list[[1]])
  p <- ncol(x_list[[1]])

  responsibilities <- matrix(0, n, g)
  colnames(responsibilities) <- paste0("component_", seq_len(g))

  for (iteration in seq_len(max_iter)) {
    log_density <- matrix(NA_real_, nrow = n, ncol = g)

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

    for (i in seq_len(n)) {
      row_log_densities <- log_density[i, ]
      normalizer <- matrix_log_sum_exp(row_log_densities)
      responsibilities[i, ] <- exp(row_log_densities - normalizer)
    }

    component_sizes <- colSums(responsibilities)
    new_params <- params

    for (component in seq_len(g)) {
      if (component_sizes[component] <= 0) {
        next
      }

      weights <- responsibilities[, component]
      weights_sum <- component_sizes[component]
      v_for_row <- make_spd(params$V[[component]])

      mean_matrix <- matrix(0, r, p)
      for (i in seq_len(n)) {
        mean_matrix <- mean_matrix + weights[i] * x_list[[i]]
      }
      mean_matrix <- mean_matrix / weights_sum

      row_cov <- matrix(0, r, r)
      for (i in seq_len(n)) {
        centered <- x_list[[i]] - mean_matrix
        row_cov <- row_cov + weights[i] * (centered %*% solve(v_for_row, t(centered)))
      }
      row_cov <- row_cov / (p * weights_sum)
      row_cov <- make_spd(row_cov)

      row_scale <- r / sum(diag(row_cov))
      row_cov <- make_spd(row_cov * row_scale)

      col_cov <- matrix(0, p, p)
      for (i in seq_len(n)) {
        centered <- x_list[[i]] - mean_matrix
        col_cov <- col_cov + weights[i] * (t(centered) %*% solve(row_cov, centered))
      }
      col_cov <- col_cov / (r * weights_sum)
      col_cov <- make_spd(col_cov)

      new_params$pi[component] <- weights_sum / n
      new_params$M[[component]] <- mean_matrix
      new_params$U[[component]] <- row_cov
      new_params$V[[component]] <- col_cov
    }

    new_params$pi <- new_params$pi / sum(new_params$pi)
    params <- new_params
  }

  list(
    pi = params$pi,
    M = params$M,
    U = params$U,
    V = params$V,
    cluster = max.col(responsibilities, ties.method = "first")
  )
}

