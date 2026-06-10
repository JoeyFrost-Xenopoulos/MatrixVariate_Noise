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