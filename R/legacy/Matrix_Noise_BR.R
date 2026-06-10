#' Construct a BR-Style Convex Hull for Matrix Noise
#'
#' @param x_list A list of same-sized matrices.
#' @param jitter Positive padding retained for compatibility.
#'
#' @return A list with the vectorized points, the convex hull object, and the
#'   log-volume of the hull.
#' @keywords internal
matrix_noise_convex_hull_support <- function(x_list, jitter = 1e-08) {
	x_list <- matrix_validate_x_list(x_list)
	noise_matrix <- do.call(rbind, lapply(x_list, function(x) as.vector(x)))
	unique_points <- unique(noise_matrix)
	point_dimension <- ncol(unique_points)

	if (nrow(unique_points) <= point_dimension) {
		stop("At least d + 1 unique vectorized matrices are required to form a convex hull.", call. = FALSE)
	}

	hull <- tryCatch(
		geometry::convhulln(unique_points, output.options = TRUE),
		error = function(e) {
			stop("Failed to build the BR noise convex hull: ", conditionMessage(e), call. = FALSE)
		}
	)

	log_volume <- if (!is.null(hull$vol)) {
		log(as.numeric(hull$vol))
	} else if (!is.null(hull$volume)) {
		log(as.numeric(hull$volume))
	} else {
		stop("The BR noise convex hull did not return a volume.", call. = FALSE)
	}

	list(
		points = unique_points,
		hull = hull,
		log_volume = log_volume,
		jitter = jitter
	)
}

#' Log Density for BR Matrix Noise
#'
#' @param x_list List of matrices to evaluate.
#' @param support Convex-hull support from `matrix_noise_convex_hull_support()`.
#'
#' @return Numeric vector of log-densities.
#' @keywords internal
matrix_noise_br_log_density <- function(x_list, support) {
	vapply(x_list, function(x) {
		vec_x <- as.vector(x)
		inside <- isTRUE(geometry::inhulln(support$hull, matrix(vec_x, nrow = 1)))
		if (inside) {
			return(-support$log_volume)
		}
		-Inf
	}, numeric(1))
}