#' Validate Matrix List Input
#'
#' @param x_list List of matrices to validate.
#' @return A list of same-sized matrices.
#' @keywords internal
matrix_validate_x_list <- function(x_list) {
	if (!is.list(x_list) || length(x_list) == 0) {
		stop("x_list must be a non-empty list of matrices.")
	}

	r <- nrow(x_list[[1]])
	p <- ncol(x_list[[1]])
	for (x in x_list) {
		if (!is.matrix(x) || nrow(x) != r || ncol(x) != p) {
			stop("All elements of x_list must be matrices with the same dimensions.")
		}
	}

	x_list
}

#' Stable Log-Sum-Exp
#'
#' @param values Numeric vector.
#' @return Numeric scalar.
#' @keywords internal
matrix_log_sum_exp <- function(values) {
	finite_values <- values[is.finite(values)]
	if (length(finite_values) == 0) {
		return(-Inf)
	}
	max_value <- max(finite_values)
	max_value + log(sum(exp(finite_values - max_value)))
}