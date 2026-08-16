#' Enforce Positive Definiteness on a Matrix
#'
#' Converts a matrix to symmetric positive definite form using iterative jittering
#' of the diagonal. This is necessary for numerical stability when computing
#' Cholesky decompositions and matrix inverses.
#'
#' @param mat A numeric matrix to be made positive definite
#' @param jitter Initial jitter amount added to diagonal (default: 1e-8)
#' @param max_tries Maximum number of jittering attempts (default: 8)
#'
#' @return A symmetric positive definite matrix
#'
#' @details
#' The function:
#' 1. Symmetrizes the matrix by averaging with its transpose
#' 2. Attempts Cholesky decomposition with increasing jitter amounts
#' 3. Returns the first successful candidate or errors if max_tries exceeded
#'
#' @export
#' @noRd
make_spd <- function(mat, jitter = 1e-8, max_tries = 8) {
	if (!is.matrix(mat) || !is.numeric(mat)) {
		stop("'mat' must be a numeric matrix.")
	}
	if (nrow(mat) != ncol(mat)) {
		stop("'mat' must be a square matrix.")
	}
	mat <- (mat + t(mat)) / 2
	for (k in 0:max_tries) {
		j <- jitter * (10^k)
		candidate <- mat + diag(j, nrow(mat))
		ok <- tryCatch({
			chol(candidate)
			TRUE
		}, error = function(e) FALSE)
		if (ok) return(candidate)
	}
	stop("Could not make covariance matrix positive definite.")
}
