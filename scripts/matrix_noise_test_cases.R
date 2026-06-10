set.seed(20260610)

simulate_matrix_group <- function(n, mean_matrix, row_sd = 0.5, col_sd = 0.5, seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}
	row_cov <- diag(row_sd, nrow(mean_matrix))
	col_cov <- diag(col_sd, ncol(mean_matrix))
	lapply(seq_len(n), function(i) {
		mean_matrix + row_cov %*% matrix(rnorm(nrow(mean_matrix) * ncol(mean_matrix)), nrow(mean_matrix), ncol(mean_matrix)) %*% col_cov
	})
}

contaminate_matrix_list <- function(x_list, contam_n, seed = NULL, min_value = -15, max_value = 15) {
	if (!is.null(seed)) {
		set.seed(seed)
	}
	if (contam_n <= 0) {
		return(x_list)
	}
	contam_n <- min(contam_n, length(x_list))
	contam_idx <- sort(sample.int(length(x_list), contam_n))
	for (i in contam_idx) {
		x <- x_list[[i]]
		column_id <- sample.int(ncol(x), 1)
		x[, column_id] <- runif(nrow(x), min_value, max_value)
		x_list[[i]] <- x
	}
	x_list
}

build_matrix_noise_test_case <- function(r,
									 p,
									 n_group = 18,
									 contam_n = 4,
									 seed = 1,
									 row_sd = 0.5,
									 col_sd = 0.5) {
	set.seed(seed)
	mean_one <- matrix(1, r, p)
	mean_two <- matrix(-1, r, p)

	group_one <- simulate_matrix_group(
		n = n_group,
		mean_matrix = mean_one,
		row_sd = row_sd,
		col_sd = col_sd,
		seed = seed + 1L
	)
	group_two <- simulate_matrix_group(
		n = n_group,
		mean_matrix = mean_two,
		row_sd = row_sd,
		col_sd = col_sd,
		seed = seed + 2L
	)

	clean_list <- c(group_one, group_two)
	contaminated_list <- contaminate_matrix_list(
		x_list = clean_list,
		contam_n = contam_n,
		seed = seed + 3L
	)

	list(
		r = r,
		p = p,
		n_group = n_group,
		contam_n = contam_n,
		seed = seed,
		x_list = contaminated_list,
		true_groups = c(rep(1L, n_group), rep(2L, n_group))
	)
}

matrix_noise_test_cases <- function() {
	case_specs <- list(
		list(name = "r2_p3", r = 2, p = 3, seed = 101),
		list(name = "r3_p5", r = 3, p = 5, seed = 202),
		list(name = "r4_p6", r = 4, p = 6, seed = 303)
	)
	cases <- lapply(case_specs, function(spec) {
		build_matrix_noise_test_case(
			r = spec$r,
			p = spec$p,
			seed = spec$seed,
			n_group = 18,
			contam_n = 4
		)
	})
	names(cases) <- vapply(case_specs, `[[`, character(1), "name")
	cases
}
