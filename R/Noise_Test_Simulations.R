# Tomarchio and Viroli noise simulation helpers
#' @keywords internal

# Tomarchio-style helpers
generate_tomarchio_noise_means <- function() {
	M1 <- matrix(c(
		-2.60, -1.10, -0.50, -0.20,
		 1.30,  0.60,  0.30,  0.10
	), nrow = 2, byrow = TRUE)

	M2 <- matrix(c(
		 1.50,  1.70,  1.90,  2.20,
		-3.70, -2.70, -2.00, -1.50
	), nrow = 2, byrow = TRUE)

	list(M1 = M1, M2 = M2)
}

generate_tomarchio_noise_covariances <- function() {
	U1 <- matrix(c(2.0, 0.0,
	               0.0, 1.0), nrow = 2, byrow = TRUE)
	U2 <- matrix(c(1.70, 0.5,
	               0.5, 1.30), nrow = 2, byrow = TRUE)
	V <- matrix(c(
		1.00, 0.50, 0.25, 0.13,
		0.50, 1.00, 0.50, 0.25,
		0.25, 0.50, 1.00, 0.50,
		0.13, 0.25, 0.50, 1.00
	), nrow = 4, byrow = TRUE)

	list(U1 = make_spd(U1), U2 = make_spd(U2), V = make_spd(V))
}

simulate_tomarchio_noise_data <- function(n = 200, seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	means <- generate_tomarchio_noise_means()
	cov_list <- generate_tomarchio_noise_covariances()
	group_sizes <- c(floor(n / 2), n - floor(n / 2))
	true_groups <- rep(1:2, times = group_sizes)
	x_list <- vector("list", n)
	component_covs <- list(cov_list$U1, cov_list$U2)
	chol_V <- chol(cov_list$V)

	idx <- 1
	for (g in 1:2) {
		chol_U <- chol(component_covs[[g]])
		for (i in seq_len(group_sizes[g])) {
			z <- matrix(rnorm(2 * 4), nrow = 2, ncol = 4)
			x_list[[idx]] <- means[[g]] + t(chol_U) %*% z %*% chol_V
			idx <- idx + 1
		}
	}

	list(X = x_list, true_groups = true_groups)
}

contaminate_selected_columns <- function(x_list, n_contam = 30, low = -15, high = 15, seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	n <- length(x_list)
	if (n_contam > n) {
		stop("n_contam cannot exceed the number of matrices in x_list")
	}

	contam_idx <- sort(sample(seq_len(n), n_contam))
	contam_cols <- sample(seq_len(ncol(x_list[[1]])), n_contam, replace = TRUE)

	for (j in seq_along(contam_idx)) {
		i <- contam_idx[j]
		col_idx <- contam_cols[j]
		x_list[[i]][, col_idx] <- runif(nrow(x_list[[i]]), low, high)
	}

	list(X = x_list, contam_idx = contam_idx, contam_cols = contam_cols)
}

two_group_accuracy <- function(pred, truth) {
	if (length(pred) == 0) {
		return(NA_real_)
	}
	acc1 <- mean(pred == truth)
	acc2 <- mean(3 - pred == truth)
	max(acc1, acc2)
}

run_tomarchio_noise_simulation <- function(N_sim = 100, n = 200, contam_n = 0, verbose = TRUE) {
	metrics <- vector("list", N_sim * 2)
	metric_idx <- 1

	for (sim in seq_len(N_sim)) {
		if (verbose && sim %% 10 == 0) {
			cat(sprintf("Tomarchio noise simulation %d/%d\n", sim, N_sim))
		}

		dataset <- simulate_tomarchio_noise_data(n = n, seed = 9000 + sim)
		if (contam_n > 0) {
			tmp <- contaminate_selected_columns(
				dataset$X,
				n_contam = contam_n,
				low = -15,
				high = 15,
				seed = 12000 + sim
			)
			dataset$X <- tmp$X
			dataset$contam_idx <- tmp$contam_idx
			dataset$contam_cols <- tmp$contam_cols
		} else {
			dataset$contam_idx <- integer(0)
			dataset$contam_cols <- integer(0)
		}

		for (method in c("hc", "br")) {
			if (method == "hc") {
				fit <- matrix_variate_noise_fit(
					x_list = dataset$X,
					g = 2,
					noise_type = method,
					select_noise_k = TRUE,
					noise_jitter = 1e-08,
					nstart = 50,
					max_iter = 50,
					verbose = FALSE
				)
			} else {
				fit <- matrix_variate_noise_fit(
					x_list = dataset$X,
					g = 2,
					noise_type = method,
					noise_k = 1e-06,
					noise_jitter = 1e-08,
					nstart = 50,
					max_iter = 50,
					verbose = FALSE
				)
			}

			keep_idx <- which(fit$cluster > 0)
			if (length(keep_idx) > 0) {
				cluster_accuracy <- two_group_accuracy(fit$cluster[keep_idx], dataset$true_groups[keep_idx])
			} else {
				cluster_accuracy <- NA_real_
			}

			noise_rate <- mean(fit$cluster == 0)
			false_noise_rate <- if (contam_n > 0) {
				clean_idx <- setdiff(seq_len(n), dataset$contam_idx)
				mean(fit$cluster[clean_idx] == 0)
			} else {
				noise_rate
			}
			contam_recall <- if (contam_n > 0) {
				mean(fit$cluster[dataset$contam_idx] == 0)
			} else {
				NA_real_
			}
			ari <- mclust::adjustedRandIndex(fit$cluster, dataset$true_groups)

			metrics[[metric_idx]] <- data.frame(
				sim = sim,
				method = method,
				noise_rate = noise_rate,
				false_noise_rate = false_noise_rate,
				contam_recall = contam_recall,
				cluster_accuracy = cluster_accuracy,
				ari = ari
			)
			metric_idx <- metric_idx + 1
		}
	}

	metrics_df <- do.call(rbind, metrics)
	method_levels <- sort(unique(metrics_df$method))

	summary_df <- data.frame(
		mean_noise_rate = sapply(method_levels, function(m) mean(metrics_df$noise_rate[metrics_df$method == m])),
		sd_noise_rate = sapply(method_levels, function(m) sd(metrics_df$noise_rate[metrics_df$method == m])),
		mean_false_noise_rate = sapply(method_levels, function(m) mean(metrics_df$false_noise_rate[metrics_df$method == m])),
		mean_contam_recall = sapply(method_levels, function(m) {
			vals <- metrics_df$contam_recall[metrics_df$method == m]
			if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
		}),
		mean_cluster_accuracy = sapply(method_levels, function(m) mean(metrics_df$cluster_accuracy[metrics_df$method == m], na.rm = TRUE)),
		mean_ari = sapply(method_levels, function(m) mean(metrics_df$ari[metrics_df$method == m], na.rm = TRUE))
	)

	rownames(summary_df) <- method_levels

	list(
		metrics = metrics_df,
		summary = summary_df,
		contam_n = contam_n
	)
}

# Viroli-style helpers
generate_viroli_noise_means <- function() {
	M1 <- matrix(0, nrow = 3, ncol = 5)
	M1[1, 1] <- 0.5
	M1[2, 1] <- 0.5

	M2 <- matrix(0, nrow = 3, ncol = 5)

	M3 <- matrix(0, nrow = 3, ncol = 5)
	M3[1, 1] <- -0.5
	M3[2, 1] <- 0.5

	list(M1 = M1, M2 = M2, M3 = M3)
}

generate_viroli_noise_corrmatrix <- function(p) {
	make_spd(clusterGeneration::rcorrmatrix(p))
}

generate_viroli_noise_dataset <- function(n = 300, group_props = c(0.3, 0.4, 0.3), seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	means_list <- generate_viroli_noise_means()
	group_sizes <- as.integer(n * group_props)
	if (sum(group_sizes) != n) {
		stop("group_props must multiply to integer group sizes that sum to n")
	}

	true_groups <- rep(seq_along(group_sizes), times = group_sizes)
	x_list <- vector("list", n)
	U_list <- lapply(seq_along(group_sizes), function(g) generate_viroli_noise_corrmatrix(3))
	V_list <- lapply(seq_along(group_sizes), function(g) generate_viroli_noise_corrmatrix(5))

	idx <- 1
	for (g in seq_along(group_sizes)) {
		chol_U <- chol(U_list[[g]])
		chol_V <- chol(V_list[[g]])
		for (i in seq_len(group_sizes[g])) {
			z <- matrix(rnorm(3 * 5), nrow = 3, ncol = 5)
			x_list[[idx]] <- means_list[[g]] + t(chol_U) %*% z %*% chol_V
			idx <- idx + 1
		}
	}

	list(
		X = x_list,
		true_groups = true_groups,
		group_sizes = group_sizes,
		U_list = U_list,
		V_list = V_list
	)
}

permute_matrix_entries <- function(x) {
	matrix(sample(as.vector(x)), nrow = nrow(x), ncol = ncol(x))
}

contaminate_viroli_matrices <- function(x_list, n_contam = 15, seed = NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	n <- length(x_list)
	if (n_contam > n) {
		stop("n_contam cannot exceed the number of matrices in x_list")
	}

	contam_idx <- sort(sample(seq_len(n), n_contam))
	for (i in contam_idx) {
		x_list[[i]] <- permute_matrix_entries(x_list[[i]])
	}

	list(X = x_list, contam_idx = contam_idx)
}

best_three_group_accuracy <- function(pred, truth) {
	keep_idx <- which(pred > 0)
	if (length(keep_idx) == 0) {
		return(NA_real_)
	}

	pred_keep <- pred[keep_idx]
	truth_keep <- truth[keep_idx]
	perms <- list(
		c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
		c(2, 3, 1), c(3, 1, 2), c(3, 2, 1)
	)
	max(sapply(perms, function(perm) mean(perm[pred_keep] == truth_keep)))
}

run_viroli_noise_simulation <- function(N_sim = 100, n = 300, contam_n = 0, verbose = TRUE) {
	results <- vector("list", N_sim * 2)
	result_idx <- 1

	for (sim in seq_len(N_sim)) {
		if (verbose && sim %% 10 == 0) {
			cat(sprintf("Viroli noise simulation %d/%d\n", sim, N_sim))
		}

		dataset <- generate_viroli_noise_dataset(n = n, seed = 14000 + sim)
		if (contam_n > 0) {
			tmp <- contaminate_viroli_matrices(dataset$X, n_contam = contam_n, seed = 17000 + sim)
			dataset$X <- tmp$X
			dataset$contam_idx <- tmp$contam_idx
		} else {
			dataset$contam_idx <- integer(0)
		}

		for (method in c("hc", "br")) {
			if (method == "hc") {
				fit <- matrix_variate_noise_fit(
					x_list = dataset$X,
					g = 3,
					noise_type = method,
					select_noise_k = TRUE,
					noise_jitter = 1e-08,
					nstart = 50,
					max_iter = 50,
					verbose = FALSE
				)
			} else {
				fit <- matrix_variate_noise_fit(
					x_list = dataset$X,
					g = 3,
					noise_type = method,
					noise_k = 1e-06,
					noise_jitter = 1e-08,
					nstart = 50,
					max_iter = 50,
					verbose = FALSE
				)
			}

			noise_rate <- mean(fit$cluster == 0)
			false_noise_rate <- if (contam_n > 0) {
				clean_idx <- setdiff(seq_len(n), dataset$contam_idx)
				mean(fit$cluster[clean_idx] == 0)
			} else {
				noise_rate
			}
			contam_recall <- if (contam_n > 0) {
				mean(fit$cluster[dataset$contam_idx] == 0)
			} else {
				NA_real_
			}
			cluster_accuracy <- best_three_group_accuracy(fit$cluster, dataset$true_groups)
			ari <- mclust::adjustedRandIndex(fit$cluster, dataset$true_groups)

			results[[result_idx]] <- data.frame(
				sim = sim,
				method = method,
				noise_rate = noise_rate,
				false_noise_rate = false_noise_rate,
				contam_recall = contam_recall,
				cluster_accuracy = cluster_accuracy,
				ari = ari
			)
			result_idx <- result_idx + 1
		}
	}

	results_df <- do.call(rbind, results)
	method_levels <- sort(unique(results_df$method))

	summary_df <- data.frame(
		mean_noise_rate = sapply(method_levels, function(m) mean(results_df$noise_rate[results_df$method == m])),
		sd_noise_rate = sapply(method_levels, function(m) sd(results_df$noise_rate[results_df$method == m])),
		mean_false_noise_rate = sapply(method_levels, function(m) mean(results_df$false_noise_rate[results_df$method == m])),
		mean_contam_recall = sapply(method_levels, function(m) {
			vals <- results_df$contam_recall[results_df$method == m]
			if (all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
		}),
		mean_cluster_accuracy = sapply(method_levels, function(m) mean(results_df$cluster_accuracy[results_df$method == m], na.rm = TRUE)),
		mean_ari = sapply(method_levels, function(m) mean(results_df$ari[results_df$method == m], na.rm = TRUE))
	)

	rownames(summary_df) <- method_levels

	list(metrics = results_df, summary = summary_df, contam_n = contam_n)
}
