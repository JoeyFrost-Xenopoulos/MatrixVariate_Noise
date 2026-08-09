#' E-Step for RIMLE ECM Algorithm
#'
#' Computes pseudo-posterior probabilities z_g for all components including
#' the improper noise component.
#'
#' @param x_list List of matrices.
#' @param params Current parameters.
#' @param g Number of Gaussian components.
#' @param n Number of observations.
#' @param k Noise constant density.
#' @return Updated params with z and Zg.
#' @noRd
rimle_e_step <- function(x_list, params, g, n, k) {
	log_density_gauss <- matrix(NA_real_, n, g)
	for (comp in seq_len(g)) {
		for (i in seq_len(n)) {
			log_density_gauss[i, comp] <- log(params$pi[comp]) +
				rimle_mv_log_density(x_list[[i]], params$M[[comp]],
						     params$U[[comp]], params$V[[comp]])
		}
	}

	log_density_noise <- log(params$pi0) + log(k)
	log_density <- cbind(log_density_gauss, matrix(log_density_noise, n, 1))

	z <- matrix(0, n, g + 1)
	for (i in seq_len(n)) {
		normalizer <- rimle_log_sum_exp(log_density[i, ])
		z[i, ] <- exp(log_density[i, ] - normalizer)
	}

	params$z <- z
	params$Zg <- colSums(z)
	params
}

#' Update V Holding U Fixed (CM1 Inner Step)
#'
#' @param x_list List of matrices.
#' @param params Current parameters.
#' @param g Number of components.
#' @param n Number of observations.
#' @param gamma_V Conditional gamma bound for V.
#' @return Updated V list.
#' @noRd
rimle_cm1_update_v <- function(x_list, params, g, n, gamma_V) {
	V_list <- vector("list", g)
	Zg <- params$Zg[seq_len(g)]

	for (comp in seq_len(g)) {
		S_V <- matrix(0, ncol(x_list[[1]]), ncol(x_list[[1]]))
		for (i in seq_len(n)) {
			centered <- x_list[[i]] - params$M[[comp]]
			S_V <- S_V + params$z[i, comp] * t(centered) %*%
				solve(make_spd(params$U[[comp]]), centered)
		}
		S_V <- S_V / (nrow(x_list[[1]]) * Zg[comp])
		V_list[[comp]] <- make_spd(S_V)
	}

	all_V_eigs <- unlist(lapply(V_list, function(V) {
		eigen(V, symmetric = TRUE, only.values = TRUE)$values
	}))

	if (max(all_V_eigs, na.rm = TRUE) / min(all_V_eigs, na.rm = TRUE) <= gamma_V) {
		return(V_list)
	}

	eig_values_list <- lapply(V_list, function(V) {
		eigen(V, symmetric = TRUE, only.values = TRUE)$values
	})

	m_star <- rimle_compute_m_star(eig_values_list, Zg, gamma_V)

	for (comp in seq_len(g)) {
		eig <- eigen(V_list[[comp]], symmetric = TRUE)
		sh <- rimle_shrinkage(eig$values, m_star, gamma_V)
		V_list[[comp]] <- eig$vectors %*% diag(sh) %*% t(eig$vectors)
	}

	V_list
}

#' Update U Holding V Fixed (CM1 Inner Step)
#'
#' @param x_list List of matrices.
#' @param params Current parameters.
#' @param g Number of components.
#' @param n Number of observations.
#' @param gamma_U Conditional gamma bound for U.
#' @return Updated U list.
#' @noRd
rimle_cm1_update_u <- function(x_list, params, g, n, gamma_U) {
	U_list <- vector("list", g)
	Zg <- params$Zg[seq_len(g)]

	for (comp in seq_len(g)) {
		S_U <- matrix(0, nrow(x_list[[1]]), nrow(x_list[[1]]))
		for (i in seq_len(n)) {
			centered <- x_list[[i]] - params$M[[comp]]
			S_U <- S_U + params$z[i, comp] *
				(centered %*% solve(make_spd(params$V[[comp]]), t(centered)))
		}
		S_U <- S_U / (ncol(x_list[[1]]) * Zg[comp])
		U_list[[comp]] <- make_spd(S_U)
	}

	all_U_eigs <- unlist(lapply(U_list, function(U) {
		eigen(U, symmetric = TRUE, only.values = TRUE)$values
	}))

	if (max(all_U_eigs, na.rm = TRUE) / min(all_U_eigs, na.rm = TRUE) <= gamma_U) {
		return(U_list)
	}

	eig_values_list <- lapply(U_list, function(U) {
		eigen(U, symmetric = TRUE, only.values = TRUE)$values
	})

	m_star <- rimle_compute_m_star(eig_values_list, Zg, gamma_U)

	for (comp in seq_len(g)) {
		eig <- eigen(U_list[[comp]], symmetric = TRUE)
		sh <- rimle_shrinkage(eig$values, m_star, gamma_U)
		U_list[[comp]] <- eig$vectors %*% diag(sh) %*% t(eig$vectors)
	}

	U_list
}

#' CM1 Step: Full Covariance Update with Flip-Flop
#'
#' Performs the nested alternating conditional maximization for U and V,
#' followed by identifiability rescaling.
#'
#' @param x_list List of matrices.
#' @param params Current parameters.
#' @param g Number of components.
#' @param n Number of observations.
#' @param gamma Eigenratio bound.
#' @param max_flip_flop Maximum flip-flop iterations.
#' @param tol Convergence tolerance.
#' @return Updated params.
#' @noRd
rimle_cm1_step <- function(x_list, params, g, n, gamma,
			    max_flip_flop = 50, tol = 1e-6) {
	r <- nrow(x_list[[1]])

	for (comp in seq_len(g)) {
		Zg <- params$Zg[comp]
		if (Zg <= 0) next
		M_new <- matrix(0, r, ncol(x_list[[1]]))
		for (i in seq_len(n)) {
			M_new <- M_new + params$z[i, comp] * x_list[[i]]
		}
		params$M[[comp]] <- M_new / Zg
	}

	prev_Q <- -Inf
	for (ff_iter in seq_len(max_flip_flop)) {
		V_old <- params$V
		U_old <- params$U

		all_U_eigs <- unlist(lapply(params$U, function(U) {
			eigen(U, symmetric = TRUE, only.values = TRUE)$values
		}))
		u_max <- max(all_U_eigs, na.rm = TRUE)
		u_min <- min(all_U_eigs, na.rm = TRUE)
		gamma_V <- gamma * u_min / u_max

		params$V <- rimle_cm1_update_v(x_list, params, g, n, gamma_V)

		all_V_eigs <- unlist(lapply(params$V, function(V) {
			eigen(V, symmetric = TRUE, only.values = TRUE)$values
		}))
		v_max <- max(all_V_eigs, na.rm = TRUE)
		v_min <- min(all_V_eigs, na.rm = TRUE)
		gamma_U <- gamma * v_min / v_max

		params$U <- rimle_cm1_update_u(x_list, params, g, n, gamma_U)

		max_change <- 0
		for (comp in seq_len(g)) {
			max_change <- max(max_change, max(abs(params$U[[comp]] - U_old[[comp]])))
			max_change <- max(max_change, max(abs(params$V[[comp]] - V_old[[comp]])))
		}

		if (max_change < tol) break

		current_Q <- rimle_compute_q1(x_list, params, g, n)
		if (ff_iter > 5 && abs(current_Q - prev_Q) < tol) break
		prev_Q <- current_Q
	}

	for (comp in seq_len(g)) {
		row_scale <- r / sum(diag(params$U[[comp]]))
		params$U[[comp]] <- make_spd(params$U[[comp]] * row_scale)
		params$V[[comp]] <- make_spd(params$V[[comp]] / row_scale)
	}

	params
}

#' CM2 Step: Update Mixing Proportions
#'
#' Updates pi_0 and pi_g subject to the noise proportion constraint.
#'
#' @param x_list List of matrices.
#' @param params Current parameters.
#' @param g Number of Gaussian components.
#' @param n Number of observations.
#' @param k Noise constant.
#' @param pi_max Maximum noise proportion.
#' @return Updated params with pi.
#' @noRd
rimle_cm2_step <- function(x_list, params, g, n, k, pi_max) {
	Zg <- params$Zg
	Z0 <- Zg[g + 1]

	if (Z0 <= n * pi_max + 1e-8) {
		params$pi0 <- Z0 / n
		params$pi <- Zg[seq_len(g)] / n
	} else {
		omega_star <- rimle_compute_omega_star(x_list, params, g, n, k, Z0, pi_max)
		params$pi0 <- omega_star
		params$pi <- (1 - omega_star) / (n - Z0) * Zg[seq_len(g)]
	}

	params$pi0 <- max(0, params$pi0)
	params$pi <- pmax(0, params$pi)

	total_pi <- params$pi0 + sum(params$pi)
	if (total_pi > 0) {
		params$pi0 <- params$pi0 / total_pi
		params$pi <- params$pi / total_pi
	}

	params
}
