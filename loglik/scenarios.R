library(Ampharos)
library(clusterGeneration)

# ──────────────────────────────────────────────────────────────────────────────
# Simulation helpers
# ──────────────────────────────────────────────────────────────────────────────
viroli_simulation <- function(r, p, n, n1, n2, n3, M1, M2, M3,
                               U1, V1, U2, V2, U3, V3, n_outliers,
                               row_cov_scale = 1, col_cov_scale = 1) {
  U1 <- row_cov_scale * U1; V1 <- col_cov_scale * V1
  U2 <- row_cov_scale * U2; V2 <- col_cov_scale * V2
  U3 <- row_cov_scale * U3; V3 <- col_cov_scale * V3
  simulate_matrix_group <- function(n, mean_mat, row_cov, col_cov) {
    lapply(seq_len(n), function(i) {
      mean_mat + row_cov %*% matrix(rnorm(r * p), r, p) %*% col_cov
    })
  }
  x_list <- c(
    simulate_matrix_group(n1, M1, U1, V1),
    simulate_matrix_group(n2, M2, U2, V2),
    simulate_matrix_group(n3, M3, U3, V3)
  )
  outlier_idx <- sample(length(x_list), n_outliers)
  for (idx in outlier_idx) {
    x_list[[idx]] <- matrix(sample(x_list[[idx]]), r, p)
  }
  list(x_list = x_list, outlier_idx = outlier_idx)
}

scaled_viroli_simulation <- function(r, p, n, n1, n2, n3, n_outliers,
                                     signal_strength = 0.5, cov_scale = 1,
                                     outlier_type = c("perm", "row_spike", "col_out"),
                                     seed = NULL) {
  outlier_type <- match.arg(outlier_type)
  if (!is.null(seed)) set.seed(seed)
  M1 <- matrix(0, r, p)
  M2 <- matrix(0, r, p)
  M3 <- matrix(0, r, p)
  signal_col <- min(2, p)
  for (j in seq_len(signal_col)) {
    M1[1, j] <- signal_strength
    M3[1, j] <- -signal_strength
  }
  U1 <- cov_scale * rcorrmatrix(r)
  U2 <- cov_scale * rcorrmatrix(r)
  U3 <- cov_scale * rcorrmatrix(r)
  V1 <- cov_scale * rcorrmatrix(p)
  V2 <- cov_scale * rcorrmatrix(p)
  V3 <- cov_scale * rcorrmatrix(p)
  result <- viroli_simulation(r, p, n, n1, n2, n3,
                              M1, M2, M3, U1, V1, U2, V2, U3, V3,
                              n_outliers = n_outliers)
  if (outlier_type == "row_spike") {
    spike_idx <- sample(seq_along(result$x_list), n_outliers)
    for (idx in spike_idx) {
      row_id <- sample.int(r, 1)
      result$x_list[[idx]][row_id, ] <- runif(p, -8, 8)
    }
  }
  if (outlier_type == "col_out") {
    col_out_idx <- sample(seq_along(result$x_list), n_outliers)
    for (idx in col_out_idx) {
      col_id <- sample.int(p, 1)
      result$x_list[[idx]][, col_id] <- runif(r, -10, 10)
    }
  }
  list(x_list = result$x_list, outlier_idx = result$outlier_idx)
}



# ══════════════════════════════════════════════════════════════════════════════
# Base scenario definitions
# ══════════════════════════════════════════════════════════════════════════════

# ─── Viroli-style scenarios (3×5, 3 groups, permutation outliers) ─────────────
v_r <- 3; v_p <- 5; v_n <- 300

M1_v_base <- matrix(c( 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = v_r, ncol = v_p, byrow = FALSE)
M2_v_base <- matrix(0, v_r, v_p)
M3_v_base <- matrix(c(-0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = v_r, ncol = v_p, byrow = FALSE)

U1v <- rcorrmatrix(v_r); V1v <- rcorrmatrix(v_p)
U2v <- rcorrmatrix(v_r); V2v <- rcorrmatrix(v_p)
U3v <- rcorrmatrix(v_r); V3v <- rcorrmatrix(v_p)

# Scenario V1: small signal (first-column means 0.05, 0 −0.05)
sim_v1 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = matrix(c(0.05, 0.05, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  M2 = M2_v_base, M3 = matrix(c(-0.05, 0.05, 0, rep(0, 12)),
                               nrow = v_r, ncol = v_p, byrow = FALSE),
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V2: reduced covariance (row/col noise scaled to 0.3)
U1v2 <- 0.3 * rcorrmatrix(v_r); V1v2 <- 0.3 * rcorrmatrix(v_p)
U2v2 <- 0.3 * rcorrmatrix(v_r); V2v2 <- 0.3 * rcorrmatrix(v_p)
U3v2 <- 0.3 * rcorrmatrix(v_r); V3v2 <- 0.3 * rcorrmatrix(v_p)
sim_v2 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v2, V1 = V1v2, U2 = U2v2, V2 = V2v2, U3 = U3v2, V3 = V3v2,
  n_outliers = 15
)

# Scenario V3: fewer permutation outliers (5 instead of 15)
set.seed(42 + 1)
sim_v3 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 5
)

# Scenario V4: higher overlap (smaller means + larger covariance)
U1v4 <- 1.2 * rcorrmatrix(v_r); V1v4 <- 1.2 * rcorrmatrix(v_p)
U2v4 <- 1.2 * rcorrmatrix(v_r); V2v4 <- 1.2 * rcorrmatrix(v_p)
U3v4 <- 1.2 * rcorrmatrix(v_r); V3v4 <- 1.2 * rcorrmatrix(v_p)
sim_v4 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = matrix(c(0.15, 0.15, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  M2 = M2_v_base,
  M3 = matrix(c(-0.15, 0.15, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  U1 = U1v4, V1 = V1v4, U2 = U2v4, V2 = V2v4, U3 = U3v4, V3 = V3v4,
  n_outliers = 15
)

# ─── Viroli-style: highly unequal 3-group proportions ────────────────────────
set.seed(42 + 200)
v_n_unequal <- 300
n1_uneq <- round(0.10 * v_n_unequal)   # 30
n2_uneq <- round(0.10 * v_n_unequal)   # 30
n3_uneq <- v_n_unequal - n1_uneq - n2_uneq   # 240
sim_v5 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n_unequal,
  n1 = n1_uneq, n2 = n2_uneq, n3 = n3_uneq,
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V6: Viroli with row-spike outliers (instead of permutation)
set.seed(42 + 201)
sim_v6 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 0
)
spike_idx <- sample(seq_along(sim_v6$x_list), 15)
for (idx in spike_idx) {
  row_id <- sample.int(v_r, 1)
  sim_v6$x_list[[idx]][row_id, ] <- runif(v_p, -8, 8)
}
sim_v6$outlier_idx <- spike_idx

# Scenario V7: Viroli with column outliers (one column replaced per outlier)
set.seed(42 + 202)
sim_v7 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 0
)
col_out_idx <- sample(seq_along(sim_v7$x_list), 15)
for (idx in col_out_idx) {
  col_id <- sample.int(v_p, 1)
  sim_v7$x_list[[idx]][, col_id] <- runif(v_r, -10, 10)
}
sim_v7$outlier_idx <- col_out_idx

# ─── Viroli extensions: signal strength, covariance scale, outlier counts ─────
# Scenario V8: strong signal (means 0.8 instead of 0.5)
set.seed(42 + 203)
sim_v8 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = matrix(c(0.8, 0.8, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  M2 = M2_v_base,
  M3 = matrix(c(-0.8, 0.8, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V9: weak signal (means 0.2)
set.seed(42 + 204)
sim_v9 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = matrix(c(0.2, 0.2, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  M2 = M2_v_base,
  M3 = matrix(c(-0.2, 0.2, 0, rep(0, 12)), nrow = v_r, ncol = v_p, byrow = FALSE),
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V10: very high covariance (x1.5 scale)
U1v10 <- 1.5 * rcorrmatrix(v_r); V1v10 <- 1.5 * rcorrmatrix(v_p)
U2v10 <- 1.5 * rcorrmatrix(v_r); V2v10 <- 1.5 * rcorrmatrix(v_p)
U3v10 <- 1.5 * rcorrmatrix(v_r); V3v10 <- 1.5 * rcorrmatrix(v_p)
set.seed(42 + 205)
sim_v10 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v10, V1 = V1v10, U2 = U2v10, V2 = V2v10, U3 = U3v10, V3 = V3v10,
  n_outliers = 15
)

# Scenario V11: very low covariance (x0.3 scale)
U1v11 <- 0.3 * rcorrmatrix(v_r); V1v11 <- 0.3 * rcorrmatrix(v_p)
U2v11 <- 0.3 * rcorrmatrix(v_r); V2v11 <- 0.3 * rcorrmatrix(v_p)
U3v11 <- 0.3 * rcorrmatrix(v_r); V3v11 <- 0.3 * rcorrmatrix(v_p)
set.seed(42 + 206)
sim_v11 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v11, V1 = V1v11, U2 = U2v11, V2 = V2v11, U3 = U3v11, V3 = V3v11,
  n_outliers = 15
)

# Scenario V12: many outliers (25)
set.seed(42 + 207)
sim_v12 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 25
)

# Scenario V13: few outliers (3)
set.seed(42 + 208)
sim_v13 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 3
)

# Scenario V14: mixed outlier types (5 perm + 5 row spike + 5 col)
set.seed(42 + 209)
sim_v14 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 0
)
perm_idx_v14 <- sample(seq_along(sim_v14$x_list), 5)
for (idx in perm_idx_v14) {
  sim_v14$x_list[[idx]] <- matrix(sample(sim_v14$x_list[[idx]]), v_r, v_p)
}
spike_idx_v14 <- sample(setdiff(seq_along(sim_v14$x_list), perm_idx_v14), 5)
for (idx in spike_idx_v14) {
  row_id <- sample.int(v_r, 1)
  sim_v14$x_list[[idx]][row_id, ] <- runif(v_p, -8, 8)
}
col_idx_v14 <- sample(setdiff(seq_along(sim_v14$x_list), c(perm_idx_v14, spike_idx_v14)), 5)
for (idx in col_idx_v14) {
  col_id <- sample.int(v_p, 1)
  sim_v14$x_list[[idx]][, col_id] <- runif(v_r, -10, 10)
}
sim_v14$outlier_idx <- sort(union(union(perm_idx_v14, spike_idx_v14), col_idx_v14))

# Scenario V15: extreme group imbalance (50/50/200)
set.seed(42 + 210)
sim_v15 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = 50, n2 = 50, n3 = 200,
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V16: larger sample size (n=300)
set.seed(42 + 211)
sim_v16 <- viroli_simulation(
  r = v_r, p = v_p, n = 300,
  n1 = round(0.3 * 300), n2 = round(0.4 * 300),
  n3 = 300 - round(0.3 * 300) - round(0.4 * 300),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V17: smaller sample size (n=300)
set.seed(42 + 212)
sim_v17 <- viroli_simulation(
  r = v_r, p = v_p, n = 300,
  n1 = round(0.3 * 300), n2 = round(0.4 * 300),
  n3 = 300 - round(0.3 * 300) - round(0.4 * 300),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v, V1 = V1v, U2 = U2v, V2 = V2v, U3 = U3v, V3 = V3v,
  n_outliers = 15
)

# Scenario V18: heteroscedastic groups (different cov scales)
set.seed(42 + 213)
sim_v18 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = 0.5 * U1v, V1 = 0.5 * V1v,
  U2 = 1.5 * U2v, V2 = 1.5 * V2v,
  U3 = 1.0 * U3v, V3 = 1.0 * V3v,
  n_outliers = 15
)

# Scenario V19: anti-correlated row covariances
set.seed(42 + 214)
U1v19 <- rcorrmatrix(v_r); V1v19 <- rcorrmatrix(v_p)
U2v19 <- rcorrmatrix(v_r); V2v19 <- rcorrmatrix(v_p)
U3v19 <- rcorrmatrix(v_r); V3v19 <- rcorrmatrix(v_p)
U1v19[1, 2] <- -U1v19[1, 2]; U1v19[2, 1] <- U1v19[1, 2]
V1v19[1, 2] <- -V1v19[1, 2]; V1v19[2, 1] <- V1v19[1, 2]
sim_v19 <- viroli_simulation(
  r = v_r, p = v_p, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v_base, M2 = M2_v_base, M3 = M3_v_base,
  U1 = U1v19, V1 = V1v19, U2 = U2v19, V2 = V2v19, U3 = U3v19, V3 = V3v19,
  n_outliers = 15
)

# Scenario V20: larger matrix dimensions (4x6)
set.seed(42 + 215)
v_r20 <- 4; v_p20 <- 6
M1_v20 <- matrix(0, v_r20, v_p20)
M2_v20 <- matrix(0, v_r20, v_p20)
M3_v20 <- matrix(0, v_r20, v_p20)
M1_v20[1, 1] <- 0.5; M3_v20[1, 1] <- -0.5
U1v20 <- rcorrmatrix(v_r20); V1v20 <- rcorrmatrix(v_p20)
U2v20 <- rcorrmatrix(v_r20); V2v20 <- rcorrmatrix(v_p20)
U3v20 <- rcorrmatrix(v_r20); V3v20 <- rcorrmatrix(v_p20)
sim_v20 <- viroli_simulation(
  r = v_r20, p = v_p20, n = v_n,
  n1 = round(0.3 * v_n), n2 = round(0.4 * v_n),
  n3 = v_n - round(0.3*v_n) - round(0.4*v_n),
  M1 = M1_v20, M2 = M2_v20, M3 = M3_v20,
  U1 = U1v20, V1 = V1v20, U2 = U2v20, V2 = V2v20, U3 = U3v20, V3 = V3v20,
  n_outliers = 15
)

# ══════════════════════════════════════════════════════════════════════════════
# Run all scenarios
# ══════════════════════════════════════════════════════════════════════════════

results <- list()

results$V1_small_signal <- future(run_loglik_scenario(
  x_list      = sim_v1$x_list,
  g           = 3,
  scenario_name = "Viroli_small_signal",
  subtitle    = "Viroli 3×5 | 3 groups | first-col means 0.05,0,−0.05 | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v1$outlier_idx
  ), seed = TRUE)

results$V2_low_covariance <- future(run_loglik_scenario(
  x_list      = sim_v2$x_list,
  g           = 3,
  scenario_name = "Viroli_low_cov",
  subtitle    = "Viroli 3×5 | 3 groups | row/col SD = 0.3 instead of 0.5 | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v2$outlier_idx
  ), seed = TRUE)

results$V3_few_outliers <- future(run_loglik_scenario(
  x_list      = sim_v3$x_list,
  g           = 3,
  scenario_name = "Viroli_few_outliers",
  subtitle    = "Viroli 3×5 | 3 groups | 5 permuted-outliers instead of 15",
  r           = v_r, p = v_p,
  outlier_idx = sim_v3$outlier_idx
  ), seed = TRUE)

results$V4_high_overlap <- future(run_loglik_scenario(
  x_list      = sim_v4$x_list,
  g           = 3,
  scenario_name = "Viroli_high_overlap",
  subtitle    = "Viroli 3×5 | weaker means (±0.15) + larger covariance (x1.2) | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v4$outlier_idx
  ), seed = TRUE)

results$V5_extreme_imbalance <- future(run_loglik_scenario(
  x_list      = sim_v5$x_list,
  g           = 3,
  scenario_name = "Viroli_extreme_imbalance",
  subtitle    = "Viroli 3x5 | 3 groups with π=(0.10,0.10,0.80) | 15 permuted outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v5$outlier_idx
  ), seed = TRUE)

results$V6_row_spike_outliers <- future(run_loglik_scenario(
  x_list      = sim_v6$x_list,
  g           = 3,
  scenario_name = "Viroli_row_spike",
  subtitle    = "Viroli 3x5 | 3 groups | 15 row-spike outliers (one noisy row) | 0 permuted",
  r           = v_r, p = v_p,
  outlier_idx = sim_v6$outlier_idx
  ), seed = TRUE)

results$V7_column_outliers <- future(run_loglik_scenario(
  x_list      = sim_v7$x_list,
  g           = 3,
  scenario_name = "Viroli_column_outliers",
  subtitle    = "Viroli 3x5 | 3 groups | 15 column-outliers (one noisy column) | 0 permuted",
  r           = v_r, p = v_p,
  outlier_idx = sim_v7$outlier_idx
  ), seed = TRUE)

results$V8_strong_signal <- future(run_loglik_scenario(
  x_list      = sim_v8$x_list,
  g           = 3,
  scenario_name = "Viroli_strong_signal",
  subtitle    = "Viroli 3x5 | 3 groups | strong means (0.8,0,−0.8) | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v8$outlier_idx
  ), seed = TRUE)

results$V9_weak_signal <- future(run_loglik_scenario(
  x_list      = sim_v9$x_list,
  g           = 3,
  scenario_name = "Viroli_weak_signal",
  subtitle    = "Viroli 3x5 | 3 groups | weak means (0.2,0,−0.2) | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v9$outlier_idx
  ), seed = TRUE)

results$V10_very_high_cov <- future(run_loglik_scenario(
  x_list      = sim_v10$x_list,
  g           = 3,
  scenario_name = "Viroli_very_high_cov",
  subtitle    = "Viroli 3x5 | 3 groups | covariance scaled x1.5 | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v10$outlier_idx
  ), seed = TRUE)

results$V11_very_low_cov <- future(run_loglik_scenario(
  x_list      = sim_v11$x_list,
  g           = 3,
  scenario_name = "Viroli_very_low_cov",
  subtitle    = "Viroli 3x5 | 3 groups | covariance scaled x0.3 | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v11$outlier_idx
  ), seed = TRUE)

results$V12_many_outliers <- future(run_loglik_scenario(
  x_list      = sim_v12$x_list,
  g           = 3,
  scenario_name = "Viroli_many_outliers",
  subtitle    = "Viroli 3x5 | 3 groups | 25 permuted outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v12$outlier_idx
  ), seed = TRUE)

results$V13_few_outliers <- future(run_loglik_scenario(
  x_list      = sim_v13$x_list,
  g           = 3,
  scenario_name = "Viroli_few_outliers_v2",
  subtitle    = "Viroli 3x5 | 3 groups | only 3 permuted outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v13$outlier_idx
  ), seed = TRUE)

results$V14_mixed_outliers <- future(run_loglik_scenario(
  x_list      = sim_v14$x_list,
  g           = 3,
  scenario_name = "Viroli_mixed_outliers",
  subtitle    = "Viroli 3x5 | 3 groups | 5 perm + 5 row-spike + 5 col outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v14$outlier_idx
  ), seed = TRUE)

results$V15_extreme_imbalance <- future(run_loglik_scenario(
  x_list      = sim_v15$x_list,
  g           = 3,
  scenario_name = "Viroli_extreme_imbalance_v2",
  subtitle    = "Viroli 3x5 | 3 groups π=(0.17,0.17,0.67) | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v15$outlier_idx
  ), seed = TRUE)

results$V16_large_n <- future(run_loglik_scenario(
  x_list      = sim_v16$x_list,
  g           = 3,
  scenario_name = "Viroli_large_n",
  subtitle    = "Viroli 3x5 | 3 groups | n=300 | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v16$outlier_idx
  ), seed = TRUE)

results$V17_small_n <- future(run_loglik_scenario(
  x_list      = sim_v17$x_list,
  g           = 3,
  scenario_name = "Viroli_small_n",
  subtitle    = "Viroli 3x5 | 3 groups | n=300 | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v17$outlier_idx
  ), seed = TRUE)

results$V18_heteroscedastic <- future(run_loglik_scenario(
  x_list      = sim_v18$x_list,
  g           = 3,
  scenario_name = "Viroli_heteroscedastic",
  subtitle    = "Viroli 3x5 | 3 groups | heteroscedastic cov scales (0.5,1.5,1.0) | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v18$outlier_idx
  ), seed = TRUE)

results$V19_anticorrelated <- future(run_loglik_scenario(
  x_list      = sim_v19$x_list,
  g           = 3,
  scenario_name = "Viroli_anticorrelated",
  subtitle    = "Viroli 3x5 | 3 groups | anti-correlated row/col covariances | 15 outliers",
  r           = v_r, p = v_p,
  outlier_idx = sim_v19$outlier_idx
  ), seed = TRUE)

results$V20_larger_dimensions <- future(run_loglik_scenario(
  x_list      = sim_v20$x_list,
  g           = 3,
  scenario_name = "Viroli_larger_dimensions",
  subtitle    = "Viroli 4x6 | 3 groups | larger r/p | 15 outliers",
  r           = v_r20, p = v_p20,
  outlier_idx = sim_v20$outlier_idx
  ), seed = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Scaled size-sensitivity scenarios
# ══════════════════════════════════════════════════════════════════════════════

size_grid <- list(
  c(2, 3), c(2, 4), c(3, 3), c(3, 4), c(3, 5),
  c(4, 6), c(5, 7), c(5, 8), c(6, 9), c(7, 10),
  c(8, 11), c(9, 12), c(10, 13)
)

cat("\n===== Running scaled size-sensitivity grids =====\n")

for (s in seq_along(size_grid)) {
  rg <- size_grid[[s]][1]
  pg <- size_grid[[s]][2]
  n_total <- 300

  plots_dir <- file.path("loglik", "plots", sprintf("size_%dx%d", rg, pg))

  cat(sprintf("\n  >>> Size %d: r=%d, p=%d (n=%d) -> %s\n",
              s - 1, rg, pg, n_total, plots_dir))

  for (noise_prop in seq(0, 30, by = 1)) {
    prop <- noise_prop / 100
    n_outliers <- min(max(0, round(n_total * prop)), n_total - 3)

    set.seed(42 + 700 + s + noise_prop)
    sim_viroli_prop <- scaled_viroli_simulation(
      r = rg, p = pg, n = n_total,
      n1 = round(0.3 * n_total), n2 = round(0.4 * n_total),
      n3 = n_total - round(0.3 * n_total) - round(0.4 * n_total),
      n_outliers = n_outliers, signal_strength = 0.5, cov_scale = 1,
      outlier_type = "perm"
    )
    results[[sprintf("V1_prop%d_sz%d", noise_prop, s - 1)]] <- future(run_loglik_scenario(
      x_list = sim_viroli_prop$x_list, g = 3,
      scenario_name = sprintf("V1_prop%d_sz%d", noise_prop, s - 1),
      subtitle = sprintf("Scaling %d: Viroli %dx%d | 3 groups | %d permuted outliers (%d%% of n=%d)",
                         s - 1, rg, pg, n_outliers, noise_prop, n_total),
      r = rg, p = pg, true_k_noise = n_outliers,
      outlier_idx = sim_viroli_prop$outlier_idx,
      plots_dir = plots_dir
    ), seed = TRUE)
  }
}

cat("  Resolving all futures ...\n")
for (nm in names(results)) {
  results[[nm]] <- value(results[[nm]])
}

all_metrics <- rbindlist(lapply(seq_along(results), function(i) {
  nm <- names(results)[[i]]
  r  <- results[[nm]]
  if (!is.null(r$metrics)) r$metrics[, scenario_idx := i]
}), fill = TRUE)

base_metrics <- all_metrics[!grepl("^V1_prop\\d+_sz\\d+$", scenario)]
scaled_metrics <- all_metrics[grepl("^V1_prop\\d+_sz\\d+$", scenario)]

if (!dir.exists("loglik/plots")) dir.create("loglik/plots", recursive = TRUE)
if (nrow(base_metrics) > 0) {
  fwrite(base_metrics, file.path("loglik", "plots", "summary_metrics.csv"), row.names = FALSE)
  cat(sprintf("Saved base metrics to %s\n", file.path("loglik", "plots", "summary_metrics.csv")))
}

for (s in seq_along(size_grid)) {
  rg <- size_grid[[s]][1]; pg <- size_grid[[s]][2]
  sz_metrics <- scaled_metrics[grepl(sprintf("_sz%d$", s - 1), scenario)]
  if (nrow(sz_metrics) > 0) {
    plots_dir <- file.path("loglik", "plots", sprintf("size_%dx%d", rg, pg))
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
    fwrite(sz_metrics, file.path(plots_dir, "metrics.csv"), row.names = FALSE)
    cat(sprintf("Saved metrics for size %dx%d to %s\n", rg, pg, file.path(plots_dir, "metrics.csv")))
  }
}

summary_data <- rbindlist(lapply(seq_along(results), function(i) {
  nm <- names(results)[[i]]
  r  <- results[[nm]]
  data.table(
    idx           = i,
    Scenario      = nm,
    g             = 3,
    Selected_k    = r$best_k,
    Max_logLik_k  = ifelse(any(r$plot_df$logLik[complete.cases(r$plot_df$logLik)] > -Inf),
                          r$plot_df$k[complete.cases(r$plot_df$logLik)][
                            which.max(r$plot_df$logLik[complete.cases(r$plot_df$logLik)])
                          ], NA_real_),
    Spearman_k_vs_logLik = r$correlation,
    N_candidates  = nrow(r$plot_df),
    N_valid_logLik= sum(complete.cases(r$plot_df$logLik))
  )
}), fill = TRUE)

setcolorder(summary_data, c("idx", "Scenario", "g", "Selected_k", "Max_logLik_k",
                             "Spearman_k_vs_logLik", "N_candidates", "N_valid_logLik"))

cat("\n===== Cross-scenario summary =====\n")
print(summary_data)

cat(sprintf("\nSaved %d scenario trend plots under %s/\n",
            length(results), normalizePath("loglik/plots")))

cat("Plots organized by size in subfolders:\n")
for (s in seq_along(size_grid)) {
  rg <- size_grid[[s]][1]; pg <- size_grid[[s]][2]
  cat(sprintf("  size_%dx%d/\n", rg, pg))
}
