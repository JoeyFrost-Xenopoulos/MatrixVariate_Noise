#' Parametric KS bootstrap vs. actual KS score (Viroli simulation)
#'
#' For each candidate k in the estimate_k grid:
#'   1. Fit the HC noise model, remove noise, refit the clean mixture.
#'   2. Collect the actual KS statistic (observation-derived distances vs chi-square).
#'   3. Generate B parametric bootstrap samples from the refitted mixture and
#'      re-run the entire k-selection chain on each bootstrap replicate,
#'      collecting the conditional bootstrap KS score for each candidate k.
#'
#' NOTE — the bootstrap null is **conditional** on the fitted mixture at each k:
#' it assesses calibration of the KS as a goodness-of-fit test for that specific
#' model, not the uncertainty in the overall k-selection procedure itself.
#' A large gap between the actual KS curve and the bootstrap null either means
#' the model is misspecified at that k or that the KS statistic has heavy tails
#' under the fitted model; both have the same diagnostic consequence: don't
#' treat the bootstrap curve as a parametric p-value for k.
#'
#' @param B Number of bootstrap replicates per k candidate
#' @param n_cores Workers for bootstrap layer
#' @export
library(Ampharos)
library(ggplot2)
library(clusterGeneration)
library(data.table)
library(future.apply)

# ──────────────────────────────────────────────────────────────────────────────
# Local matrix-variate simulation helpers (no dependency on Ampharos internals)
#
# IMPORTANT: Ampharos stores U and V as covariance matrices (not Cholesky
# factors), as confirmed by mv_log_density in Matrix_Base.R which calls
# chol(U), chol(V) and uses chol2inv / backsolve to compute the quadratic
# form.  Therefore the correct draw from MN_rp(M, U, V) is:
#
#   X = M + U^{1/2}  Z  (V^{1/2})^T,   Z ~ N(0, I_{r×p})
#
# NOT  M + U Z V,  which applies U and V as raw transformation matrices and
# produces X with an entirely wrong covariance unless U = V = I.
# ──────────────────────────────────────────────────────────────────────────────
rmv_matrix_one <- function(mean_mat, row_cov, col_cov) {
  r <- nrow(mean_mat)
  p <- ncol(mean_mat)
  Z <- matrix(rnorm(r * p), r, p)
  mean_mat + chol(row_cov) %*% Z %*% t(chol(col_cov))
}

rmv_mixture <- function(n, pi, M, U, V, rows, cols) {
  stopifnot(length(pi) == length(M))
  stopifnot(all(sapply(M, nrow) == rows))
  stopifnot(all(sapply(M, ncol) == cols))
  components <- seq_along(pi)
  assignments <- sample(components, size = n, replace = TRUE, prob = pi)
  lapply(seq_len(n), function(i) {
    comp <- assignments[i]
    rmv_matrix_one(M[[comp]], U[[comp]], V[[comp]])
  })
}

set.seed(42)

r_viroli  <- 3
p_viroli  <- 5
n_viroli  <- 300

M1_viroli <- matrix(c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = r_viroli, ncol = p_viroli, byrow = FALSE)
M2_viroli <- matrix(0, r_viroli, p_viroli)
M3_viroli <- matrix(c(-0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    nrow = r_viroli, ncol = p_viroli, byrow = FALSE)

U1_viroli <- rcorrmatrix(r_viroli); V1_viroli <- rcorrmatrix(p_viroli)
U2_viroli <- rcorrmatrix(r_viroli); V2_viroli <- rcorrmatrix(p_viroli)
U3_viroli <- rcorrmatrix(r_viroli); V3_viroli <- rcorrmatrix(p_viroli)

n1_v <- round(0.3 * n_viroli)
n2_v <- round(0.4 * n_viroli)
n3_v <- n_viroli - n1_v - n2_v

simulate_matrix_group <- function(n, mean_mat, row_cov, col_cov) {
  lapply(seq_len(n), function(i) {
    mean_mat + row_cov %*% matrix(rnorm(r_viroli * p_viroli),
                                  r_viroli, p_viroli) %*% col_cov
  })
}

x_list_viroli <- c(
  simulate_matrix_group(n1_v, M1_viroli, U1_viroli, V1_viroli),
  simulate_matrix_group(n2_v, M2_viroli, U2_viroli, V2_viroli),
  simulate_matrix_group(n3_v, M3_viroli, U3_viroli, V3_viroli)
)

outlier_idx_viroli <- sample(length(x_list_viroli), 15)
for (idx in outlier_idx_viroli) {
  x_list_viroli[[idx]] <- matrix(sample(x_list_viroli[[idx]]),
                                  r_viroli, p_viroli)
}

g   <- 3
max_iter <- 100
tol  <- 1e-6
n_boot   <- 200
verbose_boot <- FALSE

cat(sprintf("Viroli simulation: n = %d, groups = %d %d %d\n",
            length(x_list_viroli), n1_v, n2_v, n3_v))

# ──────────────────────────────────────────────────────────────────────────────
# Step 1: run estimate_k to obtain the adaptive k_grid
# ──────────────────────────────────────────────────────────────────────────────
cat("Running mv_noise_fit with estimate_k = TRUE ...\n")
fit_hc_viroli <- mv_noise_fit(
  x_list       = x_list_viroli,
  g            = g,
  noise_type   = "hc",
  nstart       = 10,
  estimate_k   = TRUE,
  init         = "kmeans",
  verbose      = TRUE,
  use_parallel = FALSE
)

k_selection <- fit_hc_viroli$k_selection
best_k      <- k_selection$selected_k

cat(sprintf("\nKS-selected k = %.4e (KS = %.4f)\n",
            best_k, k_selection$ks_scores[which.min(k_selection$ks_scores)]))

# ──────────────────────────────────────────────────────────────────────────────
# Step 2: per-candidate evaluation ─ actual KS + B bootstrap KS values
# ──────────────────────────────────────────────────────────────────────────────
k_grid <- k_selection$k_grid
cat(sprintf("Grid: %d candidates, bootstrap replicates per candidate: %d\n\n",
            length(k_grid), n_boot))

init_plan <- tryCatch(
  future::plan(future::multisession, workers = 1),
  error = function(e) future::plan(future::sequential)
)

single_candidate <- function(idx) {
  current_k <- k_grid[idx]

  fit_noise <- Ampharos:::mv_noise_fit_impl(
    x_list        = x_list_viroli,
    g             = g,
    noise_type    = "hc",
    max_iter      = max_iter,
    tol           = tol,
    nstart        = 10,
    noise_k       = current_k,
    noise_pi_init = 0.05,
    init          = "kmeans",
    verbose       = FALSE
  )

  keep_idx <- fit_noise$cluster != 0
  x_clean  <- x_list_viroli[keep_idx]
  n_clean  <- length(x_clean)

  actual_ks <- NA_real_
  actual_p  <- NA_real_
  boot_mean <- NA_real_
  boot_se   <- NA_real_
  boot_ci_lo <- NA_real_
  boot_ci_hi <- NA_real_
  boot_stats <- rep(NA_real_, n_boot)

  if (n_clean > g) {
    fit_clean <- tryCatch(
      mv_mixture_fit(x_list = x_clean, g = g,
                     max_iter = max_iter, tol = tol, verbose = FALSE),
      error = function(e) NULL
    )

    if (!is.null(fit_clean)) {
      ks_res <- suppressWarnings(
        tryCatch(
          Ampharos:::mv_noise_ks_score(fit_clean, x_clean),
          error = function(e) list(statistic = NA_real_, p.value = NA_real_,
                                    n_used = n_clean)
        )
      )
      actual_ks <- ks_res$statistic
      actual_p  <- ks_res$p.value

      for (b in seq_len(n_boot)) {
        x_boot <- tryCatch(
          rmv_mixture(
            n     = n_clean,
            pi    = fit_clean$pi,
            M     = fit_clean$M,
            U     = fit_clean$U,
            V     = fit_clean$V,
            rows  = r_viroli,
            cols  = p_viroli
          ),
          error = function(e) NULL
        )

        if (is.null(x_boot)) {
          boot_stats[b] <- NA_real_
          next
        }

        fit_boot_clean <- tryCatch(
          mv_mixture_fit(x_list = x_boot, g = g,
                         max_iter = max_iter, tol = tol, verbose = FALSE),
          error = function(e) NULL
        )

        if (is.null(fit_boot_clean)) {
          boot_stats[b] <- NA_real_
          next
        }

        ks_boot <- suppressWarnings(
          tryCatch(
            Ampharos:::mv_noise_ks_score(fit_boot_clean, x_boot),
            error = function(e) list(statistic = NA_real_, p.value = NA_real_)
          )
        )
        boot_stats[b] <- ks_boot$statistic
      }

      ok_b <- is.finite(boot_stats)
      if (any(ok_b)) {
        boot_mean  <- mean(boot_stats[ok_b])
        boot_se    <- sd(boot_stats[ok_b])
        qs <- quantile(boot_stats[ok_b], c(0.025, 0.975), names = FALSE)
        boot_ci_lo <- qs[1]
        boot_ci_hi <- qs[2]
      }
    }
  }

  list(
    k           = current_k,
    actual_ks   = actual_ks,
    actual_p    = actual_p,
    bootstrap_mean = boot_mean,
    bootstrap_se   = boot_se,
    bootstrap_ci_lo = boot_ci_lo,
    bootstrap_ci_hi = boot_ci_hi,
    n_clean     = n_clean
  )
}

cat("Computing actual KS + parametric bootstrap null for each k ...\n")
res_list <- lapply(seq_along(k_grid), single_candidate)
dt <- rbindlist(res_list)

# ──────────────────────────────────────────────────────────────────────────────
# Step 3: summarise bootstrap null by quantiles across replicates
# ──────────────────────────────────────────────────────────────────────────────
boot_mat <- dt[, .(k, actual_ks, bootstrap_mean, bootstrap_se,
                    bootstrap_ci_lo, bootstrap_ci_hi)]

boot_mat[, type := "Bootstrap mean"]
boot_mat[!is.na(actual_ks), type := ""]  # placeholder for actual curve

# Dynamic annotation
best_idx <- which.min(dt$actual_ks)

# ──────────────────────────────────────────────────────────────────────────────
# Plot 1: actual KS + bootstrap mean ± 95 % CI ribbon
# ──────────────────────────────────────────────────────────────────────────────
p_compare <- ggplot(dt, aes(x = k)) +
  geom_ribbon(aes(ymin = bootstrap_ci_lo, ymax = bootstrap_ci_hi),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = bootstrap_mean, color = "Bootstrap null (refit-inside-bootstrap)"),
            linewidth = 1) +
  geom_line(aes(y = actual_ks,   color = "Actual KS on data"), linewidth = 1.2) +
  geom_vline(xintercept = best_k, color = "black", linetype = "dashed",
             linewidth = 0.8) +
  scale_color_manual(
    name = "",
    values = c("Actual KS on data" = "darkorange2",
               "Bootstrap null (refit-inside-bootstrap)" = "steelblue")
  ) +
  scale_x_log10() +
  labs(
    title    = "Actual KS vs. conditional bootstrap null per k (Viroli simulation)",
    subtitle = sprintf(
      "Actual KS minimised at k = %.2e (KS = %.4f) | %d bootstrap replicates\n
      Bootstrap: draw from fitted mixture → refit clean model inside each replicate → recompute KS",
      dt$k[best_idx], dt$actual_ks[best_idx], n_boot
    ),
    x        = expression(k~("noise height, log scale")),
    y        = "KS statistic (lower = better)",
    caption  = "Conditional null: bootstrap tracks estimation noise, not k-selection uncertainty.\n
                A large gap ⇒ model misspecification or heavy KS tails at that k."
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.caption    = element_text(size = 8, hjust = 0)
  )

# ──────────────────────────────────────────────────────────────────────────────
# Plot 2: selected-k diagnostics — actual KS vs. bootstrap null
# ──────────────────────────────────────────────────────────────────────────────
best_row <- dt[best_idx]

p_hist <- ggplot(data.frame(ks_val = c(
  best_row$bootstrap_ci_lo, best_row$bootstrap_mean, best_row$bootstrap_ci_hi,
  best_row$actual_ks)), aes(x = ks_val)) +
  geom_segment(
    data = data.frame(xlo = best_row$bootstrap_ci_lo,
                      xhi = best_row$bootstrap_ci_hi,
                      y   = 0),
    aes(x = xlo, xend = xhi, y = y, yend = y),
    color = "steelblue", linewidth = 10, alpha = 0.3
  ) +
  geom_vline(xintercept = best_row$actual_ks,
             color = "darkorange2", linetype = "dashed", linewidth = 1.5) +
  geom_vline(xintercept = best_row$bootstrap_mean,
             color = "steelblue", linetype = "solid", linewidth = 1.2) +
  annotate("text", x = best_row$actual_ks, y = 0,
            label = sprintf("Actual KS = %.4f", best_row$actual_ks),
            color = "darkorange2", hjust = -0.05, vjust = -1, size = 3.5) +
  annotate("text", x = best_row$bootstrap_mean, y = 0,
            label = sprintf("Boot mean = %.4f", best_row$bootstrap_mean),
            color = "steelblue", hjust = -0.05, vjust =  1.1, size = 3.5) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_y_continuous(labels = NULL, breaks = NULL) +
  labs(
    title = sprintf(
      "Actual KS vs. conditional bootstrap null at selected k = %.2e",
      best_k
    ),
    subtitle = sprintf(
      "Bootstrap mean = %.4f | 95 %% CI = [%.4f, %.4f] | Actual = %.4f | p = %.4f\n
      Bootstrap replicates each refit the clean mixture before scoring.",
      best_row$bootstrap_mean,
      best_row$bootstrap_ci_lo,
      best_row$bootstrap_ci_hi,
      best_row$actual_ks,
      best_row$actual_p
    ),
    x = "KS statistic (lower = better)",
    y = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# ──────────────────────────────────────────────────────────────────────────────
# Diagnostics
# ──────────────────────────────────────────────────────────────────────────────
cat("\n--- Per-k summary at selected k ---\n")
print(dt[dt$k == best_k, ])

cat(sprintf("\nActual KS minimised at  k = %.4e\n", dt$k[best_idx]))
ok_boot <- !is.na(dt$bootstrap_mean)
if (any(ok_boot)) {
  idx_boot_best <- which.min(dt$bootstrap_mean[ok_boot])
  cat(sprintf("Bootstrap mean min at  k = %.4e\n", dt$k[ok_boot][idx_boot_best]))
}

print(p_compare)
print(p_hist)
