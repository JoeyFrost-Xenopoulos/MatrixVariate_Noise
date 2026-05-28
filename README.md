# AMPH

<img src="data/mde_icon_2.png" align="right" width="140" alt="Drifloon hex sticker" />

This repository contains R code for fitting matrix-variate Gaussian mixture
models with an explicit background noise component (HC; Hennig–Coretto and
BR; Banfield–Raftery styles). The implementation supports EM fitting,
HC noise `k` grid search and Mahalanobis-based diagnostics.

Fit a model to a list of matrix observations (`x_list`):

```r
# x_list is a list of r x p numeric matrices
fit <- matrix_variate_noise_fit(x_list, g = 3, noise_type = "hc")
print(fit$pi)
table(fit$cluster)
plot(fit$logLik, type = "b")
```

Usage notes and tips
- Use `noise_type = "hc"` for a constant improper baseline (set `noise_k`).
- Use `noise_type = "br"` to restrict noise to the convex hull of the data
	(requires `geometry` and may be expensive for large `r * p`).
- To automatically select an HC `k` value, set `select_noise_k = TRUE`.
- The core utilities are in `R/Matrix.R` (`make_spd()`, `matrix_mahalanobis()`, `matrix_variate_log_density()`, and initialization via `matrix_mixture_kmeans_init()`).
