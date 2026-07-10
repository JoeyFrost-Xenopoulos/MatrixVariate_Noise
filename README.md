# Ampharos

<img src="data/mde_icon_2.png" align="right" width="140" alt="Drifloon hex sticker" />

Matrix-variate Gaussian mixture models with explicit background noise
components. Supports EM and MM estimation, four initialization schemes
(k-means, k-means++, random, emrefine), two noise types (HC improper constant, BR
convex-hull uniform), and automatic noise-level selection via KS
goodness-of-fit.

## Quick Start

```r
# x_list is a list of r × p numeric matrices
fit <- matrix_variate_noise_fit(x_list, g = 3, noise_type = "hc")
print(fit$pi)
table(fit$cluster)
plot(fit$logLik, type = "b")
```

## Core Functions

| Function | Description |
|----------|-------------|
| `matrix_variate_mixture_fit()` | Standard EM mixture model (no noise) |
| `matrix_variate_noise_fit()` | Noise mixture with HC or BR noise component |

## Initialization Schemes

All fitting functions accept an `init` argument:

- **`"kmeans++"`** — K-means++ seeding (Arthur & Vassilvitskii, 2007): selects initial centers via D^2 weighting for better spread, then runs k-means. More robust than plain k-means on overlapping clusters.
- **`"emrefine"`** — Starts with random assignment, then runs a few EM iterations to refine. Balances speed and quality.

See `scripts/compare_initializations.R` for a benchmark comparing these schemes.

## Noise Types

- **HC (Hennig–Coretto):** Constant improper background density. Set `noise_k` or use `estimate_k = TRUE` for automatic selection via KS test.
- **BR (Banfield–Raftery):** Uniform density within the convex hull of the data. Requires the `geometry` package; expensive for large `r * p`.

## Usage Notes

- Use `noise_type = "hc"` for a constant improper baseline (set `noise_k`).
- Use `noise_type = "br"` to restrict noise to the convex hull of the data
  (requires `geometry`; may be expensive for large `r * p`).
- To automatically select an HC `k` value, set `estimate_k = TRUE`.