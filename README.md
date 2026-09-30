# SSLfmm

`SSLfmm` is an R package for semi-supervised Gaussian finite mixture models with partially observed class labels. It supports complete-case, MCAR, entropy-dependent MAR, and mixed MCAR/MAR analyses. In the mixed formulation, the source of a missing label may be known or unknown. The package provides a common workflow for model fitting, simulation, prediction, classification performance assessment, and entropy-based diagnostics.

## User-facing API

- `fit_sslfmm()` — fit `cc`, `mcar`, `mar`, or `mixed` models.
- `initialize_sslfmm()` — stable parameter initialization.
- `rmix()` — simple Gaussian finite-mixture generator.
- `simulate_sslfmm()` — simulate all four label-observation mechanisms.
- `simulate_mixed_missingness()` — convenience wrapper for mixed missingness.
- `predict()` — classes, posterior probabilities, entropy, or all three for an `SSLfmm` fit.
- `classification_performance()` — classification metrics and confusion matrix.

Low-level likelihood, parameter-packing, Cholesky, and entropy helpers are internal and intentionally not exported.

## Covariance input

Simulation accepts either:

- one shared `p x p` symmetric positive-definite matrix; or
- one `p x p x g` array of symmetric positive-definite component covariance matrices.

For `p = 1`, a length-one scalar is also accepted as a shared variance. Matrix and array inputs are validated explicitly, including dimensions, finite values, symmetry, and positive definiteness.

Fitting supports `covariance_type = "equal"` and `covariance_type = "unequal"` throughout initialization, likelihood fitting, and prediction.

## Stable simulation return format

`simulate_sslfmm()` and `simulate_mixed_missingness()` always return exactly five top-level components:

```r
c("data", "true_setup", "groups", "probs", "raw")
```

The leading `data` columns are kept in a stable documented order:

```r
x1, ..., xp, en, missing, label, truth
```

The current package then adds explicit fields:

```r
observed_missing, latent_missing, missing_source, prob_mar, entropy
```

`en` is identical to `entropy`, and `missing` is identical to `observed_missing`. In simulation, `latent_missing` is the true MCAR-channel trigger.

`groups` begins with:

```r
mar_group, obs_group, mcar_in_mar, mcar_in_obs
```

and additionally includes directly useful `observed`, `mcar`, `mar`, and `missing` row indices.

## Mixed missingness indicators

For `fit_sslfmm(method = "mixed")`:

- `indicator = "unknown"`: the MCAR/MAR source of each missing label is unknown and `alpha` is estimated jointly.
- `indicator = "known"`: supply the known source of each missing label via `missing_source` (`"mcar"` / `"mar"`, or a logical/0-1 MCAR indicator).

An unknown-source fit stores `latent_missing_probability`, the fitted posterior probability that a missing label came through the MCAR channel.

## Minimal example

```r
mu <- matrix(c(-1, 1), nrow = 1, ncol = 2)
sim <- simulate_mixed_missingness(
  n = 200,
  pi = c(0.5, 0.5),
  mu = mu,
  sigma = matrix(1, 1, 1),
  seed = 1
)

x <- as.matrix(sim$data["x1"])
fit <- fit_sslfmm(
  x, sim$data$label,
  g = 2,
  method = "mixed",
  covariance_type = "equal",
  indicator = "unknown",
  n_starts = 5,
  seed = 2
)

predict(fit, x[1:10, , drop = FALSE], type = "posterior")
classification_performance(
  sim$data$truth,
  predict(fit, x),
  predict(fit, x, type = "posterior")
)
```


## Included case-study data

Version 0.2.2 includes the semi-synthetic `blood_transfusion`
data set used in the software-paper application. It can be loaded directly
from the package:

```r
library(SSLfmm)

data("blood_transfusion")
head(blood_transfusion)
table(blood_transfusion$missing_indicator)
```

The complete reference labels are retained for evaluation only; the partially
observed response is stored in `observed`.

## Development repository and issue tracker

Development is version controlled in the public GitHub repository:

- https://github.com/wujrtudou/SSLfmm
- Issues: https://github.com/wujrtudou/SSLfmm/issues

The CRAN release remains the recommended installation source for most users.

## Installation and checking

Install a built source tarball with:

```r
install.packages("SSLfmm_0.2.2.tar.gz", repos = NULL, type = "source")
```

Or install an unpacked source directory from a shell with:

```sh
R CMD INSTALL SSLfmm
```

For formal validation:

```sh
R CMD build SSLfmm
R CMD check SSLfmm_0.2.2.tar.gz --as-cran
```

