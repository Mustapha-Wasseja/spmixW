# spmixW

[![CRAN status](https://www.r-pkg.org/badges/version/spmixW)](https://cran.r-project.org/package=spmixW)
[![License: GPL v3+](https://img.shields.io/badge/License-GPL%20v3%2B-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Bayesian Spatial Panel Data Models with Convex Combinations of Weight Matrices.**

`spmixW` provides Bayesian Markov chain Monte Carlo (MCMC) estimation of
spatial panel data models including Spatial Autoregressive (SAR), Spatial
Durbin Model (SDM), Spatial Error Model (SEM), Spatial Durbin Error Model
(SDEM), and Spatial Lag of X (SLX) specifications with fixed effects.
The package supports convex combinations of multiple spatial weight matrices
and Bayesian Model Averaging (BMA) over subsets of weight matrices,
following Debarsy and LeSage (2021).

A Python port is available at
[Mustapha-Wasseja/spmixw-py](https://github.com/Mustapha-Wasseja/spmixw-py)
(install with `pip install spmixw`).

## Installation

Install the released version from CRAN:

```r
install.packages("spmixW")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Mustapha-Wasseja/spmixW")
```

## Quick start

```r
library(spmixW)

# Create coordinates and a weight matrix
set.seed(42)
N <- 80
coords <- cbind(runif(N), runif(N))
W <- as.matrix(normw(make_knw(coords, k = 5, row_normalise = FALSE)))

# Simulate panel data from a SAR data-generating process
T <- 10
Wbig <- kronecker(diag(T), W)
X <- matrix(rnorm(N * T * 2), ncol = 2)
y <- as.numeric(
  solve(diag(N * T) - 0.5 * Wbig) %*% (X %*% c(1.5, -0.8) + rnorm(N * T))
)

# Estimate
res <- spmodel(
  y ~ x1 + x2,
  data = data.frame(
    region = rep(1:N, T), year = rep(1:T, each = N),
    y = y, x1 = X[, 1], x2 = X[, 2]
  ),
  W = W,
  model = "sar",
  id = "region", time = "year",
  effects = "twoway",
  ndraw = 5000, nomit = 1000
)

print(res)
```

## Models supported

| Model | Function | Spatial parameter | Effects decomposition |
|---|---|---|---|
| OLS | `ols_panel()` | none | β only |
| SAR | `sar_panel()` | ρ (lag) | LeSage-Pace direct/indirect/total |
| SDM | `sdm_panel()` | ρ (lag) | LeSage-Pace direct/indirect/total |
| SEM | `sem_panel()` | λ (error) | β only |
| SDEM | `sdem_panel()` | λ (error) | direct = β, indirect = θ |
| SLX | `slx_panel()` | none | direct = β, indirect = θ |

All models support fixed-effects specifications (`"none"`, `"region"`,
`"time"`, `"twoway"`) and optional heteroscedastic errors via the
Geweke (1993) Student-t mixture.

## Convex combination models

When the appropriate weight matrix is uncertain, `spmixW` can estimate a
**convex combination** of multiple weight matrices,
W_c = Σ γ_m W_m, with the γ weights estimated jointly with the
spatial parameter:

```r
W1 <- as.matrix(normw(make_knw(coords, k = 4)))
W2 <- as.matrix(normw(make_knw(coords, k = 8)))

res <- sar_conv_panel(
  y, X, list(W1, W2), N, T,
  ndraw = 20000, nomit = 5000,
  prior = list(model = 3, taylor_order = 6)
)
```

Available convex variants: `sar_conv_panel()`, `sdm_conv_panel()`,
`sem_conv_panel()`, `sdem_conv_panel()`.

## Bayesian Model Averaging

When several candidate weight matrices are available, BMA averages over
all 2^M − 1 non-empty subsets:

```r
res_bma <- sar_conv_bma(
  y, X, list(W1, W2, W3), N, T,
  ndraw = 15000, nomit = 5000
)
print(res_bma)  # shows model probabilities and BMA-averaged estimates
```

## Model comparison

`lmarginal_panel()` computes log-marginal likelihoods for SLX, SDM, and
SDEM specifications:

```r
dm <- demean_panel(y, X, N, T, model = 3)
lm_res <- lmarginal_panel(dm$ywith, dm$xwith, W, N, T)
model_probs(lm_res$lmarginal)
```

## License

GPL (>= 3)

## References

- Debarsy, N. and LeSage, J.P. (2021). "Bayesian model averaging for spatial
  autoregressive models based on convex combinations of different types of
  connectivity matrices." *Journal of Business & Economic Statistics*,
  40(2), 547-558. <doi:10.1080/07350015.2020.1840993>
- LeSage, J.P. and Pace, R.K. (2009). *Introduction to Spatial Econometrics*.
  Taylor & Francis/CRC Press.
- LeSage, J.P. (2020). "Fast MCMC estimation of multiple W-matrix spatial
  regression models and Metropolis-Hastings Monte Carlo log-marginal
  likelihoods." *Journal of Geographical Systems*, 22(1), 47-75.
- Geweke, J. (1993). "Bayesian Treatment of the Independent Student-t
  Linear Model." *Journal of Applied Econometrics*, 8(S1), S19-S40.
- Pace, R.K. and Barry, J.P. (1997). "Quick Computation of Spatial
  Autoregressive Estimators." *Geographical Analysis*, 29(3), 232-247.
