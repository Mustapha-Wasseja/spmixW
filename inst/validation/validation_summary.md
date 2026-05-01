# spmixW Package Validation Summary

## Overview

Validation of the `spmixW` R package against the MATLAB `panelg` toolbox
by James P. LeSage, based on the demo programs in `panelg_manual.pdf`.

**Date:** 2026-03-31
**R version:** 4.5.2
**Package version:** 0.1.0

## Results

| # | Validation | Tests | Passed | Failed | Status |
|---|-----------|-------|--------|--------|--------|
| 1 | OLS Panel (homo, hetero, prior) | 9 | 9 | 0 | PASS |
| 2 | SAR Panel (rho, beta, effects) | 13 | 13 | 0 | PASS |
| 3 | SDM Panel (beta, theta, rho, effects) | 9 | 9 | 0 | PASS |
| 4 | SEM Panel (beta, lambda) | 3 | 3 | 0 | PASS |
| 5 | SDEM Panel (beta, theta, lambda, effects) | 9 | 9 | 0 | PASS |
| 6 | SLX Panel (beta, theta, effects) | 6 | 6 | 0 | PASS |
| 7 | SAR Convex Combination (rho, gamma, acc rates) | 10 | 10 | 0 | PASS |
| 8 | Model Comparison (log-marginals) | 3 | 3 | 0 | PASS |
| 9 | SAR Convex BMA (model probs, BMA beta) | 5 | 5 | 0 | PASS |
| 10 | Cross-Sectional T=1 (SAR + SAR conv) | 8 | 8 | 0 | PASS |
| **TOTAL** | | **75** | **75** | **0** | **ALL PASS** |

## Key Findings

### Standard Models (Validations 1-6)

- **OLS**: Beta estimates within 0.07 of truth. Heteroscedastic vi correctly
  inflates for outlier periods (ratio 1.70x). Tight prior pulls estimates
  toward prior mean (beta ~ 0.63 vs truth 1.0 and prior 0.5).
- **SAR**: rho = 0.589 (true 0.6), matching MATLAB's 0.599 within MCMC noise.
  Direct + indirect = total exact to machine precision. Spatial multiplier
  correctly inflates direct effects above raw beta.
- **SDM**: All four coefficients (beta, theta) recovered within 0.03 of truth.
  Indirect effects correctly negative (theta = -1 dominates).
- **SEM**: Lambda = 0.676 (true 0.7), beta within 0.023.
- **SDEM**: All parameters within tolerance. Direct effects ~ beta,
  indirect ~ theta (SDEM-specific decomposition confirmed).
- **SLX**: Exact effects decomposition confirmed (direct = beta, indirect = theta).

### Convex Combination Models (Validations 7-10)

- **SAR Conv (M=3)**: Correctly identifies gamma ordering (gamma2 > gamma1 > gamma3).
  Gamma3 (true 0.0) estimated at 0.045 (near zero). MH acceptance rates
  healthy: rho 47%, gamma 41%.
- **SAR Conv (M=2)**: Tighter estimates with fewer matrices. Gamma1 = 0.263
  (true 0.3), gamma2 = 0.737 (true 0.7).
- **Model Comparison**: SDM correctly identified with probability 1.0000 when
  data was generated from SDM DGP.
- **BMA**: W1-only model gets probability 1.0000 when DGP uses only W1.
  BMA-averaged beta within 0.06 of truth.
- **Cross-Sectional (T=1)**: All functions work correctly with single time period.
  SAR conv with T=1 recovers gamma ordering.

### Comparison with MATLAB Output

Where MATLAB printed output was available from the manual:

| Parameter | True | R Estimate | MATLAB Estimate | Agreement |
|-----------|------|-----------|----------------|-----------|
| OLS beta_1 (run a) | 1.0 | 1.002 | 0.965 | Both within tolerance |
| OLS sigma^2 (run a) | 5.0 | 4.411 | 4.353 | Close agreement |
| OLS beta_1 (prior) | 0.5* | 0.640 | 0.629 | Excellent agreement |
| SAR rho | 0.6 | 0.589 | 0.599 | Both within 0.02 of truth |
| SAR sigma^2 | 1.0 | 0.867 | 0.869 | Excellent agreement |
| SDEM beta_1 | 1.0 | 0.999 | 0.989 | Excellent agreement |
| SDEM lambda | 0.7 | 0.677 | 0.688 | Both within 0.025 |

*Note: Different RNG streams between R and MATLAB produce different
random data, so exact output match is not expected. The key validation
is that both recover truth within similar tolerances.

## Failures

**None.** All 75 checks passed.

## Notes

- Validation uses N=200 (scaled down from MATLAB's N=3111 for computational
  feasibility). Parameter recovery is slightly less precise with smaller N
  but still within reasonable MCMC tolerance.
- MCMC draws: standard models use 5000 draws / 1500 burn-in; convex models
  use 15000-20000 draws / 5000 burn-in; BMA uses 10000 / 3000 per model.
- R and MATLAB use different RNG streams (even with same seed), so data
  differs. Validation focuses on parameter recovery from truth, not
  exact replication of MATLAB output.
