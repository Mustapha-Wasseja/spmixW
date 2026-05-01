#' spmixW: Bayesian Spatial Panel Data Models with Convex Combinations
#' of Weight Matrices
#'
#' R implementation of Bayesian MCMC spatial panel data models, inspired by
#' the MATLAB Spatial Econometrics Toolbox. Supports SAR, SDM, SEM, SDEM,
#' and SLX specifications with fixed effects, convex combinations of multiple
#' spatial weight matrices, and Bayesian Model Averaging.
#'
#' @references
#' Debarsy, N. and LeSage, J.P. (2021). "Bayesian model averaging for spatial
#' autoregressive models based on convex combinations of different types of
#' connectivity matrices." \emph{Journal of Business & Economic Statistics},
#' 40(2), 547-558.
#'
#' LeSage, J.P. and Pace, R.K. (2009). \emph{Introduction to Spatial
#' Econometrics}. Taylor & Francis/CRC Press.
#'
#' LeSage, J.P. (2020). "Fast MCMC estimation of multiple W-matrix spatial
#' regression models and Metropolis-Hastings Monte Carlo log-marginal
#' likelihoods." \emph{Journal of Geographical Systems}, 22(1), 47-75.
#'
#' Geweke, J. (1993). "Bayesian Treatment of the Independent Student-t Linear
#' Model." \emph{Journal of Applied Econometrics}, 8(S1), S19-S40.
#'
#' Pace, R.K. and Barry, J.P. (1997). "Quick Computation of Spatial
#' Autoregressive Estimators." \emph{Geographical Analysis}, 29(3), 232-247.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Matrix Matrix Diagonal sparseMatrix solve t rowSums colSums summary
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rchisq runif cov sd pt quantile var dist approx density
#'   model.frame terms na.fail setNames
#' @importFrom coda mcmc effectiveSize geweke.diag
#' @importFrom graphics abline barplot legend lines par plot.new
#' @importFrom utils combn
## usethis namespace: end
NULL
