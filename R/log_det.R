#' Exact Log-Determinant via Eigenvalues
#'
#' Computes \eqn{\log|I_N - \rho W|} for a grid of \eqn{\rho} values
#' using the eigenvalues of \eqn{W}. This is exact but requires computing
#' all eigenvalues of \eqn{W}, which is \eqn{O(N^3)}.
#'
#' @param W A square (sparse or dense) spatial weight matrix (\eqn{N \times N}).
#' @param rmin Numeric scalar: lower bound of \eqn{\rho} grid (default \code{-1}).
#' @param rmax Numeric scalar: upper bound of \eqn{\rho} grid (default \code{1}).
#' @param grid_step Numeric scalar: grid spacing (default \code{0.001}).
#'
#' @return A matrix with two columns:
#'   \describe{
#'     \item{rho}{Grid of \eqn{\rho} values.}
#'     \item{lndet}{Corresponding \eqn{\log|I_N - \rho W|} values.}
#'   }
#'
#' @details
#' Given eigenvalues \eqn{\lambda_1, \ldots, \lambda_N} of \eqn{W}:
#' \deqn{\log|I_N - \rho W| = \sum_{i=1}^{N} \log(1 - \rho \lambda_i)}
#' The result is pre-computed on a fine grid and used for griddy Gibbs sampling
#' of the spatial autoregressive parameter.
#'
#' @references
#' LeSage, J.P. and Pace, R.K. (2009). \emph{Introduction to Spatial
#' Econometrics}. Taylor & Francis/CRC Press.
#'
#' @examples
#' \donttest{
#' W <- matrix(c(0, 0.5, 0.5, 0,
#'               0.5, 0, 0, 0.5,
#'               0.5, 0, 0, 0.5,
#'               0, 0.5, 0.5, 0), 4, 4)
#' detval <- log_det_exact(W)
#' head(detval)
#' }
#'
#' @export
log_det_exact <- function(W, rmin = -1, rmax = 1, grid_step = 0.001) {

  W <- as.matrix(W)
  N <- nrow(W)
  stopifnot(N == ncol(W))

  # Compute eigenvalues of W (real parts only for asymmetric W)
  eig_vals <- Re(eigen(W, only.values = TRUE)$values)

  rho_grid <- seq(rmin, rmax, by = grid_step)
  # For each rho, sum log(1 - rho * lambda_i) over all eigenvalues
  # Vectorised: outer product approach
  # log_terms is length(rho_grid) x N
  log_terms <- outer(rho_grid, eig_vals, function(r, lam) log(1 - r * lam))
  lndet <- rowSums(log_terms)

  cbind(rho = rho_grid, lndet = lndet)
}


#' Monte Carlo Log-Determinant Approximation (Pace-Barry 1999)
#'
#' Approximates \eqn{\log|I_N - \rho W|} using the stochastic trace estimator
#' of Barry and Pace (1999). Suitable for large spatial weight matrices where
#' exact eigenvalue computation is infeasible.
#'
#' @param W A square sparse spatial weight matrix (\eqn{N \times N}).
#' @param rmin Numeric scalar: lower bound of \eqn{\rho} grid (default \code{-1}).
#' @param rmax Numeric scalar: upper bound of \eqn{\rho} grid (default \code{1}).
#' @param grid_step Numeric scalar: grid spacing (default \code{0.001}).
#' @param order Integer: number of terms in the Taylor expansion (default \code{50}).
#' @param iter Integer: number of random vectors for trace estimation (default \code{30}).
#'
#' @return A matrix with two columns:
#'   \describe{
#'     \item{rho}{Grid of \eqn{\rho} values.}
#'     \item{lndet}{Corresponding approximate \eqn{\log|I_N - \rho W|} values.}
#'   }
#'
#' @details
#' Uses the identity:
#' \deqn{\log|I - \rho W| = -\sum_{j=1}^{\infty} \frac{\rho^j}{j} \text{tr}(W^j)}
#' The traces \eqn{\text{tr}(W^j)} are estimated stochastically:
#' \deqn{\text{tr}(W^j) \approx \frac{1}{q} \sum_{l=1}^{q} u_l' W^j u_l}
#' where \eqn{u_l} are random \eqn{N(0, I)} vectors.
#'
#' The second-order trace \eqn{\text{tr}(W^2)} is computed exactly as
#' \eqn{\text{tr}(W' W) = \sum_{ij} w_{ij}^2} for greater accuracy.
#'
#' @references
#' Pace, R.K. and Barry, J.P. (1997). "Quick Computation of Spatial
#' Autoregressive Estimators." \emph{Geographical Analysis}, 29(3), 232-247.
#'
#' Barry, R.P. and Pace, R.K. (1999). "Monte Carlo estimates of the log
#' determinant of large sparse matrices." \emph{Linear Algebra and its
#' Applications}, 289, 41-54.
#'
#' @examples
#' \donttest{
#' library(Matrix)
#' N <- 100
#' W <- sparseMatrix(
#'   i = c(1:N, 1:(N-1)),
#'   j = c(c(2:N, 1), 2:N),
#'   x = rep(0.5, 2*N - 1),
#'   dims = c(N, N)
#' )
#' detval <- log_det_mc(W)
#' head(detval)
#' }
#'
#' @export
log_det_mc <- function(W, rmin = -1, rmax = 1, grid_step = 0.001,
                       order = 50L, iter = 30L) {

  W <- Matrix::Matrix(W, sparse = TRUE)
  N <- nrow(W)

  # Estimate tr(W^j) for j = 1, ..., order using random vectors
  traces <- numeric(order)

  # Generate random vectors
  U <- matrix(rnorm(N * iter), nrow = N, ncol = iter)

  # Iteratively compute W^j * U and estimate traces
  WjU <- U  # W^0 * U = U
  for (j in seq_len(order)) {
    WjU <- as.matrix(W %*% WjU)  # W^j * U
    # tr(W^j) ≈ mean(u' * W^j * u) = mean(colSums(U * WjU))
    traces[j] <- mean(colSums(U * WjU))
  }

  # Override trace(W) = 0 for row-standardised W (diagonal is 0)
  traces[1] <- 0
  # Exact computation of tr(W^2) = tr(W'W) = sum of squared elements
  traces[2] <- sum(W * W)  # element-wise product, then sum

  # Compute log-det on grid: log|I - rho*W| = -sum_{j=1}^{order} (rho^j / j) * tr(W^j)
  rho_grid <- seq(rmin, rmax, by = grid_step)
  j_seq <- seq_len(order)

  # Outer product: rho_grid^j / j, then multiply by traces
  # rho_powers is length(rho_grid) x order
  rho_powers <- outer(rho_grid, j_seq, "^")
  coeffs <- traces / j_seq  # order-length vector

  lndet <- -as.numeric(rho_powers %*% coeffs)

  cbind(rho = rho_grid, lndet = lndet)
}
