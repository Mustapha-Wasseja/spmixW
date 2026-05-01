#' Log-Marginal Likelihoods for Static Spatial Panel Models
#'
#' Computes log-marginal likelihoods for three spatial panel specifications
#' (SLX, SDM, SDEM) under diffuse priors on \eqn{\beta, \sigma^2} and a
#' uniform prior on \eqn{\rho / \lambda}. Used for Bayesian model comparison.
#'
#' @param y Numeric vector of length NT (should be demeaned if FE are used).
#' @param X Numeric matrix NT x k (WITHOUT intercept for FE models).
#' @param W Spatial weight matrix (N x N).
#' @param N Integer: number of cross-sectional units.
#' @param Time Integer: number of time periods.
#' @param prior Optional list with fields:
#'   \describe{
#'     \item{lflag}{0 = exact log-det, 1 = MC approximation (default).}
#'     \item{order}{MC order (default 50).}
#'     \item{iter}{MC iterations (default 30).}
#'     \item{rmin, rmax}{Bounds for rho grid.}
#'   }
#'
#' @return A list with:
#'   \describe{
#'     \item{logm_slx}{Log-marginal for SLX model.}
#'     \item{logm_sdm}{Log-marginal for SDM model (integrated over rho).}
#'     \item{logm_sdem}{Log-marginal for SDEM model (integrated over lambda).}
#'     \item{lmarginal}{Vector of all three log-marginals.}
#'     \item{probs}{Posterior model probabilities (assuming equal priors).}
#'   }
#'
#' @details
#' The SLX model has no spatial parameter, so its marginal is analytic.
#' The SDM and SDEM marginals require numerical integration over
#' \eqn{\rho} (or \eqn{\lambda}) using the trapezoid rule on a fine grid.
#'
#' For SDM, the concentrated log-marginal profile is:
#' \deqn{-(dof) \log(e_0'e_0 - 2\rho e_0'e_d + \rho^2 e_d'e_d) + T \log|I - \rho W|}
#' where \eqn{e_0, e_d} are OLS residuals from y and Wy on the SDM design matrix.
#'
#' @references
#' LeSage, J.P. (2014). "What Regional Scientists Need to Know about Spatial
#' Econometrics." \emph{Review of Regional Studies}, 44(1), 13-32.
#'
#' LeSage, J.P. and Pace, R.K. (2009). \emph{Introduction to Spatial
#' Econometrics}. Taylor & Francis/CRC Press.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30; Time <- 10
#' coords <- cbind(runif(N), runif(N))
#' W <- as.matrix(normw(make_knw(coords, k = 5, row_normalise = FALSE)))
#' X <- matrix(rnorm(N * Time * 2), ncol = 2)
#' y <- rnorm(N * Time)
#' res <- lmarginal_panel(y, X, W, N, Time)
#' res$probs
#' }
#'
#' @export
lmarginal_panel <- function(y, X, W, N, Time, prior = list()) {

  y <- as.numeric(y)
  X <- as.matrix(X)
  nobs <- N * Time
  nx <- ncol(X)

  # Parse options
  lflag   <- prior$lflag %||% 1L
  mc_order <- prior$order %||% 50L
  mc_iter  <- prior$iter  %||% 30L
  rmin    <- prior$rmin  %||% -0.9999
  rmax    <- prior$rmax  %||% 0.9999

  # Build large W
  W_sp <- Matrix::Matrix(W, sparse = TRUE)
  Wbig <- kronecker(Matrix::Diagonal(Time), W_sp)

  # SDM design matrix: [1, X, WX]
  WX <- as.matrix(Wbig %*% X)
  xsdm <- cbind(rep(1, nobs), X, WX)

  # ---- Degrees of freedom ----
  dof <- (nobs - 1) / 2

  # Uniform prior width
  D <- 1 - 1 / rmin

  # ---- SLX model (no spatial parameter) ----
  xpx_sdm <- crossprod(xsdm)
  lndetx_sdm <- log(det(xpx_sdm))

  logC_slx <- lgamma(dof) - dof * log(2 * pi) - 0.5 * lndetx_sdm

  bo <- solve(xpx_sdm, crossprod(xsdm, y))
  eo <- y - xsdm %*% bo
  epeo <- as.numeric(crossprod(eo))

  logm_slx <- -dof * log(epeo) + logC_slx

  # ---- SDM model (integrate over rho) ----
  logC_sdm <- -log(D) + lgamma(dof) - dof * log(2 * pi) - 0.5 * lndetx_sdm

  Wy <- as.numeric(Wbig %*% y)
  bd <- solve(xpx_sdm, crossprod(xsdm, Wy))
  ed <- Wy - xsdm %*% bd
  eped <- as.numeric(crossprod(ed))
  epeod <- as.numeric(crossprod(ed, eo))

  # Rho grid
  incr <- 0.001
  xx <- seq(rmin, rmax, by = incr)
  ngrid <- length(xx)

  # Log-det on this grid
  Wsmall <- as.matrix(W_sp)
  if (lflag == 0L) {
    lndet_vals <- log_det_exact(Wsmall, rmin, rmax, grid_step = incr)
  } else {
    lndet_vals <- log_det_mc(Matrix::Matrix(Wsmall, sparse = TRUE),
                             rmin, rmax, grid_step = incr,
                             order = mc_order, iter = mc_iter)
  }
  # Match grid lengths
  lndet_at_xx <- approx(lndet_vals[, "rho"], lndet_vals[, "lndet"],
                         xout = xx, rule = 2)$y

  # SDM profile
  z_sdm <- epeo - 2 * xx * epeod + xx^2 * eped
  logm_sdm_profile <- -dof * log(z_sdm) + Time * lndet_at_xx

  adj_sdm <- max(logm_sdm_profile)
  yy_sdm <- exp(logm_sdm_profile - adj_sdm)

  # Trapezoid integration
  isum_sdm <- sum(diff(xx) * (yy_sdm[-1] + yy_sdm[-ngrid]) / 2)
  logm_sdm <- log(isum_sdm) + adj_sdm + logC_sdm

  # ---- SDEM model (integrate over lambda) ----
  logC_sdem <- -log(D) + lgamma(dof) - dof * log(2 * pi)

  # For SDEM, need to compute at each grid point:
  # A(lambda) = (I - lambda*W), then filtered cross-products
  xpx <- crossprod(xsdm)
  xpWx <- crossprod(xsdm, as.matrix(Wbig %*% xsdm))
  Wxsdm <- as.matrix(Wbig %*% xsdm)
  xpWpx <- crossprod(Wxsdm, xsdm)
  xpWpWx <- crossprod(Wxsdm)

  xpy <- crossprod(xsdm, y)
  xpWy <- crossprod(xsdm, Wy)
  xpWpy <- crossprod(Wxsdm, y)
  xpWpWy <- crossprod(Wxsdm, Wy)

  ypy <- as.numeric(crossprod(y))
  ypWy <- as.numeric(crossprod(y, Wy))
  ypWpWy <- as.numeric(crossprod(Wy))

  Q1 <- numeric(ngrid)
  Q3 <- numeric(ngrid)

  for (i in seq_len(ngrid)) {
    rho_i <- xx[i]

    Axx <- xpx - rho_i * xpWx - rho_i * xpWpx + rho_i^2 * xpWpWx
    Axy <- xpy - rho_i * xpWy - rho_i * xpWpy + rho_i^2 * xpWpWy
    Ayy <- ypy - rho_i * ypWy - rho_i * ypWy + rho_i^2 * ypWpWy

    Q3[i] <- log(det(Axx))
    b_i <- solve(Axx, Axy)
    Q1[i] <- Ayy - as.numeric(crossprod(b_i, Axx %*% b_i))
  }

  logm_sdem_profile <- -dof * log(Q1) + Time * lndet_at_xx - 0.5 * Q3

  adj_sdem <- max(logm_sdem_profile)
  yy_sdem <- exp(logm_sdem_profile - adj_sdem)

  isum_sdem <- sum(diff(xx) * (yy_sdem[-1] + yy_sdem[-ngrid]) / 2)
  logm_sdem <- log(isum_sdem) + adj_sdem + logC_sdem

  # ---- Model probabilities ----
  lmarginal <- c(logm_slx, logm_sdm, logm_sdem)
  names(lmarginal) <- c("SLX", "SDM", "SDEM")
  probs <- model_probs(lmarginal)
  names(probs) <- names(lmarginal)

  list(
    logm_slx  = logm_slx,
    logm_sdm  = logm_sdm,
    logm_sdem = logm_sdem,
    lmarginal = lmarginal,
    probs     = probs
  )
}
