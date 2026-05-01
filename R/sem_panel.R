#' Bayesian SEM Panel Model with Fixed Effects
#'
#' Estimates a Bayesian Spatial Error Model (SEM) panel via MCMC:
#' \deqn{y = X \beta + \text{FE} + u, \quad u = \lambda W u + \varepsilon}
#' with griddy Gibbs sampling for \eqn{\lambda}.
#'
#' @inheritParams sar_panel
#'
#' @return An S3 object of class \code{"spmixW"} with fields similar to
#'   \code{\link{ols_panel}} plus:
#'   \describe{
#'     \item{rho}{Posterior mean of \eqn{\lambda} (spatial error parameter).}
#'     \item{pdraw}{\code{coda::mcmc} of \eqn{\lambda} draws.}
#'     \item{lndet}{Log-determinant grid.}
#'   }
#'
#' @details
#' The SEM differs from SAR in that the spatial parameter enters the error
#' structure. The MCMC sampler filters both y and X by \eqn{(I - \lambda W)}
#' before drawing \eqn{\beta}. The griddy Gibbs for \eqn{\lambda} evaluates
#' the concentrated log-posterior on a grid, where for each grid point the
#' filtered data \eqn{(I - \lambda W) y} and \eqn{(I - \lambda W) X} are
#' used to compute the residual sum of squares.
#'
#' @references
#' LeSage, J.P. and Pace, R.K. (2009). \emph{Introduction to Spatial
#' Econometrics}. Taylor & Francis/CRC Press.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30; Time <- 10; lambda_true <- 0.5
#' coords <- cbind(runif(N), runif(N))
#' W <- as.matrix(normw(make_knw(coords, k = 5, row_normalise = FALSE)))
#' Wbig <- kronecker(diag(Time), W)
#' X <- matrix(rnorm(N * Time * 2), ncol = 2)
#' u <- solve(diag(N*Time) - lambda_true * Wbig) %*% rnorm(N*Time)
#' y <- X %*% c(1, -0.5) + u
#' res <- sem_panel(as.numeric(y), X, W, N, Time, ndraw = 5000, nomit = 2000,
#'                  prior = list(model = 0, rval = 0))
#' print(res)
#' }
#'
#' @export
sem_panel <- function(y, X, W, N, Time, ndraw = 5500L, nomit = 1500L,
                      prior = list()) {

  start_time <- proc.time()

  # ---- Input validation ----
  y <- as.numeric(y)
  X <- as.matrix(X)
  nobs <- N * Time
  nvar <- ncol(X)

  stopifnot(length(y) == nobs, nrow(X) == nobs, ndraw > nomit)

  # ---- Handle W dimensions ----
  W_sp <- Matrix::Matrix(W, sparse = TRUE)
  if (nrow(W_sp) == N && ncol(W_sp) == N) {
    Wbig <- kronecker(Matrix::Diagonal(Time), W_sp)
  } else if (nrow(W_sp) == nobs) {
    Wbig <- W_sp
  } else {
    stop("W must be N x N or NT x NT")
  }

  # ---- Parse priors ----
  model   <- prior$model %||% 0L
  rval    <- prior$rval  %||% 4
  nu      <- prior$nu    %||% 0
  d0      <- prior$d0    %||% 0
  thin    <- prior$thin  %||% 1L
  lflag   <- prior$lflag %||% 1L
  mc_order <- prior$order %||% 50L
  mc_iter  <- prior$iter  %||% 30L
  rmin    <- prior$rmin  %||% -1
  rmax    <- prior$rmax  %||% 1

  c_beta  <- prior$beta_mean %||% rep(0, nvar)
  C_beta  <- prior$beta_var  %||% (diag(nvar) * 1e12)

  stopifnot(model %in% 0:3)

  Q   <- solve(C_beta)
  Qpc <- Q %*% c_beta
  homo <- (rval == 0)

  meth <- switch(as.character(model),
    "0" = "psem_g", "1" = "semsfe_g", "2" = "semtfe_g", "3" = "semstfe_g"
  )

  # ---- Demean for fixed effects ----
  dm <- demean_panel(y, X, N, Time, model)
  ywith <- dm$ywith
  xwith <- dm$xwith

  # Pre-compute W*ywith, W*xwith (used in each SEM iteration)
  Wy <- as.numeric(Wbig %*% ywith)
  Wx <- as.matrix(Wbig %*% xwith)

  # ---- Compute log-determinant grid ----
  Wsmall <- if (nrow(W_sp) == N) as.matrix(W_sp) else as.matrix(Wbig[1:N, 1:N])

  if (!is.null(prior$lndet)) {
    detval <- prior$lndet
  } else if (lflag == 0L) {
    detval <- log_det_exact(Wsmall, rmin, rmax)
  } else {
    detval <- log_det_mc(Matrix::Matrix(Wsmall, sparse = TRUE),
                         rmin, rmax, order = mc_order, iter = mc_iter)
  }
  # Panel scaling: T * log|I_N - lambda*W|
  detval_panel <- detval
  detval_panel[, "lndet"] <- Time * detval[, "lndet"]

  # ---- Initialise MCMC ----
  n_save <- floor((ndraw - nomit) / thin)
  bsave <- matrix(0, nrow = n_save, ncol = nvar)
  ssave <- numeric(n_save)
  psave <- numeric(n_save)

  V    <- rep(1, nobs)
  vi   <- rep(1, nobs)
  sige <- 1
  rho  <- 0.5  # lambda (spatial error parameter, stored as rho for consistency)
  vmean <- rep(0, nobs)
  save_idx <- 0L

  # ---- MCMC loop ----
  for (iter in seq_len(ndraw)) {

    if (homo) {
      # ---- Beta draw: filter by (I - lambda*W) ----
      xs <- xwith - rho * Wx
      ys <- ywith - rho * Wy
      XtX <- crossprod(xs)
      prec <- XtX + sige * Q
      AI   <- solve(prec)
      b0   <- AI %*% (crossprod(xs, ys) + sige * Qpc)
      bhat <- b0 + t(chol(sige * AI)) %*% rnorm(nvar)

      # ---- Sige draw ----
      e  <- ywith - xwith %*% bhat
      es <- e - rho * as.numeric(Wbig %*% e)
      nu1 <- nobs + 2 * nu
      d1  <- 2 * d0 + as.numeric(crossprod(es))
      sige <- d1 / rchisq(1, df = nu1)

      # ---- Lambda draw (griddy Gibbs) ----
      rho <- .draw_rho_sem(detval_panel, ywith, xwith, Wy, Wx,
                           V, nobs, nvar, rmin, rmax)

    } else {
      # ---- Heteroscedastic beta draw ----
      sqrtV <- sqrt(V)
      xs_raw <- xwith * sqrtV
      ys_raw <- ywith * sqrtV
      Wxs <- as.matrix(Wbig %*% xs_raw)
      Wys <- as.numeric(Wbig %*% ys_raw)
      xss <- xs_raw - rho * Wxs
      yss <- ys_raw - rho * Wys
      XtX <- crossprod(xss)
      prec <- XtX + sige * Q
      AI   <- solve(prec)
      b0   <- AI %*% (crossprod(xss, yss) + sige * Qpc)
      bhat <- b0 + t(chol(sige * AI)) %*% rnorm(nvar)

      # ---- Sige draw ----
      e  <- ys_raw - xs_raw %*% bhat
      ed <- e - rho * as.numeric(Wbig %*% e)
      nu1 <- nobs + 2 * nu
      d1  <- 2 * d0 + as.numeric(crossprod(ed))
      sige <- d1 / rchisq(1, df = nu1)

      # ---- Vi draw ----
      ev <- as.numeric(ywith * sqrtV - (xwith * sqrtV) %*% bhat)
      chiv <- rchisq(nobs, df = rval + 1)
      vi <- as.numeric((ev^2 / sige + rval) / chiv)
      V  <- 1 / vi

      # ---- Lambda draw ----
      rho <- .draw_rho_sem(detval_panel, ywith, xwith, Wy, Wx,
                           V, nobs, nvar, rmin, rmax)
    }

    # ---- Save draws ----
    if (iter > nomit && ((iter - nomit) %% thin == 0L)) {
      save_idx <- save_idx + 1L
      bsave[save_idx, ] <- as.numeric(bhat)
      ssave[save_idx]   <- sige
      psave[save_idx]   <- rho
      vmean <- vmean + vi
    }
  }

  vmean <- vmean / n_save

  # ---- Post-processing ----
  bmean <- colMeans(bsave)
  smean <- mean(ssave)
  rmean <- mean(psave)
  sige  <- smean

  bstd <- apply(bsave, 2, sd)
  pstd <- sd(psave)
  tstat <- c(bmean / bstd, rmean / pstd)

  # Fixed effects (SEM: no spatial lag term in FE recovery)
  fe_result <- .compute_ols_fe(
    y, X, ywith, xwith, bmean, sige,
    dm$meanny, dm$meannx, dm$meanty, dm$meantx,
    N, Time, nobs, model
  )

  # Fit statistics
  resid <- y - fe_result$yhat
  yme <- y - mean(y)
  rsqr <- 1 - sum(resid^2) / sum(yme^2)

  ywhat <- xwith %*% bmean
  res1  <- ywith - mean(ywith)
  res2  <- ywhat - mean(ywith)
  corr2 <- (sum(res1 * res2))^2 / (sum(res1^2) * sum(res2^2))

  # Log-likelihood at posterior mean
  e_lik <- ywith - xwith %*% bmean
  es_lik <- e_lik - rmean * as.numeric(Wbig %*% e_lik)
  lik <- -(nobs / 2) * log(2 * pi * sige) - sum(es_lik^2) / (2 * sige) +
    Time * approx(detval[, "rho"], detval[, "lndet"], xout = rmean)$y

  # ---- Assemble output ----
  result <- list(
    beta   = bmean,
    rho    = rmean,
    sige   = smean,
    bdraw  = coda::mcmc(bsave, start = nomit + 1, thin = thin),
    pdraw  = coda::mcmc(matrix(psave, ncol = 1), start = nomit + 1, thin = thin),
    sdraw  = coda::mcmc(matrix(ssave, ncol = 1), start = nomit + 1, thin = thin),
    vmean  = vmean,
    tstat  = tstat,
    yhat   = fe_result$yhat,
    resid  = resid,
    rsqr   = rsqr,
    corr2  = corr2,
    lik    = lik,
    intercept = fe_result$intercept,
    sfe    = fe_result$sfe,
    tfe    = fe_result$tfe,
    tsfe   = fe_result$tsfe,
    ttfe   = fe_result$ttfe,
    lndet  = detval,
    nobs   = nobs,
    nvar   = nvar,
    N      = N,
    Time   = Time,
    model  = model,
    meth   = meth,
    ndraw  = ndraw,
    nomit  = nomit,
    rval   = rval,
    cov    = cov(cbind(bsave, psave)),
    time   = (proc.time() - start_time)[3]
  )

  result <- .add_field_aliases(result)
  class(result) <- "spmixW"
  result
}


# ============================================================
# Internal: Griddy Gibbs for lambda in SEM
# Different from SAR: must recompute filtered residuals at each grid point
# Griddy Gibbs for lambda in SEM (LeSage & Pace 2009, Ch. 5)
# ============================================================
.draw_rho_sem <- function(detval, y, x, Wy, Wx, V, n, k, rmin, rmax) {

  nmk <- (n - k) / 2
  nrho <- nrow(detval)
  rho_grid <- detval[, 1]
  lndet_vals <- detval[, 2]

  # For SEM, we need epe(rho) = e(rho)'e(rho) where
  # e(rho) = (I - rho*W)*y - (I - rho*W)*X * b_hat(rho)
  # and b_hat(rho) = ((I-rho*W)X)'((I-rho*W)X))^{-1} ((I-rho*W)X)'((I-rho*W)y)
  #
  # Compute on a coarse grid, then interpolate to fine grid

  sqrtV <- sqrt(V)
  coarse_grid <- seq(rmin + 0.01, rmax - 0.01, by = 0.01)
  ng <- length(coarse_grid)
  epet <- numeric(ng)

  for (i in seq_len(ng)) {
    rho_i <- coarse_grid[i]
    xs <- (x - rho_i * Wx) * sqrtV
    ys <- (y - rho_i * Wy) * sqrtV
    bs <- solve(crossprod(xs), crossprod(xs, ys))
    e <- ys - xs %*% bs
    epet[i] <- as.numeric(crossprod(e))
  }

  # Interpolate to the fine rho grid
  epe <- approx(coarse_grid, epet, xout = rho_grid, rule = 2)$y

  # Concentrated log-posterior
  den <- lndet_vals - nmk * log(epe)

  # Normalise and sample
  adj <- max(den)
  den <- den - adj
  x_den <- exp(den)

  nr <- length(x_den)
  isum <- sum((rho_grid[2:nr] + rho_grid[1:(nr-1)]) *
              (x_den[2:nr] - x_den[1:(nr-1)]) / 2)
  z <- abs(x_den / isum)
  cdf <- cumsum(z)

  rnd <- runif(1) * sum(z)
  ind <- which(cdf <= rnd)
  idraw <- if (length(ind) > 0) max(ind) else 1L
  if (idraw > 0 && idraw < nrho) {
    return(rho_grid[idraw])
  }
  return(rho_grid[max(1, idraw)])
}
