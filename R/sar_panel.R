#' Bayesian SAR Panel Model with Fixed Effects
#'
#' Estimates a Bayesian Spatial Autoregressive (SAR) panel model via MCMC:
#' \deqn{y = \rho W y + X \beta + \text{FE} + \varepsilon}
#' with griddy Gibbs sampling for \eqn{\rho} and LeSage-Pace scalar
#' effects estimates (direct, indirect, total).
#'
#' @param y Numeric vector of length \eqn{NT}.
#' @param X Numeric matrix \eqn{NT \times k}.
#' @param W Spatial weight matrix (\eqn{N \times N} or \eqn{NT \times NT}).
#' @param N Integer: number of cross-sectional units.
#' @param Time Integer: number of time periods.
#' @param ndraw Integer: total MCMC draws (including burn-in).
#' @param nomit Integer: burn-in draws to discard.
#' @param prior List of prior hyperparameters (see \code{\link{ols_panel}} for
#'   common fields). Additional fields:
#'   \describe{
#'     \item{lflag}{0 = exact log-det, 1 = MC approximation (default).}
#'     \item{order}{Taylor order for MC log-det (default 50).}
#'     \item{iter}{MC iterations for log-det (default 30).}
#'     \item{rmin, rmax}{Bounds for \eqn{\rho} grid (default \code{c(-1,1)}).}
#'     \item{lndet}{Pre-computed log-det grid (2-column matrix) to reuse.}
#'   }
#'
#' @return An S3 object of class \code{"spmixW"} with all fields from
#'   \code{\link{ols_panel}} plus:
#'   \describe{
#'     \item{rho}{Posterior mean of \eqn{\rho}.}
#'     \item{pdraw}{\code{coda::mcmc} of \eqn{\rho} draws.}
#'     \item{direct}{Matrix (p x 5): direct effects (mean, t-stat, p-value, lower, upper).}
#'     \item{indirect}{Matrix (p x 5): indirect effects.}
#'     \item{total}{Matrix (p x 5): total effects.}
#'     \item{direct_draws, indirect_draws, total_draws}{MCMC draws of effects.}
#'     \item{lndet}{Log-determinant grid (for reuse in subsequent calls).}
#'   }
#'
#' @details
#' The griddy Gibbs sampler for \eqn{\rho} evaluates the conditional
#' log-posterior on a fine grid using the concentrated likelihood
#' (beta integrated out), then uses inverse-CDF sampling via the
#' trapezoid rule.
#'
#' Effects are computed using the scalar summary measures of
#' LeSage and Pace (2009, Ch. 4) based on stochastic trace estimates
#' of powers of \eqn{W}.
#'
#' @references
#' LeSage, J.P. and Pace, R.K. (2009). \emph{Introduction to Spatial
#' Econometrics}. Taylor & Francis/CRC Press.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30; Time <- 10; rho_true <- 0.5
#' coords <- cbind(runif(N), runif(N))
#' W <- normw(make_knw(coords, k = 5, row_normalise = FALSE))
#' W <- as.matrix(W)
#' Wbig <- kronecker(diag(Time), W)
#' X <- matrix(rnorm(N * Time * 2), ncol = 2)
#' y <- solve(diag(N*Time) - rho_true * Wbig) %*% (X %*% c(1, -0.5) + rnorm(N*Time))
#' res <- sar_panel(as.numeric(y), X, W, N, Time, ndraw = 5000, nomit = 2000,
#'                  prior = list(model = 0, rval = 0))
#' print(res)
#' }
#'
#' @export
sar_panel <- function(y, X, W, N, Time, ndraw = 5500L, nomit = 1500L,
                      prior = list()) {

  start_time <- proc.time()

  # ---- Input validation ----
  y <- as.numeric(y)
  X <- as.matrix(X)
  nobs <- N * Time
  nvar <- ncol(X)

  stopifnot(length(y) == nobs, nrow(X) == nobs, ndraw > nomit)

  # ---- Handle W dimensions ----
  W <- Matrix::Matrix(W, sparse = TRUE)
  if (nrow(W) == N && ncol(W) == N) {
    Wbig <- kronecker(Matrix::Diagonal(Time), W)
  } else if (nrow(W) == nobs) {
    Wbig <- W
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
    "0" = "psar_g", "1" = "sarsfe_g", "2" = "sartfe_g", "3" = "sarstfe_g"
  )

  # ---- Compute Wy BEFORE demeaning ----
  Wy <- as.numeric(Wbig %*% y)

  # ---- Demean for fixed effects (including Wy) ----
  dm  <- demean_panel(y, X, N, Time, model)
  dmw <- demean_panel(Wy, matrix(0, nobs, 1), N, Time, model)
  ywith  <- dm$ywith
  xwith  <- dm$xwith
  wywith <- dmw$ywith

  # ---- Compute log-determinant grid ----
  # Use the N x N W (not the NT x NT Wbig)
  Wsmall <- if (nrow(W) == N) as.matrix(W) else as.matrix(Wbig[1:N, 1:N])

  if (!is.null(prior$lndet)) {
    detval <- prior$lndet
  } else if (lflag == 0L) {
    detval <- log_det_exact(Wsmall, rmin, rmax)
  } else {
    detval <- log_det_mc(Matrix::Matrix(Wsmall, sparse = TRUE),
                         rmin, rmax, order = mc_order, iter = mc_iter)
  }
  # Scale by T for panel: the SAR log-likelihood has T * log|I_N - rho*W|
  # We store the per-unit log-det and multiply by T inside draw_rho
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
  rho  <- 0.5
  vmean <- rep(0, nobs)
  save_idx <- 0L

  # ---- MCMC loop ----
  for (iter in seq_len(ndraw)) {

    if (homo) {
      # ---- Beta draw (homoscedastic) ----
      ys <- ywith - rho * wywith
      XtX <- crossprod(xwith)
      prec <- XtX + sige * Q
      AI   <- solve(prec)
      b0   <- AI %*% (crossprod(xwith, ys) + sige * Qpc)
      bhat <- b0 + t(chol(sige * AI)) %*% rnorm(nvar)

      # ---- Sige draw ----
      e <- ys - xwith %*% bhat
      nu1 <- nu + nobs
      d1  <- d0 + as.numeric(crossprod(e))
      sige <- d1 / rchisq(1, df = nu1)

      # ---- Rho draw (griddy Gibbs) ----
      # Concentrated residuals for griddy Gibbs
      xs <- xwith; ys_raw <- ywith; Wys_raw <- wywith
      b0g <- solve(crossprod(xs) + sige * Q, crossprod(xs, ys_raw) + sige * Qpc)
      bdg <- solve(crossprod(xs) + sige * Q, crossprod(xs, Wys_raw) + sige * Qpc)
      e0 <- ys_raw - xs %*% b0g
      ed <- Wys_raw - xs %*% bdg
      epe0  <- as.numeric(crossprod(e0))
      eped  <- as.numeric(crossprod(ed))
      epe0d <- as.numeric(crossprod(ed, e0))
      rho <- .draw_rho(detval_panel, epe0, eped, epe0d, nobs, nvar)

    } else {
      # ---- Beta draw (heteroscedastic) ----
      sqrtV <- sqrt(V)
      xs  <- xwith * sqrtV
      ys  <- ywith * sqrtV
      Wys <- wywith * sqrtV
      XtX <- crossprod(xs)
      prec <- XtX + sige * Q
      AI   <- solve(prec)
      yss  <- ys - rho * Wys
      b0   <- AI %*% (crossprod(xs, yss) + sige * Qpc)
      bhat <- b0 + t(chol(sige * AI)) %*% rnorm(nvar)

      # ---- Sige draw ----
      e <- yss - xs %*% bhat
      nu1 <- nu + nobs
      d1  <- d0 + as.numeric(crossprod(e))
      sige <- d1 / rchisq(1, df = nu1)

      # ---- Vi draw ----
      ev <- as.numeric(ywith - rho * wywith - xwith %*% bhat)
      chiv <- rchisq(nobs, df = rval + 1)
      vi <- as.numeric((ev^2 / sige + rval) / chiv)
      V  <- 1 / vi

      # ---- Rho draw (griddy Gibbs on V-scaled data) ----
      b0g <- solve(XtX + sige * Q, crossprod(xs, ys) + sige * Qpc)
      bdg <- solve(XtX + sige * Q, crossprod(xs, Wys) + sige * Qpc)
      e0 <- ys - xs %*% b0g
      ed <- Wys - xs %*% bdg
      epe0  <- as.numeric(crossprod(e0))
      eped  <- as.numeric(crossprod(ed))
      epe0d <- as.numeric(crossprod(ed, e0))
      rho <- .draw_rho(detval_panel, epe0, eped, epe0d, nobs, nvar)
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

  # ---- Detect number of non-intercept variables for effects ----
  # If first column of X is a constant, skip it for effects
  p <- nvar
  cflag <- FALSE
  if (all(abs(X[, 1] - X[1, 1]) < 1e-10)) {
    cflag <- TRUE
    p <- nvar - 1
  }

  # ---- Compute effects (LeSage-Pace scalar summary) ----
  effects <- .compute_sar_effects(bsave, psave, Wbig, nobs, nvar, p, cflag)

  # t-statistics
  bstd <- apply(bsave, 2, sd)
  pstd <- sd(psave)
  tstat <- c(bmean / bstd, rmean / pstd)

  # ---- Fixed effects recovery ----
  fe_result <- .compute_sar_fe(
    y, X, Wy, ywith, xwith, bmean, rmean, sige,
    dm$meanny, dm$meannx, dm$meanty, dm$meantx,
    dmw$meanny, dmw$meanty,
    N, Time, nobs, model
  )

  # R-squared
  resid <- y - rmean * Wy - fe_result$yhat_xb
  yme <- y - mean(y)
  rsqr <- 1 - sum(resid^2) / sum(yme^2)

  # Corr-squared
  ywhat <- xwith %*% bmean
  res1  <- ywith - mean(ywith)
  res2  <- ywhat - mean(ywith)
  corr2 <- (sum(res1 * res2))^2 / (sum(res1^2) * sum(res2^2))

  # Log-likelihood at posterior mean
  e_lik <- ywith - rmean * wywith - xwith %*% bmean
  lik <- -(nobs / 2) * log(2 * pi * sige) - sum(e_lik^2) / (2 * sige) +
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
    direct = effects$direct_out,
    indirect = effects$indirect_out,
    total  = effects$total_out,
    direct_draws   = effects$direct_draws,
    indirect_draws = effects$indirect_draws,
    total_draws    = effects$total_draws,
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
    p      = p,
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
# Internal: Griddy Gibbs draw for rho
# Griddy Gibbs for rho (LeSage & Pace 2009, Ch. 5)
# ============================================================
.draw_rho <- function(detval, epe0, eped, epe0d, nt, k) {
  # detval[,1] = rho grid, detval[,2] = log-det values (already *T for panel)
  nmk <- (nt - k) / 2
  nrho <- nrow(detval)
  rho_grid <- detval[, 1]
  lndet_vals <- detval[, 2]

  # Concentrated log-posterior on the grid
  z <- epe0 - 2 * rho_grid * epe0d + rho_grid^2 * eped
  z <- -nmk * log(z)
  den <- lndet_vals + z

  # Normalise (subtract max for numerical stability)
  adj <- max(den)
  den <- den - adj
  x <- exp(den)

  # Trapezoid rule for CDF
  n <- length(x)
  isum <- sum((rho_grid[2:n] + rho_grid[1:(n-1)]) * (x[2:n] - x[1:(n-1)]) / 2)
  z <- abs(x / isum)
  cdf <- cumsum(z)

  # Inverse CDF sampling
  rnd <- runif(1) * sum(z)
  ind <- which(cdf <= rnd)
  idraw <- if (length(ind) > 0) max(ind) else 1L
  if (idraw > 0 && idraw < nrho) {
    return(rho_grid[idraw])
  }
  return(rho_grid[max(1, idraw)])
}


# ============================================================
# Internal: SAR effects (LeSage-Pace scalar summary)
# ============================================================
.compute_sar_effects <- function(bsave, psave, Wbig, nobs, nvar, p, cflag) {

  # Stochastic trace estimation: tr(W^j) for j = 1, ..., maxorder
  uiter <- 50L
  maxorderu <- 100L
  rv <- matrix(rnorm(nobs * uiter), nrow = nobs, ncol = uiter)
  tracew <- numeric(maxorderu)
  wjjju <- rv
  for (j in seq_len(maxorderu)) {
    wjjju <- as.matrix(Wbig %*% wjjju)
    tracew[j] <- mean(colMeans(rv * wjjju))
  }

  # Override first two traces
  tracew[1] <- 0
  tracew[2] <- sum(Wbig * Wbig) / nobs  # exact tr(W^2)/N

  trs <- c(1, tracew)
  ntrs <- length(trs)
  ree <- 0:(ntrs - 1)

  ndrawsg <- nrow(bsave)
  # Variable indices (skip intercept if present)
  var_idx <- if (cflag) 2:nvar else 1:nvar
  var_idx <- var_idx[1:p]

  total_draws    <- matrix(0, ndrawsg, p)
  direct_draws   <- matrix(0, ndrawsg, p)
  indirect_draws <- matrix(0, ndrawsg, p)

  for (i in seq_len(ndrawsg)) {
    rmat <- psave[i]^ree
    for (j in seq_len(p)) {
      beta_j <- bsave[i, var_idx[j]]
      # Total effect for variable j: sum over traces of beta * rho^r
      total_val    <- sum(beta_j * rmat)
      # Direct effect: sum of (beta * tr(W^r)/N) * rho^r
      direct_val   <- sum(beta_j * trs * rmat)
      indirect_val <- total_val - direct_val

      total_draws[i, j]    <- total_val
      direct_draws[i, j]   <- direct_val
      indirect_draws[i, j] <- indirect_val
    }
  }

  # Summary statistics for each variable
  .effects_summary <- function(draws) {
    p_vars <- ncol(draws)
    out <- matrix(0, p_vars, 5)
    colnames(out) <- c("Mean", "t-stat", "p-value", "Lower05", "Upper95")
    for (j in seq_len(p_vars)) {
      m <- mean(draws[, j])
      s <- sd(draws[, j])
      t_val <- m / s
      p_val <- 2 * (1 - pt(abs(t_val), df = nrow(draws) - 1))
      bounds <- quantile(draws[, j], c(0.025, 0.975))
      out[j, ] <- c(m, t_val, p_val, bounds[1], bounds[2])
    }
    out
  }

  list(
    direct_out      = .effects_summary(direct_draws),
    indirect_out    = .effects_summary(indirect_draws),
    total_out       = .effects_summary(total_draws),
    direct_draws    = direct_draws,
    indirect_draws  = indirect_draws,
    total_draws     = total_draws
  )
}


# ============================================================
# Internal: SAR fixed effects recovery
# ============================================================
.compute_sar_fe <- function(y, X, Wy, ywith, xwith, beta, rho, sige,
                            meanny, meannx, meanty, meantx,
                            meannwy, meantwy,
                            N, Time, nobs, model) {
  intercept <- sfe <- tfe <- tsfe <- ttfe <- NULL

  XtXi <- tryCatch(solve(crossprod(xwith)), error = function(e) {
    MASS::ginv(as.matrix(crossprod(xwith)))
  })

  if (model == 0L) {
    yhat_xb <- as.numeric(X %*% beta)
    yhat <- yhat_xb
  } else if (model == 1L) {
    intercept <- mean(y) - mean(Wy) * rho - as.numeric(colMeans(X) %*% beta)
    sfe <- meanny - meannwy * rho - meannx %*% beta - intercept
    yhat_xb <- as.numeric(X %*% beta) + rep(as.numeric(sfe), times = Time) + intercept
    yhat <- yhat_xb
    tsfe <- sfe / sqrt(sige / Time + diag(sige * meannx %*% XtXi %*% t(meannx)))
  } else if (model == 2L) {
    intercept <- mean(y) - mean(Wy) * rho - as.numeric(colMeans(X) %*% beta)
    tfe <- meanty - meantwy * rho - meantx %*% beta - intercept
    yhat_xb <- as.numeric(X %*% beta) + rep(as.numeric(tfe), each = N) + intercept
    yhat <- yhat_xb
    ttfe <- tfe / sqrt(sige / N + diag(sige * meantx %*% XtXi %*% t(meantx)))
  } else {
    intercept <- mean(y) - mean(Wy) * rho - as.numeric(colMeans(X) %*% beta)
    sfe <- meanny - meannwy * rho - meannx %*% beta - intercept
    tfe <- meanty - meantwy * rho - meantx %*% beta - intercept
    yhat_xb <- as.numeric(X %*% beta) +
      rep(as.numeric(sfe), times = Time) +
      rep(as.numeric(tfe), each = N) + intercept
    yhat <- yhat_xb
    tsfe <- sfe / sqrt(sige / Time + diag(sige * meannx %*% XtXi %*% t(meannx)))
    ttfe <- tfe / sqrt(sige / N + diag(sige * meantx %*% XtXi %*% t(meantx)))
  }

  list(yhat = as.numeric(yhat), yhat_xb = as.numeric(yhat_xb),
       intercept = intercept, sfe = as.numeric(sfe), tfe = as.numeric(tfe),
       tsfe = as.numeric(tsfe), ttfe = as.numeric(ttfe))
}
