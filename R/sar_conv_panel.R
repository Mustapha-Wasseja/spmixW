#' Bayesian SAR Panel with Convex Combination of Weight Matrices
#'
#' Estimates a SAR panel model where the spatial weight matrix is a convex
#' combination \eqn{W_c = \sum_{m=1}^{M} \gamma_m W_m} with \eqn{\gamma}
#' on the unit simplex. Uses Metropolis-Hastings for both \eqn{\rho} and
#' \eqn{\gamma}, with Taylor series log-determinant approximation.
#'
#' @param y Numeric vector of length NT.
#' @param X Numeric matrix NT x k.
#' @param Wlist A list of M spatial weight matrices (each N x N).
#' @param N Integer: number of cross-sectional units.
#' @param Time Integer: number of time periods.
#' @param ndraw Integer: total MCMC draws (recommend >= 10000).
#' @param nomit Integer: burn-in draws.
#' @param prior List of prior hyperparameters. See Details.
#'
#' @return An S3 object of class \code{"spmixW"} with:
#'   \describe{
#'     \item{beta}{Posterior mean of beta.}
#'     \item{rho}{Posterior mean of rho.}
#'     \item{gamma}{Posterior mean of gamma (M x 1).}
#'     \item{bdraw, pdraw, sdraw, gdraw}{MCMC draws for beta, rho, sigma, gamma.}
#'     \item{rho_acc_rate}{Acceptance rate for rho MH sampler.}
#'     \item{gamma_acc_rate}{Acceptance rate for gamma MH sampler.}
#'     \item{direct, indirect, total}{Effects estimates.}
#'     \item{traces}{Pre-computed Taylor traces (for reuse).}
#'   }
#'
#' @details
#' The model is:
#' \deqn{y = \rho W_c(\gamma) y + X \beta + \text{FE} + \varepsilon, \quad
#' \varepsilon \sim N(0, \sigma^2 I_{NT})}
#'
#' The MCMC sampler draws:
#' \enumerate{
#'   \item \eqn{\beta | \text{rest}} from multivariate normal (conjugate).
#'   \item \eqn{\sigma^2 | \text{rest}} from inverse-gamma.
#'   \item \eqn{\rho | \text{rest}} via Metropolis-Hastings random walk
#'     with adaptive tuning (target acceptance 40-60\%).
#'   \item \eqn{\gamma | \text{rest}} via Metropolis-Hastings with reversible
#'     jump proposal on the simplex (matching LeSage's coin-flip method).
#' }
#'
#' Pre-computes \eqn{\tilde{y} = [y, -W_1 y, \ldots, -W_M y]} so that
#' for any \eqn{(\rho, \gamma)}: \eqn{(I - \rho W_c) y = \tilde{y} \omega}
#' where \eqn{\omega = (1, \rho \gamma_1, \ldots, \rho \gamma_M)'}.
#'
#' @section Taylor approximation accuracy:
#' The convex combination models use a Taylor series approximation for the
#' log-determinant. The default order is 6 (extending the original order 4
#' of Debarsy and LeSage 2021). Approximation accuracy improves with larger N: at N=500 and
#' rho=0.6, expect approximately 0.13 upward bias in rho; at N=3000+ the
#' bias is negligible (<0.02). Users can control this via
#' \code{prior$taylor_order}. For small-N applications where precise rho
#' estimation is critical, compare results against \code{\link{sar_panel}}
#' with a known single W matrix as a benchmark.
#'
#' @references
#' Debarsy, N. and LeSage, J.P. (2021). "Bayesian model averaging for spatial
#' autoregressive models based on convex combinations of different types of
#' connectivity matrices." \emph{Journal of Business & Economic Statistics},
#' 40(2), 547-558.
#'
#' LeSage, J.P. (2020). "Fast MCMC estimation of multiple W-matrix spatial
#' regression models and Metropolis-Hastings Monte Carlo log-marginal
#' likelihoods." \emph{Journal of Geographical Systems}, 22(1), 47-75.
#'
#' @export
sar_conv_panel <- function(y, X, Wlist, N, Time, ndraw = 25000L,
                           nomit = 5000L, prior = list()) {

  start_time <- proc.time()

  y <- as.numeric(y)
  X <- as.matrix(X)
  nobs <- N * Time
  nvar <- ncol(X)
  M <- length(Wlist)

  stopifnot(length(y) == nobs, nrow(X) == nobs, ndraw > nomit, M >= 1)

  # ---- Parse priors ----
  model  <- prior$model %||% 0L
  thin   <- prior$thin  %||% 1L
  rmin   <- prior$rmin  %||% -0.9999
  rmax   <- prior$rmax  %||% 0.9999
  taylor_order <- prior$taylor_order %||% 6L

  # ---- Build big W matrices ----
  Wbigs <- lapply(Wlist, function(Wm) {
    Wm <- Matrix::Matrix(Wm, sparse = TRUE)
    if (nrow(Wm) == N) kronecker(Matrix::Diagonal(Time), Wm) else Wm
  })

  # ---- Demean ----
  dm <- demean_panel(y, X, N, Time, model)
  ywith <- dm$ywith
  xwith <- dm$xwith

  # ---- Pre-compute Taylor traces ----
  # Use N x N matrices for traces (more efficient)
  Wsmalls <- lapply(Wlist, function(Wm) {
    Wm <- as.matrix(Wm)
    if (nrow(Wm) == N) Wm else Wm[1:N, 1:N]
  })

  if (!is.null(prior$traces)) {
    traces <- prior$traces
  } else {
    traces <- log_det_taylor(Wsmalls, max_order = taylor_order)
  }
  # Scale for panel: multiply all trace terms by T
  # Scale ALL trace terms by T for panel log-det: T * log|I_N - rho*W|
  traces_panel <- traces
  for (p_name in names(traces$traces)) {
    traces_panel$traces[[p_name]] <- Time * traces$traces[[p_name]]
  }
  # Also scale backward-compat fields
  if (!is.null(traces_panel$T2)) traces_panel$T2 <- Time * traces$T2
  if (!is.null(traces_panel$T3)) traces_panel$T3 <- Time * traces$T3
  if (!is.null(traces_panel$T4)) traces_panel$T4 <- Time * traces$T4

  # ---- Pre-compute ys: [ywith, -W1*ywith, ..., -WM*ywith] ----
  ys <- matrix(0, nobs, M + 1)
  ys[, 1] <- ywith
  for (m in seq_len(M)) {
    ys[, m + 1] <- -as.numeric(Wbigs[[m]] %*% ywith)
  }

  # ---- Pre-compute xs: array [nobs, M+1, k] ----
  xs <- array(0, dim = c(nobs, M + 1, nvar))
  for (j in seq_len(nvar)) {
    xs[, 1, j] <- xwith[, j]
    for (m in seq_len(M)) {
      xs[, m + 1, j] <- -as.numeric(Wbigs[[m]] %*% xwith[, j])
    }
  }

  # ---- Detect intercept for effects ----
  cflag <- all(abs(X[, 1] - X[1, 1]) < 1e-10)
  p <- if (cflag) nvar - 1 else nvar

  # ---- Initialise MCMC ----
  n_save <- floor((ndraw - nomit) / thin)
  bsave <- matrix(0, n_save, nvar)
  ssave <- numeric(n_save)
  psave <- numeric(n_save)
  gsave <- matrix(0, n_save, M)

  rho   <- 0.5
  sige  <- 1
  gamma <- rep(1 / M, M)
  save_idx <- 0L

  # MH tuning
  cc <- 0.1  # rho proposal sd
  rho_acc <- 0L
  gamma_acc <- 0L
  gam_std <- rep(0.1, M)
  dd <- 3.0
  noo <- 1000L  # adaptation interval

  # ---- Helper: evaluate log-conditional for (rho, gamma) ----
  .eval_cond <- function(rho_val, gamma_val, bhat_val) {
    omega <- c(1, rho_val * gamma_val)
    # Reconstruct filtered data
    xo <- matrix(0, nobs, nvar)
    for (j in seq_len(nvar)) xo[, j] <- xs[, , j] %*% omega
    e <- ys %*% omega - xo %*% bhat_val
    epe <- as.numeric(crossprod(e))
    nmk <- (nobs - nvar) / 2
    xpx <- crossprod(xo)
    detx <- 0.5 * log(det(xpx))
    # Taylor log-det (already scaled by T)
    detm <- eval_taylor_lndet(traces_panel, rho_val, gamma_val)
    detm - detx - nmk * log(epe)
  }

  # ---- MCMC loop ----
  for (iter in seq_len(ndraw)) {

    # ---- 1. Draw beta ----
    omega <- c(1, rho * gamma)
    xo <- matrix(0, nobs, nvar)
    for (j in seq_len(nvar)) xo[, j] <- xs[, , j] %*% omega
    AI <- solve(crossprod(xo))
    b0 <- AI %*% crossprod(xo, ys %*% omega)
    bhat <- b0 + t(chol(sige * AI)) %*% rnorm(nvar)

    # ---- 2. Draw sigma^2 ----
    e <- ys %*% omega - xo %*% bhat
    d1 <- as.numeric(crossprod(e))
    sige <- d1 / rchisq(1, df = nobs)

    # ---- 3. Draw rho (MH random walk) ----
    rho2 <- rho + cc * rnorm(1)
    # Reflect into bounds
    if (rho2 > rmin && rho2 < rmax) {
      alpMH <- .eval_cond(rho2, gamma, bhat) - .eval_cond(rho, gamma, bhat)
      if (log(runif(1)) < alpMH) {
        rho <- rho2
        rho_acc <- rho_acc + 1L
      }
    }

    # Adapt cc
    acc_rate_rho <- rho_acc / iter
    if (acc_rate_rho < 0.4) cc <- cc / 1.1
    if (acc_rate_rho > 0.6) cc <- cc * 1.1
    if (cc < 0.001) cc <- 0.1
    if (cc > 1.0) cc <- 0.1

    # ---- 4. Draw gamma (MH reversible jump on simplex) ----
    # Skip gamma sampling if M = 1 (gamma is always 1)
    gnew <- NULL
    if (M >= 2) {
      if (iter <= noo) {
        gnew <- .gamma_proposal_uniform(gamma, M)
      } else {
        gnew <- .gamma_proposal_adapted(gamma, M, gam_std, dd)
      }
    }

    if (!is.null(gnew)) {
      alpMH <- .eval_cond(rho, gnew, bhat) - .eval_cond(rho, gamma, bhat)
      if (log(runif(1)) < alpMH) {
        gamma <- gnew
        gamma_acc <- gamma_acc + 1L
      }
    }

    # Adapt dd
    acc_rate_gamma <- gamma_acc / iter
    if (acc_rate_gamma < 0.10) dd <- dd / 1.1
    if (acc_rate_gamma > 0.40) dd <- dd * 1.1
    if (dd > 3.0) dd <- 3.0
    if (dd < 1.0) dd <- 1.0

    # Update gam_std periodically
    if (iter == noo || (iter > noo && iter %% noo == 0L)) {
      start_i <- max(1L, save_idx - noo + 1L)
      if (save_idx > 10) {
        gam_std <- pmax(apply(gsave[start_i:save_idx, , drop = FALSE], 2, sd), 0.01)
      }
    }

    # ---- Save draws ----
    if (iter > nomit && ((iter - nomit) %% thin == 0L)) {
      save_idx <- save_idx + 1L
      bsave[save_idx, ] <- as.numeric(bhat)
      ssave[save_idx]   <- sige
      psave[save_idx]   <- rho
      gsave[save_idx, ] <- gamma
    }
  }

  # ---- Post-processing ----
  bmean <- colMeans(bsave)
  rmean <- mean(psave)
  smean <- mean(ssave)
  gmean <- colMeans(gsave)

  # ---- Log-marginal likelihood (Chib-style) ----
  # Evaluate log-posterior at each saved draw
  drawpost <- numeric(n_save)
  for (i in seq_len(n_save)) {
    drawpost[i] <- .eval_cond(psave[i], gsave[i, ], bsave[i, ])
  }

  # Normalising constant at posterior mean
  omega_post <- c(1, rmean * gmean)
  xo_post <- matrix(0, nobs, nvar)
  for (j in seq_len(nvar)) xo_post[, j] <- xs[, , j] %*% omega_post
  xpx_post <- crossprod(xo_post)
  lndetx <- log(det(xpx_post))
  dof <- (nobs - nvar) / 2
  D <- 1 - 1 / rmin
  logC <- -log(D) + lgamma(dof) - dof * log(2 * pi) - 0.5 * lndetx

  logmarginal <- mean(drawpost) + logC

  # Build Wc at posterior mean gamma for effects
  Wc <- Matrix::Matrix(0, nobs, nobs, sparse = TRUE)
  for (m in seq_len(M)) Wc <- Wc + gmean[m] * Wbigs[[m]]

  # Effects (using posterior mean Wc)
  effects <- .compute_sar_effects(bsave, psave, Wc, nobs, nvar, p, cflag)

  # t-statistics
  bstd <- apply(bsave, 2, sd)
  pstd <- sd(psave)
  tstat <- c(bmean / bstd, rmean / pstd)

  # Fixed effects recovery
  Wy_orig <- as.numeric(Wc %*% y)
  dmw <- demean_panel(Wy_orig, matrix(0, nobs, 1), N, Time, model)

  fe_result <- .compute_sar_fe(
    y, X, Wy_orig, ywith, xwith, bmean, rmean, smean,
    dm$meanny, dm$meannx, dm$meanty, dm$meantx,
    dmw$meanny, dmw$meanty,
    N, Time, nobs, model
  )

  resid <- y - rmean * Wy_orig - fe_result$yhat_xb
  yme <- y - mean(y)
  rsqr <- 1 - sum(resid^2) / sum(yme^2)

  meth <- switch(as.character(model),
    "0" = "psar_conv_g", "1" = "sarsfe_conv_g",
    "2" = "sartfe_conv_g", "3" = "sarstfe_conv_g"
  )

  result <- list(
    beta   = bmean,
    rho    = rmean,
    sige   = smean,
    gamma  = gmean,
    bdraw  = coda::mcmc(bsave, start = nomit + 1, thin = thin),
    pdraw  = coda::mcmc(matrix(psave, ncol = 1), start = nomit + 1, thin = thin),
    sdraw  = coda::mcmc(matrix(ssave, ncol = 1), start = nomit + 1, thin = thin),
    gdraw  = coda::mcmc(gsave, start = nomit + 1, thin = thin),
    rho_acc_rate   = rho_acc / ndraw,
    gamma_acc_rate = gamma_acc / ndraw,
    direct   = effects$direct_out,
    indirect = effects$indirect_out,
    total    = effects$total_out,
    direct_draws   = effects$direct_draws,
    indirect_draws = effects$indirect_draws,
    total_draws    = effects$total_draws,
    logmarginal = logmarginal,
    traces = traces,
    yhat   = fe_result$yhat,
    resid  = resid,
    rsqr   = rsqr,
    intercept = fe_result$intercept,
    sfe    = fe_result$sfe,
    tfe    = fe_result$tfe,
    nobs   = nobs,
    nvar   = nvar,
    p      = p,
    N      = N,
    Time   = Time,
    M      = M,
    model  = model,
    tstat  = tstat,
    meth   = meth,
    ndraw  = ndraw,
    nomit  = nomit,
    time   = (proc.time() - start_time)[3]
  )

  result <- .add_field_aliases(result)
  class(result) <- "spmixW"
  result
}


# ---- Internal: gamma proposal (uniform, early iterations) ----
.gamma_proposal_uniform <- function(gamma, M) {
  gtst <- numeric(M - 1)
  for (j in seq_len(M - 1)) {
    coin <- runif(1)
    if (coin <= 1/3) {
      gtst[j] <- runif(1, 0, gamma[j])
    } else if (coin <= 2/3) {
      gtst[j] <- gamma[j]
    } else {
      gtst[j] <- runif(1, gamma[j], 1)
    }
  }
  gnew <- c(gtst, 1 - sum(gtst))

  # Reject if any negative (retry up to 50 times)
  attempts <- 0L
  while (any(gnew < 0) && attempts < 50L) {
    for (j in seq_len(M - 1)) {
      coin <- runif(1)
      if (coin <= 1/3) {
        gtst[j] <- runif(1, 0, gamma[j])
      } else if (coin <= 2/3) {
        gtst[j] <- gamma[j]
      } else {
        gtst[j] <- runif(1, gamma[j], 1)
      }
    }
    gnew <- c(gtst, 1 - sum(gtst))
    attempts <- attempts + 1L
  }

  if (any(gnew < 0)) return(NULL)
  gnew
}


# ---- Internal: gamma proposal (adapted, later iterations) ----
.gamma_proposal_adapted <- function(gamma, M, gam_std, dd) {
  gtst <- numeric(M - 1)
  for (j in seq_len(M - 1)) {
    coin <- runif(1)
    if (coin <= 1/3) {
      gtst[j] <- runif(1, gamma[j] - dd * gam_std[j], gamma[j])
    } else if (coin <= 2/3) {
      gtst[j] <- gamma[j]
    } else {
      gtst[j] <- runif(1, gamma[j], gamma[j] + dd * gam_std[j])
    }
  }
  gnew <- c(gtst, 1 - sum(gtst))

  attempts <- 0L
  while (any(gnew < 0) && attempts < 50L) {
    for (j in seq_len(M - 1)) {
      coin <- runif(1, 0.0001, 0.9999)
      if (coin <= 1/3) {
        gtst[j] <- runif(1, gamma[j] - dd * gam_std[j], gamma[j])
      } else if (coin <= 2/3) {
        gtst[j] <- gamma[j]
      } else {
        gtst[j] <- runif(1, gamma[j], gamma[j] + dd * gam_std[j])
      }
    }
    gnew <- c(gtst, 1 - sum(gtst))
    attempts <- attempts + 1L
  }

  if (any(gnew < 0)) return(NULL)
  gnew
}
