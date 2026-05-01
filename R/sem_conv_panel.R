#' Bayesian SEM Panel with Convex Combination of Weight Matrices
#'
#' Estimates a SEM panel model with \eqn{W_c = \sum \gamma_m W_m}:
#' \deqn{y = X\beta + u, \quad u = \lambda W_c u + \varepsilon}
#'
#' @inheritParams sar_conv_panel
#'
#' @return An S3 object of class \code{"spmixW"} with gamma, rho (lambda),
#'   acceptance rates, and MCMC draws.
#'
#' @details
#' Follows the same MH structure as \code{\link{sar_conv_panel}} but the
#' spatial parameter enters the error structure. The pre-computed arrays
#' \code{ys} and \code{xs} represent \eqn{(I - \lambda W_c)} filtering.
#'
#' @section Taylor approximation accuracy:
#' See \code{\link{sar_conv_panel}} for details on Taylor order and accuracy.
#'
#' @references
#' Debarsy, N. and LeSage, J.P. (2021). "Bayesian model averaging for spatial
#' autoregressive models based on convex combinations of different types of
#' connectivity matrices." \emph{Journal of Business & Economic Statistics},
#' 40(2), 547-558.
#'
#' @export
sem_conv_panel <- function(y, X, Wlist, N, Time, ndraw = 25000L,
                           nomit = 5000L, prior = list()) {

  start_time <- proc.time()

  y <- as.numeric(y)
  X <- as.matrix(X)
  nobs <- N * Time
  nvar <- ncol(X)
  M <- length(Wlist)

  stopifnot(length(y) == nobs, nrow(X) == nobs, ndraw > nomit, M >= 1)

  model  <- prior$model %||% 0L
  thin   <- prior$thin  %||% 1L
  taylor_order <- prior$taylor_order %||% 6L
  rmin   <- prior$rmin  %||% -0.9999
  rmax   <- prior$rmax  %||% 0.9999

  # Build big W matrices
  Wbigs <- lapply(Wlist, function(Wm) {
    Wm <- Matrix::Matrix(Wm, sparse = TRUE)
    if (nrow(Wm) == N) kronecker(Matrix::Diagonal(Time), Wm) else Wm
  })

  # Demean
  dm <- demean_panel(y, X, N, Time, model)
  ywith <- dm$ywith
  xwith <- dm$xwith

  # Taylor traces
  Wsmalls <- lapply(Wlist, function(Wm) {
    Wm <- as.matrix(Wm)
    if (nrow(Wm) == N) Wm else Wm[1:N, 1:N]
  })
  if (!is.null(prior$traces)) {
    traces <- prior$traces
  } else {
    traces <- log_det_taylor(Wsmalls, max_order = taylor_order)
  }
  # Scale ALL trace terms by T for panel log-det
  traces_panel <- traces
  for (p_name in names(traces$traces)) {
    traces_panel$traces[[p_name]] <- Time * traces$traces[[p_name]]
  }
  if (!is.null(traces_panel$T2)) traces_panel$T2 <- Time * traces$T2
  if (!is.null(traces_panel$T3)) traces_panel$T3 <- Time * traces$T3
  if (!is.null(traces_panel$T4)) traces_panel$T4 <- Time * traces$T4

  # Pre-compute ys and xs for SEM: omega = (1, lambda*gamma)
  # (I - lambda*Wc)*y = y - lambda*Wc*y = ys * omega where ys = [y, -W1*y, ..., -WM*y]
  ys <- matrix(0, nobs, M + 1)
  ys[, 1] <- ywith
  for (m in seq_len(M)) ys[, m + 1] <- -as.numeric(Wbigs[[m]] %*% ywith)

  xs <- array(0, dim = c(nobs, M + 1, nvar))
  for (j in seq_len(nvar)) {
    xs[, 1, j] <- xwith[, j]
    for (m in seq_len(M)) xs[, m + 1, j] <- -as.numeric(Wbigs[[m]] %*% xwith[, j])
  }

  # Initialise
  n_save <- floor((ndraw - nomit) / thin)
  bsave <- matrix(0, n_save, nvar)
  ssave <- numeric(n_save)
  psave <- numeric(n_save)
  gsave <- matrix(0, n_save, M)

  lam   <- 0.5
  sige  <- 1
  gamma <- rep(1 / M, M)
  save_idx <- 0L

  cc <- 0.1; lam_acc <- 0L; gamma_acc <- 0L
  gam_std <- rep(0.1, M); dd <- 3.0; noo <- 1000L

  .eval_cond <- function(lam_val, gamma_val, bhat_val) {
    omega <- c(1, lam_val * gamma_val)
    xo <- matrix(0, nobs, nvar)
    for (j in seq_len(nvar)) xo[, j] <- xs[, , j] %*% omega
    e <- ys %*% omega - xo %*% bhat_val
    epe <- as.numeric(crossprod(e))
    nmk <- (nobs - nvar) / 2
    xpx <- crossprod(xo)
    detx <- 0.5 * log(det(xpx))
    detm <- eval_taylor_lndet(traces_panel, lam_val, gamma_val)
    detm - detx - nmk * log(epe)
  }

  for (iter in seq_len(ndraw)) {
    # 1. Draw beta (SEM: filter by (I - lam*Wc))
    omega <- c(1, lam * gamma)
    xo <- matrix(0, nobs, nvar)
    for (j in seq_len(nvar)) xo[, j] <- xs[, , j] %*% omega
    AI <- solve(crossprod(xo))
    b0 <- AI %*% crossprod(xo, ys %*% omega)
    bhat <- b0 + t(chol(sige * AI)) %*% rnorm(nvar)

    # 2. Draw sigma^2
    e <- ys %*% omega - xo %*% bhat
    d1 <- as.numeric(crossprod(e))
    sige <- d1 / rchisq(1, df = nobs)

    # 3. Draw lambda (MH)
    lam2 <- lam + cc * rnorm(1)
    if (lam2 > rmin && lam2 < rmax) {
      alpMH <- .eval_cond(lam2, gamma, bhat) - .eval_cond(lam, gamma, bhat)
      if (log(runif(1)) < alpMH) { lam <- lam2; lam_acc <- lam_acc + 1L }
    }
    acc_rate_lam <- lam_acc / iter
    if (acc_rate_lam < 0.4) cc <- cc / 1.1
    if (acc_rate_lam > 0.6) cc <- cc * 1.1
    if (cc < 0.001 || cc > 1.0) cc <- 0.1

    # 4. Draw gamma (MH) — skip if M = 1
    gnew <- NULL
    if (M >= 2) {
      gnew <- if (iter <= noo) .gamma_proposal_uniform(gamma, M) else
        .gamma_proposal_adapted(gamma, M, gam_std, dd)
    }
    if (!is.null(gnew)) {
      alpMH <- .eval_cond(lam, gnew, bhat) - .eval_cond(lam, gamma, bhat)
      if (log(runif(1)) < alpMH) { gamma <- gnew; gamma_acc <- gamma_acc + 1L }
    }
    acc_rate_g <- gamma_acc / iter
    if (acc_rate_g < 0.10) dd <- dd / 1.1
    if (acc_rate_g > 0.40) dd <- dd * 1.1
    dd <- max(1.0, min(3.0, dd))

    if (iter == noo || (iter > noo && iter %% noo == 0L)) {
      start_i <- max(1L, save_idx - noo + 1L)
      if (save_idx > 10) gam_std <- pmax(apply(gsave[start_i:save_idx, , drop = FALSE], 2, sd), 0.01)
    }

    if (iter > nomit && ((iter - nomit) %% thin == 0L)) {
      save_idx <- save_idx + 1L
      bsave[save_idx, ] <- as.numeric(bhat)
      ssave[save_idx]   <- sige
      psave[save_idx]   <- lam
      gsave[save_idx, ] <- gamma
    }
  }

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

  bstd <- apply(bsave, 2, sd)
  pstd <- sd(psave)
  tstat <- c(bmean / bstd, rmean / pstd)

  fe_result <- .compute_ols_fe(
    y, X, ywith, xwith, bmean, smean,
    dm$meanny, dm$meannx, dm$meanty, dm$meantx,
    N, Time, nobs, model
  )
  resid <- y - fe_result$yhat
  yme <- y - mean(y)
  rsqr <- 1 - sum(resid^2) / sum(yme^2)

  meth <- switch(as.character(model),
    "0" = "psem_conv_g", "1" = "semsfe_conv_g",
    "2" = "semtfe_conv_g", "3" = "semstfe_conv_g"
  )

  result <- list(
    beta = bmean, rho = rmean, sige = smean, gamma = gmean,
    bdraw = coda::mcmc(bsave, start = nomit + 1, thin = thin),
    pdraw = coda::mcmc(matrix(psave, ncol = 1), start = nomit + 1, thin = thin),
    sdraw = coda::mcmc(matrix(ssave, ncol = 1), start = nomit + 1, thin = thin),
    gdraw = coda::mcmc(gsave, start = nomit + 1, thin = thin),
    rho_acc_rate = lam_acc / ndraw,
    gamma_acc_rate = gamma_acc / ndraw,
    logmarginal = logmarginal,
    traces = traces,
    yhat = fe_result$yhat, resid = resid, rsqr = rsqr,
    intercept = fe_result$intercept, sfe = fe_result$sfe, tfe = fe_result$tfe,
    tstat = tstat, nobs = nobs, nvar = nvar, N = N, Time = Time, M = M,
    model = model, meth = meth, ndraw = ndraw, nomit = nomit,
    time = (proc.time() - start_time)[3]
  )
  result <- .add_field_aliases(result)
  class(result) <- "spmixW"
  result
}
