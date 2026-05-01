#' Bayesian OLS Panel Model with Fixed Effects
#'
#' Estimates a Bayesian panel regression model (no spatial component) with
#' optional spatial and/or time-period fixed effects via MCMC sampling.
#' Supports both homoscedastic and heteroscedastic error structures.
#'
#' @param y Numeric vector of length \eqn{NT}: dependent variable.
#'   Data must be sorted by time then by spatial unit: all N regions in
#'   period 1, then all N in period 2, etc.
#' @param X Numeric matrix of dimension \eqn{NT \times k}: explanatory variables.
#' @param N Integer: number of cross-sectional units.
#' @param Time Integer: number of time periods.
#' @param ndraw Integer: total number of MCMC draws (including burn-in).
#' @param nomit Integer: number of initial draws to discard as burn-in.
#' @param prior A list of prior hyperparameters:
#'   \describe{
#'     \item{model}{Integer 0-3: fixed-effects specification
#'       (0 = pooled, 1 = region FE, 2 = time FE, 3 = both). Default 0.}
#'     \item{rval}{Numeric > 0: degrees-of-freedom for heteroscedastic errors.
#'       Each \eqn{v_i \sim \chi^2(r)/r}. Set to 0 for homoscedastic. Default 4.}
#'     \item{beta_mean}{Numeric vector of length k: prior mean for \eqn{\beta}.
#'       Default \code{rep(0, k)}.}
#'     \item{beta_var}{Numeric \eqn{k \times k} matrix: prior variance for
#'       \eqn{\beta}. Default \code{diag(k) * 1e12} (diffuse).}
#'     \item{nu}{Numeric >= 0: shape parameter for inverse-gamma prior on
#'       \eqn{\sigma^2}. Default 0 (diffuse).}
#'     \item{d0}{Numeric >= 0: scale parameter for inverse-gamma prior on
#'       \eqn{\sigma^2}. Default 0 (diffuse).}
#'     \item{thin}{Integer >= 1: thinning interval. Default 1 (no thinning).}
#'   }
#'
#' @return An S3 object of class \code{"spmixW"} containing:
#'   \describe{
#'     \item{beta}{Posterior mean of \eqn{\beta} (length k).}
#'     \item{sige}{Posterior mean of \eqn{\sigma^2}.}
#'     \item{bdraw}{\code{coda::mcmc} object: \eqn{(ndraw - nomit) \times k}
#'       matrix of retained \eqn{\beta} draws.}
#'     \item{sdraw}{\code{coda::mcmc} object: retained \eqn{\sigma^2} draws.}
#'     \item{vmean}{Posterior mean of the variance weights \eqn{v_i}
#'       (length NT; all ones if homoscedastic).}
#'     \item{yhat}{Fitted values including fixed effects.}
#'     \item{resid}{Residuals \eqn{y - \hat{y}}.}
#'     \item{rsqr}{R-squared.}
#'     \item{corr2}{Squared correlation between actual and fitted demeaned values.}
#'     \item{sfe}{Spatial fixed effects (if \code{model \%in\% c(1,3)}).}
#'     \item{tfe}{Time fixed effects (if \code{model \%in\% c(2,3)}).}
#'     \item{intercept}{Estimated intercept.}
#'     \item{tstat}{Posterior t-statistics (mean / sd of draws).}
#'     \item{nobs, nvar, N, Time, model, meth}{Metadata.}
#'     \item{time}{Wall-clock time in seconds.}
#'   }
#'
#' @details
#' The model is:
#' \deqn{y = X\beta + \iota_T \otimes \mu + \nu \otimes \iota_N + \varepsilon}
#' where \eqn{\varepsilon \sim N(0, \sigma^2 V)}, \eqn{V = \text{diag}(v_1, \ldots, v_{NT})}.
#'
#' \strong{MCMC sampling steps:}
#' \enumerate{
#'   \item Draw \eqn{\beta | \sigma^2, V, y} from the multivariate normal
#'     posterior (conjugate update).
#'   \item Draw \eqn{\sigma^2 | \beta, V, y} from the inverse-gamma posterior.
#'   \item (Heteroscedastic only) Draw each \eqn{v_i | \beta, \sigma^2, y}
#'     from \eqn{(\epsilon_i^2 / \sigma^2 + r) / \chi^2(r+1)}.
#' }
#'
#' @references
#' LeSage, J.P. and Pace, R.K. (2009). \emph{Introduction to Spatial
#' Econometrics}. Taylor & Francis/CRC Press.
#'
#' Geweke, J. (1993). "Bayesian Treatment of the Independent Student-t Linear
#' Model." \emph{Journal of Applied Econometrics}, 8(S1), S19-S40.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' N <- 20; Time <- 10; k <- 2
#' beta_true <- c(1.5, -0.8)
#' X <- matrix(rnorm(N * Time * k), ncol = k)
#' y <- X %*% beta_true + rnorm(N * Time, sd = 0.5)
#' res <- ols_panel(y, X, N, Time, ndraw = 2000, nomit = 500,
#'                  prior = list(model = 1, rval = 4))
#' print(res)
#' }
#'
#' @export
ols_panel <- function(y, X, N, Time, ndraw = 5500L, nomit = 1500L,
                      prior = list()) {


  start_time <- proc.time()

  # ---- Input validation ----
  y <- as.numeric(y)
  X <- as.matrix(X)
  nobs <- N * Time
  nvar <- ncol(X)

  stopifnot(
    length(y) == nobs,
    nrow(X) == nobs,
    ndraw > nomit,
    nomit >= 0
  )

  # ---- Parse priors ----
  model   <- prior$model %||% 0L
  rval    <- prior$rval  %||% 4
  nu      <- prior$nu    %||% 0
  d0      <- prior$d0    %||% 0
  thin    <- prior$thin  %||% 1L

  c_beta  <- prior$beta_mean %||% rep(0, nvar)
  C_beta  <- prior$beta_var  %||% (diag(nvar) * 1e12)

  stopifnot(model %in% 0:3)

  # Prior precision
  Q   <- solve(C_beta)
  Qpc <- Q %*% c_beta

  homo <- (rval == 0)

  # ---- Method label ----
  meth <- switch(as.character(model),
    "0" = "pols_g", "1" = "olssfe_g", "2" = "olstfe_g", "3" = "olsstfe_g"
  )

  # ---- Demean for fixed effects ----
  dm <- demean_panel(y, X, N, Time, model)
  ywith <- dm$ywith
  xwith <- dm$xwith

  # ---- Initialise MCMC storage ----
  n_save <- floor((ndraw - nomit) / thin)
  bsave <- matrix(0, nrow = n_save, ncol = nvar)
  ssave <- numeric(n_save)

  # State variables
  # V convention (following Geweke 1993):
  #   vi  = variance weight for obs i (large vi = MORE variance)
  #   V   = 1/vi = precision weight
  #   GLS pre-multiplication: xs = xwith * sqrt(V) = xwith * 1/sqrt(vi)
  #   This transforms y,X so that OLS on transformed data = GLS on original
  V    <- rep(1, nobs)
  vi   <- rep(1, nobs)
  sige <- 1
  vmean <- rep(0, nobs)

  # ---- MCMC loop ----
  save_idx <- 0L

  for (iter in seq_len(ndraw)) {

    if (homo) {
      # ---- Homoscedastic beta draw ----
      # Posterior: beta | y, sige ~ N(b0, sige * AI)
      #   where AI = (X'X + sige*Q)^{-1}, b0 = AI*(X'y + sige*Q*c)
      XtX <- crossprod(xwith)
      prec <- XtX + sige * Q
      AI   <- solve(prec)
      b0   <- AI %*% (crossprod(xwith, ywith) + sige * Qpc)
      bhat <- b0 + t(chol(sige * AI)) %*% rnorm(nvar)
    } else {
      # ---- Heteroscedastic beta draw ----
      # Pre-multiply data by sqrt(V) = sqrt(1/vi) for GLS
      # xs'*xs = X' * diag(V) * X = X' * diag(1/vi) * X
      sqrtV <- sqrt(V)               # sqrt(1/vi) = precision^{1/2}
      xs <- xwith * sqrtV            # element-wise: row i scaled by 1/sqrt(vi)
      ys <- ywith * sqrtV
      XtX <- crossprod(xs)
      prec <- XtX + sige * Q
      AI   <- solve(prec)
      b0   <- AI %*% (crossprod(xs, ys) + sige * Qpc)
      bhat <- b0 + t(chol(sige * AI)) %*% rnorm(nvar)
    }

    # ---- Draw sigma^2 ----
    # Posterior: sige | beta, V, y ~ IG(nu1/2, d1/2)
    # Sampled as sige = d1 / chi where chi ~ chi2(nu1)
    nu1 <- nu + nobs
    if (homo) {
      e <- ywith - xwith %*% bhat
    } else {
      e <- ys - xs %*% bhat           # residuals on V-scaled data
    }
    d1  <- d0 + as.numeric(crossprod(e))
    chi <- rchisq(1, df = nu1)
    sige <- d1 / chi

    # ---- Draw vi (heteroscedastic only) ----
    # Conditional: vi = (ei^2/sige + rval) / chi2(rval + 1)
    # where rval = prior df, "+1" from the single-observation likelihood
    # Then V = 1/vi (precision weights for next iteration's beta draw)
    if (!homo) {
      ev   <- as.numeric(ywith - xwith %*% bhat)  # residuals, force to plain vector
      chiv <- rchisq(nobs, df = rval + 1)
      vi   <- as.numeric((ev^2 / sige + rval) / chiv)
      V    <- 1 / vi
    }

    # ---- Save draws (after burn-in, with thinning) ----
    if (iter > nomit && ((iter - nomit) %% thin == 0L)) {
      save_idx <- save_idx + 1L
      bsave[save_idx, ] <- as.numeric(bhat)
      ssave[save_idx]   <- sige
      vmean <- vmean + vi
    }
  }

  vmean <- vmean / n_save

  # ---- Post-processing ----
  bmean <- colMeans(bsave)
  smean <- mean(ssave)
  sige  <- smean

  # t-statistics from posterior
  bstd  <- apply(bsave, 2, sd)
  tstat <- bmean / bstd

  # Fitted values and fixed effects
  fe_result <- .compute_ols_fe(
    y, X, ywith, xwith, bmean, sige,
    dm$meanny, dm$meannx, dm$meanty, dm$meantx,
    N, Time, nobs, model
  )

  # R-squared
  resid <- y - fe_result$yhat
  yme   <- y - mean(y)
  rsqr  <- 1 - sum(resid^2) / sum(yme^2)

  # Corr-squared (on demeaned scale)
  ywhat <- xwith %*% bmean
  res1  <- ywith - mean(ywith)
  res2  <- ywhat - mean(ywith)
  corr2 <- (sum(res1 * res2))^2 / (sum(res1^2) * sum(res2^2))

  # Log-likelihood at posterior mean
  e_lik <- ywith - xwith %*% bmean
  lik   <- -(nobs / 2) * log(2 * pi * sige) - sum(e_lik^2) / (2 * sige)

  # ---- Assemble output ----
  result <- list(
    beta      = bmean,
    sige      = smean,
    bdraw     = coda::mcmc(bsave, start = nomit + 1, thin = thin),
    sdraw     = coda::mcmc(matrix(ssave, ncol = 1), start = nomit + 1, thin = thin),
    vmean     = vmean,
    tstat     = tstat,
    yhat      = fe_result$yhat,
    resid     = resid,
    rsqr      = rsqr,
    corr2     = corr2,
    lik       = lik,
    intercept = fe_result$intercept,
    sfe       = fe_result$sfe,
    tfe       = fe_result$tfe,
    tsfe      = fe_result$tsfe,
    ttfe      = fe_result$ttfe,
    nobs      = nobs,
    nvar      = nvar,
    N         = N,
    Time      = Time,
    model     = model,
    meth      = meth,
    ndraw     = ndraw,
    nomit     = nomit,
    rval      = rval,
    cov       = cov(bsave),
    time      = (proc.time() - start_time)[3]
  )

  result <- .add_field_aliases(result)
  class(result) <- "spmixW"
  result
}


# ---- Internal: compute OLS fixed effects ----
.compute_ols_fe <- function(y, X, ywith, xwith, beta, sige,
                            meanny, meannx, meanty, meantx,
                            N, Time, nobs, model) {

  intercept <- NULL
  sfe <- tfe <- tsfe <- ttfe <- NULL

  XtXinv_diag <- function() {
    tryCatch(solve(crossprod(xwith)), error = function(e) {
      MASS::ginv(as.matrix(crossprod(xwith)))
    })
  }

  if (model == 0L) {
    yhat <- X %*% beta
  } else if (model == 1L) {
    intercept <- mean(y) - as.numeric(colMeans(X) %*% beta)
    sfe <- meanny - meannx %*% beta - intercept
    yhat <- as.numeric(X %*% beta) + rep(as.numeric(sfe), times = Time) + intercept

    XtXi <- XtXinv_diag()
    tsfe <- sfe / sqrt(sige / Time + diag(sige * meannx %*% XtXi %*% t(meannx)))
  } else if (model == 2L) {
    intercept <- mean(y) - as.numeric(colMeans(X) %*% beta)
    tfe <- meanty - meantx %*% beta - intercept
    yhat <- as.numeric(X %*% beta) + rep(as.numeric(tfe), each = N) + intercept

    XtXi <- XtXinv_diag()
    ttfe <- tfe / sqrt(sige / N + diag(sige * meantx %*% XtXi %*% t(meantx)))
  } else {
    intercept <- mean(y) - as.numeric(colMeans(X) %*% beta)
    sfe <- meanny - meannx %*% beta - intercept
    tfe <- meanty - meantx %*% beta - intercept
    yhat <- as.numeric(X %*% beta) +
      rep(as.numeric(sfe), times = Time) +
      rep(as.numeric(tfe), each = N) + intercept

    XtXi <- XtXinv_diag()
    tsfe <- sfe / sqrt(sige / Time + diag(sige * meannx %*% XtXi %*% t(meannx)))
    ttfe <- tfe / sqrt(sige / N + diag(sige * meantx %*% XtXi %*% t(meantx)))
  }

  list(
    yhat = as.numeric(yhat),
    intercept = intercept,
    sfe = as.numeric(sfe),
    tfe = as.numeric(tfe),
    tsfe = as.numeric(tsfe),
    ttfe = as.numeric(ttfe)
  )
}


# Null-coalescing operator (internal)
`%||%` <- function(x, y) if (is.null(x)) y else x
