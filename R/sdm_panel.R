#' Bayesian SDM Panel Model with Fixed Effects
#'
#' Estimates a Spatial Durbin Model (SDM) panel via MCMC. This is a thin
#' wrapper around \code{\link{sar_panel}} that augments the X matrix with
#' WX before calling the SAR estimator, then computes SDM-specific effects.
#'
#' \deqn{y = \rho W y + X \beta + W X \theta + \text{FE} + \varepsilon}
#'
#' @inheritParams sar_panel
#'
#' @return An S3 object of class \code{"spmixW"} with all SAR fields.
#'   The \code{beta} vector contains \eqn{[\beta', \theta']'} (length \eqn{2k}
#'   or \eqn{2k - 1} if X has an intercept). Effects are SDM-style: both
#'   \eqn{\beta} and \eqn{\theta} contribute to the spatial multiplier.
#'
#' @details
#' Internally, this function augments X with WX and calls
#' \code{\link{sar_panel}}, following the delegation pattern of
#' LeSage and Pace (2009, Ch. 10). The SDM effects computation differs from SAR because the
#' spatial multiplier \eqn{(I - \rho W)^{-1}} acts on both \eqn{\beta I}
#' and \eqn{\theta W}, yielding different direct/indirect decompositions.
#'
#' @references
#' LeSage, J.P. and Pace, R.K. (2009). \emph{Introduction to Spatial
#' Econometrics}. Taylor & Francis/CRC Press.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30; Time <- 10; rho_true <- 0.4
#' coords <- cbind(runif(N), runif(N))
#' W <- as.matrix(normw(make_knw(coords, k = 5, row_normalise = FALSE)))
#' Wbig <- kronecker(diag(Time), W)
#' X <- matrix(rnorm(N * Time * 2), ncol = 2)
#' WX <- Wbig %*% X
#' y <- solve(diag(N*Time) - rho_true * Wbig) %*%
#'        (X %*% c(1, -0.5) + WX %*% c(0.3, 0.2) + rnorm(N*Time))
#' res <- sdm_panel(as.numeric(y), X, W, N, Time, ndraw = 5000, nomit = 2000)
#' print(res)
#' }
#'
#' @export
sdm_panel <- function(y, X, W, N, Time, ndraw = 5500L, nomit = 1500L,
                      prior = list()) {

  X <- as.matrix(X)
  nobs <- N * Time
  nvar_orig <- ncol(X)

  # ---- Handle W dimensions ----
  W_sp <- Matrix::Matrix(W, sparse = TRUE)
  if (nrow(W_sp) == N && ncol(W_sp) == N) {
    Wbig <- kronecker(Matrix::Diagonal(Time), W_sp)
  } else if (nrow(W_sp) == nobs) {
    Wbig <- W_sp
  } else {
    stop("W must be N x N or NT x NT")
  }

  # ---- Augment X with WX ----
  WX <- as.matrix(Wbig %*% X)

  # Check for intercept: if first column is constant, don't include WX
  # for that column (WX of a constant = constant for row-normalised W)
  has_intercept <- all(abs(X[, 1] - X[1, 1]) < 1e-10)
  if (has_intercept) {
    xmat <- cbind(X, WX[, -1, drop = FALSE])
  } else {
    xmat <- cbind(X, WX)
  }

  # ---- Call sar_panel with augmented X ----
  result <- sar_panel(y, xmat, W, N, Time, ndraw, nomit, prior)

  # ---- Update method label ----
  result$meth <- switch(as.character(prior$model %||% 0L),
    "0" = "psdm_g", "1" = "sdmsfe_g", "2" = "sdmtfe_g", "3" = "sdmstfe_g"
  )

  # ---- Recompute SDM-style effects ----
  # For SDM: total(j) = sum_r rho^r * (beta_j + theta_j)
  #          direct(j) = sum_r rho^r * (beta_j * tr(W^r)/N + theta_j * tr(W^{r+1})/N)
  bsave <- as.matrix(result$bdraw)
  psave <- as.numeric(result$pdraw)
  ndrawsg <- nrow(bsave)

  # Stochastic trace estimation
  uiter <- 50L; maxorderu <- 100L
  rv <- matrix(rnorm(nobs * uiter), nrow = nobs, ncol = uiter)
  tracew <- numeric(maxorderu)
  wjjju <- rv
  for (j in seq_len(maxorderu)) {
    wjjju <- as.matrix(Wbig %*% wjjju)
    tracew[j] <- mean(colMeans(rv * wjjju))
  }
  tracew[1] <- 0
  tracew[2] <- sum(Wbig * Wbig) / nobs

  trs <- c(1, tracew)
  ntrs <- length(trs)
  ree <- 0:(ntrs - 1)

  # Figure out which columns of bsave correspond to beta vs theta
  if (has_intercept) {
    p <- nvar_orig - 1  # number of non-intercept X variables
    beta_idx <- 2:nvar_orig           # indices in bsave for beta (skip intercept)
    theta_idx <- (nvar_orig + 1):ncol(bsave)  # indices for theta
  } else {
    p <- nvar_orig
    beta_idx <- 1:nvar_orig
    theta_idx <- (nvar_orig + 1):ncol(bsave)
  }

  total_draws    <- matrix(0, ndrawsg, p)
  direct_draws   <- matrix(0, ndrawsg, p)
  indirect_draws <- matrix(0, ndrawsg, p)

  # SDM effects: for each variable j, both beta_j and theta_j contribute
  # trs[1] = 1, trs[2] = 0 (trace W / N), trs[3] = tr(W^2)/N, ...
  # Trace of W shifted by one power for theta terms
  trs_shifted <- c(tracew, 0)  # tr(W^1)/N, tr(W^2)/N, ..., 0
  # Length should match ntrs
  trs_shifted <- c(tracew[1:min(maxorderu, ntrs-1)],
                   rep(0, max(0, ntrs - maxorderu)))

  for (i in seq_len(ndrawsg)) {
    rmat <- psave[i]^ree
    for (j in seq_len(p)) {
      bj <- bsave[i, beta_idx[j]]
      tj <- bsave[i, theta_idx[j]]

      # Total: (beta_j + theta_j) * sum(rho^r)
      total_val <- sum((bj + tj) * rmat)

      # Direct: beta_j * sum(tr(W^r)/N * rho^r) + theta_j * sum(tr(W^{r+1})/N * rho^r)
      direct_val <- sum(bj * trs * rmat) + sum(tj * trs_shifted[1:ntrs] * rmat)

      indirect_val <- total_val - direct_val

      total_draws[i, j]    <- total_val
      direct_draws[i, j]   <- direct_val
      indirect_draws[i, j] <- indirect_val
    }
  }

  .effects_summary <- function(draws) {
    p_vars <- ncol(draws)
    out <- matrix(0, p_vars, 5)
    colnames(out) <- c("Mean", "t-stat", "p-value", "Lower05", "Upper95")
    for (jj in seq_len(p_vars)) {
      m <- mean(draws[, jj]); s <- sd(draws[, jj])
      t_val <- m / s
      p_val <- 2 * (1 - pt(abs(t_val), df = nrow(draws) - 1))
      bounds <- quantile(draws[, jj], c(0.025, 0.975))
      out[jj, ] <- c(m, t_val, p_val, bounds[1], bounds[2])
    }
    out
  }

  result$direct   <- .effects_summary(direct_draws)
  result$indirect <- .effects_summary(indirect_draws)
  result$total    <- .effects_summary(total_draws)
  result$direct_draws   <- direct_draws
  result$indirect_draws <- indirect_draws
  result$total_draws    <- total_draws
  result$nvar_orig <- nvar_orig
  result$p <- p

  result <- .add_field_aliases(result)
  result
}
