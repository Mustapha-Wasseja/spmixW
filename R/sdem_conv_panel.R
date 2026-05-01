#' Bayesian SDEM Panel with Convex Combination of Weight Matrices
#'
#' Estimates a SDEM panel model with \eqn{W_c = \sum \gamma_m W_m}:
#' \deqn{y = X\beta + W_1 X \theta_1 + \ldots + W_M X \theta_M + u, \quad
#' u = \lambda W_c u + \varepsilon}
#'
#' Augments X with \code{[X, W1*X, W2*X, ..., WM*X]} and calls
#' \code{\link{sem_conv_panel}}.
#'
#' @inheritParams sar_conv_panel
#'
#' @return An S3 object of class \code{"spmixW"} with convex combination fields
#'   plus SDEM-specific effects.
#'
#' @details
#' For SDEM, direct effects are the \eqn{\beta} coefficients, indirect effects
#' are the \eqn{\theta_m} coefficients summed across W matrices, and total =
#' direct + indirect. The spatial error parameter does not affect effects.
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
sdem_conv_panel <- function(y, X, Wlist, N, Time, ndraw = 25000L,
                            nomit = 5000L, prior = list()) {

  X <- as.matrix(X)
  nobs <- N * Time
  nvar_orig <- ncol(X)
  M <- length(Wlist)

  # Build big W
  Wbigs <- lapply(Wlist, function(Wm) {
    Wm <- Matrix::Matrix(Wm, sparse = TRUE)
    if (nrow(Wm) == N) kronecker(Matrix::Diagonal(Time), Wm) else Wm
  })

  # Augment X
  has_intercept <- all(abs(X[, 1] - X[1, 1]) < 1e-10)
  xmat <- X
  for (m in seq_len(M)) {
    WmX <- as.matrix(Wbigs[[m]] %*% X)
    if (has_intercept) {
      xmat <- cbind(xmat, WmX[, -1, drop = FALSE])
    } else {
      xmat <- cbind(xmat, WmX)
    }
  }

  # Call sem_conv_panel
  result <- sem_conv_panel(y, xmat, Wlist, N, Time, ndraw, nomit, prior)

  result$meth <- switch(as.character(prior$model %||% 0L),
    "0" = "psdem_conv_g", "1" = "sdemsfe_conv_g",
    "2" = "sdemtfe_conv_g", "3" = "sdemstfe_conv_g"
  )
  result$nvar_orig <- nvar_orig

  # SDEM effects: direct = beta, indirect = sum of theta_m across W's
  bsave <- as.matrix(result$bdraw)
  ndrawsg <- nrow(bsave)
  p_x <- if (has_intercept) nvar_orig - 1 else nvar_orig

  if (has_intercept) {
    beta_idx <- 2:nvar_orig
  } else {
    beta_idx <- 1:nvar_orig
  }
  # theta indices: after X columns, M blocks of p_x columns each
  theta_start <- nvar_orig + 1
  theta_end <- ncol(bsave)

  direct_draws <- bsave[, beta_idx, drop = FALSE]
  # Sum indirect effects across all M weight matrices
  indirect_draws <- matrix(0, ndrawsg, p_x)
  for (m_idx in seq_len(M)) {
    cols <- theta_start + (m_idx - 1) * p_x + (seq_len(p_x) - 1)
    cols <- cols[cols <= theta_end]
    if (length(cols) == p_x) {
      indirect_draws <- indirect_draws + bsave[, cols, drop = FALSE]
    }
  }
  total_draws <- direct_draws + indirect_draws

  .effects_summary <- function(draws) {
    pv <- ncol(draws)
    out <- matrix(0, pv, 5)
    colnames(out) <- c("Mean", "t-stat", "p-value", "Lower05", "Upper95")
    for (j in seq_len(pv)) {
      m <- mean(draws[, j]); s <- sd(draws[, j])
      tv <- m / s
      pval <- 2 * (1 - pt(abs(tv), df = nrow(draws) - 1))
      bounds <- quantile(draws[, j], c(0.025, 0.975))
      out[j, ] <- c(m, tv, pval, bounds[1], bounds[2])
    }
    out
  }

  result$direct   <- .effects_summary(direct_draws)
  result$indirect <- .effects_summary(indirect_draws)
  result$total    <- .effects_summary(total_draws)
  result$direct_draws   <- direct_draws
  result$indirect_draws <- indirect_draws
  result$total_draws    <- total_draws
  result$p <- p_x

  result <- .add_field_aliases(result)
  result
}
