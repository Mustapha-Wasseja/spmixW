#' Bayesian SDEM Panel Model with Fixed Effects
#'
#' Estimates a Spatial Durbin Error Model (SDEM) panel via MCMC. This is a
#' thin wrapper around \code{\link{sem_panel}} that augments X with WX.
#'
#' \deqn{y = X \beta + W X \theta + \text{FE} + u, \quad u = \lambda W u + \varepsilon}
#'
#' @inheritParams sar_panel
#'
#' @return An S3 object of class \code{"spmixW"} with SEM fields plus
#'   SDEM-specific effects:
#'   \describe{
#'     \item{direct}{Direct effects = coefficients on X (\eqn{\beta}).}
#'     \item{indirect}{Indirect effects = coefficients on WX (\eqn{\theta}).}
#'     \item{total}{Total = direct + indirect.}
#'   }
#'
#' @details
#' For SDEM, effects decomposition is straightforward (no spatial multiplier
#' on X, unlike SDM): direct effects are \eqn{\beta}, indirect effects are
#' \eqn{\theta}, and total = \eqn{\beta + \theta}. The spatial error
#' parameter \eqn{\lambda} does not affect the effects decomposition.
#'
#' @references
#' Debarsy, N. and LeSage, J.P. (2021). "Bayesian model averaging for
#' spatial autoregressive models based on convex combinations of different
#' types of connectivity matrices." \emph{Journal of Business & Economic
#' Statistics}, 40(2), 547-558.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' N <- 30; Time <- 10; lambda_true <- 0.5
#' coords <- cbind(runif(N), runif(N))
#' W <- as.matrix(normw(make_knw(coords, k = 5, row_normalise = FALSE)))
#' Wbig <- kronecker(diag(Time), W)
#' X <- matrix(rnorm(N * Time * 2), ncol = 2)
#' WX <- Wbig %*% X
#' u <- solve(diag(N*Time) - lambda_true * Wbig) %*% rnorm(N*Time)
#' y <- X %*% c(1, -0.5) + WX %*% c(0.3, 0.2) + u
#' res <- sdem_panel(as.numeric(y), X, W, N, Time, ndraw = 5000, nomit = 2000)
#' print(res)
#' }
#'
#' @export
sdem_panel <- function(y, X, W, N, Time, ndraw = 5500L, nomit = 1500L,
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

  has_intercept <- all(abs(X[, 1] - X[1, 1]) < 1e-10)
  if (has_intercept) {
    xmat <- cbind(X, WX[, -1, drop = FALSE])
  } else {
    xmat <- cbind(X, WX)
  }

  # ---- Call sem_panel with augmented X ----
  result <- sem_panel(y, xmat, W, N, Time, ndraw, nomit, prior)

  # ---- Update method label ----
  result$meth <- switch(as.character(prior$model %||% 0L),
    "0" = "psdem_g", "1" = "sdemsfe_g", "2" = "sdemtfe_g", "3" = "sdemstfe_g"
  )

  # ---- SDEM effects: direct = beta, indirect = theta, total = beta + theta ----
  bsave <- as.matrix(result$bdraw)
  ndrawsg <- nrow(bsave)

  if (has_intercept) {
    p <- nvar_orig - 1
    beta_idx  <- 2:nvar_orig
    theta_idx <- (nvar_orig + 1):ncol(bsave)
  } else {
    p <- nvar_orig
    beta_idx  <- 1:nvar_orig
    theta_idx <- (nvar_orig + 1):ncol(bsave)
  }

  direct_draws   <- bsave[, beta_idx, drop = FALSE]
  indirect_draws <- bsave[, theta_idx, drop = FALSE]
  total_draws    <- direct_draws + indirect_draws

  .effects_summary <- function(draws) {
    p_vars <- ncol(draws)
    out <- matrix(0, p_vars, 5)
    colnames(out) <- c("Mean", "t-stat", "p-value", "Lower05", "Upper95")
    for (j in seq_len(p_vars)) {
      m <- mean(draws[, j]); s <- sd(draws[, j])
      t_val <- m / s
      p_val <- 2 * (1 - pt(abs(t_val), df = nrow(draws) - 1))
      bounds <- quantile(draws[, j], c(0.025, 0.975))
      out[j, ] <- c(m, t_val, p_val, bounds[1], bounds[2])
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
