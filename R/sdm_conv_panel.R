#' Bayesian SDM Panel with Convex Combination of Weight Matrices
#'
#' Estimates a SDM panel model with \eqn{W_c = \sum \gamma_m W_m}:
#' \deqn{y = \rho W_c y + X\beta + W_1 X \theta_1 + \ldots + W_M X \theta_M + \varepsilon}
#'
#' Unlike the standard SDM wrapper, the WX terms are pre-computed for each
#' individual \eqn{W_m} since the convex combination \eqn{W_c} changes at
#' each MCMC iteration. The augmented X is \code{[X, W1*X, W2*X, ..., WM*X]}.
#'
#' @inheritParams sar_conv_panel
#'
#' @return An S3 object of class \code{"spmixW"} with convex combination fields.
#'
#' @details
#' Internally augments X with \code{[X, W1*X, W2*X, ..., WM*X]} and
#' then calls \code{\link{sar_conv_panel}} for estimation. The effects decomposition accounts for all M
#' sets of WX coefficients.
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
sdm_conv_panel <- function(y, X, Wlist, N, Time, ndraw = 25000L,
                           nomit = 5000L, prior = list()) {

  X <- as.matrix(X)
  nobs <- N * Time
  nvar_orig <- ncol(X)
  M <- length(Wlist)

  # Build big W matrices
  Wbigs <- lapply(Wlist, function(Wm) {
    Wm <- Matrix::Matrix(Wm, sparse = TRUE)
    if (nrow(Wm) == N) kronecker(Matrix::Diagonal(Time), Wm) else Wm
  })

  # Augment X with W_m * X for each m
  # SDM conv uses [X, W1*X, W2*X, ..., WM*X]
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

  # Call sar_conv_panel with augmented X
  result <- sar_conv_panel(y, xmat, Wlist, N, Time, ndraw, nomit, prior)

  # Update method
  result$meth <- switch(as.character(prior$model %||% 0L),
    "0" = "psdm_conv_g", "1" = "sdmsfe_conv_g",
    "2" = "sdmtfe_conv_g", "3" = "sdmstfe_conv_g"
  )
  result$nvar_orig <- nvar_orig

  result <- .add_field_aliases(result)
  result
}
