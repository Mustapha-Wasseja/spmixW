#' Panel Demeaning for Fixed Effects
#'
#' Removes fixed effects from panel data via the within transformation.
#' Data must be sorted by time then by spatial unit: all N regions in period 1,
#' then all N regions in period 2, etc. (i.e., \eqn{y_{it}} is stored at
#' position \eqn{i + (t-1) \times N}).
#'
#' @param y Numeric vector of length \eqn{NT}: dependent variable.
#' @param X Numeric matrix of dimension \eqn{NT \times k}: explanatory variables.
#' @param N Integer: number of cross-sectional units (regions).
#' @param Time Integer: number of time periods.
#' @param model Integer controlling fixed-effects specification:
#'   \describe{
#'     \item{0}{No demeaning (pooled).}
#'     \item{1}{Region (spatial) fixed effects only.}
#'     \item{2}{Time period fixed effects only.}
#'     \item{3}{Both region and time period fixed effects (two-way FE).}
#'   }
#'
#' @return A list with elements:
#'   \describe{
#'     \item{ywith}{Demeaned \code{y} (length \eqn{NT}).}
#'     \item{xwith}{Demeaned \code{X} (\eqn{NT \times k}).}
#'     \item{meanny}{Region means of y (length N); zero if \code{model \%in\% c(0,2)}.}
#'     \item{meannx}{Region means of X (\eqn{N \times k}); zero if \code{model \%in\% c(0,2)}.}
#'     \item{meanty}{Time means of y (length T); zero if \code{model \%in\% c(0,1)}.}
#'     \item{meantx}{Time means of X (\eqn{T \times k}); zero if \code{model \%in\% c(0,1)}.}
#'   }
#'
#' @details
#' The demeaning follows the standard panel FE within transformation:
#' \itemize{
#'   \item Model 1 (region FE): subtract region means across time.
#'   \item Model 2 (time FE): subtract time-period means across regions.
#'   \item Model 3 (two-way FE): subtract both region and time means,
#'         then add back the grand mean (to avoid double subtraction).
#' }
#'
#' @references
#' LeSage, J.P. and Pace, R.K. (2009). \emph{Introduction to Spatial
#' Econometrics}. Taylor & Francis/CRC Press.
#'
#' @examples
#' N <- 10; Time <- 5; k <- 2
#' set.seed(42)
#' y <- rnorm(N * Time)
#' X <- matrix(rnorm(N * Time * k), ncol = k)
#' dm <- demean_panel(y, X, N, Time, model = 3)
#' # Verify the demeaned y has (approximately) zero region and time means
#'
#' @export
demean_panel <- function(y, X, N, Time, model = 0L) {

  y <- as.numeric(y)
  X <- as.matrix(X)
  nobs <- N * Time
  nvar <- ncol(X)

  stopifnot(
    length(y) == nobs,
    nrow(X) == nobs,
    model %in% 0:3
  )

  meanny <- numeric(N)
  meannx <- matrix(0, nrow = N, ncol = nvar)
  meanty <- numeric(Time)
  meantx <- matrix(0, nrow = Time, ncol = nvar)

  # Index helper: observation for unit i at time t is at position i + (t-1)*N
  # For region means: collect all T observations for unit i
  if (model %in% c(1L, 3L)) {
    # Reshape y into N x T matrix for efficient computation
    ymat <- matrix(y, nrow = N, ncol = Time)  # columns = time periods
    meanny <- rowMeans(ymat)

    # Reshape each column of X similarly
    for (j in seq_len(nvar)) {
      xmat_j <- matrix(X[, j], nrow = N, ncol = Time)
      meannx[, j] <- rowMeans(xmat_j)
    }
  }

  # For time means: collect all N observations in period t
  if (model %in% c(2L, 3L)) {
    ymat <- matrix(y, nrow = N, ncol = Time)
    meanty <- colMeans(ymat)

    for (j in seq_len(nvar)) {
      xmat_j <- matrix(X[, j], nrow = N, ncol = Time)
      meantx[, j] <- colMeans(xmat_j)
    }
  }

  # Apply demeaning
  if (model == 0L) {
    ywith <- y
    xwith <- X
  } else if (model == 1L) {
    # Subtract region means: kron(ones_T, meanny)
    ywith <- y - rep(meanny, times = Time)
    xwith <- X - do.call(rbind, replicate(Time, meannx, simplify = FALSE))
  } else if (model == 2L) {
    # Subtract time means: kron(meanty, ones_N)
    ywith <- y - rep(meanty, each = N)
    xwith <- X - meantx[rep(seq_len(Time), each = N), , drop = FALSE]
  } else {
    # model == 3: two-way within transformation
    # y_it_within = y_it - ybar_i. - ybar_.t + ybar_..
    # where ybar_i. = region mean, ybar_.t = time mean, ybar_.. = grand mean
    grand_mean_y <- mean(y)
    grand_mean_x <- colMeans(X)

    ywith <- y -
      rep(meanny, times = Time) -       # kron(ones_T, meanny)
      rep(meanty, each = N) +           # kron(meanty, ones_N)
      grand_mean_y                       # scalar, recycled to length NT

    xwith <- X -
      do.call(rbind, replicate(Time, meannx, simplify = FALSE)) -
      meantx[rep(seq_len(Time), each = N), , drop = FALSE] +
      matrix(grand_mean_x, nrow = nobs, ncol = nvar, byrow = TRUE)
  }

  list(
    ywith  = ywith,
    xwith  = xwith,
    meanny = meanny,
    meannx = meannx,
    meanty = meanty,
    meantx = meantx
  )
}
