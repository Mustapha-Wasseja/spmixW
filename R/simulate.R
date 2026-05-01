#' Simulate Spatial Panel Data
#'
#' Generates synthetic spatial panel data from a SAR/SDM/SEM/SDEM DGP,
#' returning a ready-to-use data frame for \code{\link{spmodel}}.
#'
#' @param N Integer: number of cross-sectional units (regions).
#' @param Time Integer: number of time periods (use 1 for cross-sectional).
#' @param W A single N x N weight matrix, or a list of M weight matrices
#'   for convex combination DGPs.
#' @param gamma Numeric vector of convex weights (required when \code{W} is
#'   a list, must sum to 1). Ignored when \code{W} is a single matrix.
#' @param rho Numeric: spatial autoregressive parameter (default 0).
#' @param beta Numeric vector of true regression coefficients (length k).
#' @param theta Numeric vector of WX coefficients for SDM/SDEM DGPs
#'   (default NULL). Must have same length as \code{beta} if provided.
#' @param sigma2 Numeric: error variance (default 1).
#' @param effects Character: fixed-effects specification.
#'   \code{"none"}, \code{"region"}, \code{"time"}, or \code{"twoway"}
#'   (default \code{"twoway"}).
#' @param seed Integer or NULL: random seed for reproducibility.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{region}{Integer region identifier (1 to N).}
#'     \item{year}{Integer time period identifier (omitted if \code{Time = 1}).}
#'     \item{y}{Simulated response variable.}
#'     \item{x1, x2, ...}{Simulated predictor variables (standard normal).}
#'   }
#'   The data frame is sorted by time then region, ready for
#'   \code{\link{spmodel}}.
#'
#' @details
#' The DGP is:
#' \deqn{y = (I_{NT} - \rho W_c)^{-1} (X \beta + W_c X \theta + \text{FE} + \varepsilon)}
#' where \eqn{W_c = \sum \gamma_m W_m} (or just W if a single matrix),
#' \eqn{\theta} is omitted for SAR/SEM DGPs, and FE are generated as
#' \eqn{\mu_i = i/N} (region) and \eqn{\nu_t = t/T} (time).
#'
#' @examples
#' \donttest{
#' coords <- cbind(runif(80), runif(80))
#' W1 <- make_knw(coords, k = 4)
#' W2 <- make_knw(coords, k = 8)
#'
#' panel <- simulate_panel(
#'   N = 80, T = 10,
#'   W = list(W1, W2),
#'   gamma = c(0.7, 0.3),
#'   rho = 0.5,
#'   beta = c(1, -1),
#'   seed = 42
#' )
#'
#' res <- spmodel(y ~ x1 + x2, data = panel,
#'                W = list(geography = W1, trade = W2),
#'                model = "sar",
#'                id = "region", time = "year",
#'                effects = "twoway",
#'                ndraw = 8000, nomit = 2000)
#' print(res)
#' }
#'
#' @export
simulate_panel <- function(N, Time = 10L, W, gamma = NULL, rho = 0,
                           beta, theta = NULL, sigma2 = 1,
                           effects = "twoway", seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  k <- length(beta)
  nobs <- N * Time
  is_list_W <- is.list(W) && !inherits(W, "Matrix")

  # ---- Input validation ----
  if (is_list_W) {
    M <- length(W)
    if (is.null(gamma)) stop("gamma is required when W is a list of weight matrices.")
    stopifnot(length(gamma) == M, all(gamma >= 0))
    if (abs(sum(gamma) - 1) > 1e-10) stop("gamma must sum to 1.")
    for (m in seq_len(M)) {
      stopifnot(nrow(W[[m]]) == N, ncol(W[[m]]) == N)
    }
    # Build Wc
    Wc <- Reduce("+", mapply(function(g, Wm) g * as.matrix(Wm), gamma, W,
                              SIMPLIFY = FALSE))
  } else {
    stopifnot(nrow(W) == N, ncol(W) == N)
    Wc <- as.matrix(W)
  }

  if (!is.null(theta)) {
    stopifnot(length(theta) == k)
  }

  # ---- Generate X ----
  X <- matrix(rnorm(nobs * k), ncol = k)

  # ---- Fixed effects ----
  fe_vec <- numeric(nobs)
  if (effects %in% c("region", "twoway")) {
    fe_vec <- fe_vec + rep((1:N) / N, times = Time)
  }
  if (effects %in% c("time", "twoway")) {
    fe_vec <- fe_vec + rep((1:Time) / Time, each = N)
  }

  # ---- Build Wbig and systematic component ----
  Wbig <- kronecker(diag(Time), Wc)
  mu <- X %*% beta + fe_vec

  # Add WX*theta if SDM/SDEM

  if (!is.null(theta)) {
    WX <- Wbig %*% X
    mu <- mu + WX %*% theta
  }

  # ---- Generate y ----
  evec <- rnorm(nobs, sd = sqrt(sigma2))
  A <- diag(nobs) - rho * Wbig
  y <- as.numeric(solve(A, mu + evec))

  # ---- Assemble data frame ----
  df <- data.frame(region = rep(1:N, times = Time))
  if (Time > 1) {
    df$year <- rep(seq_len(Time), each = N)
  }
  df$y <- y
  for (j in seq_len(k)) {
    df[[paste0("x", j)]] <- X[, j]
  }

  df
}
