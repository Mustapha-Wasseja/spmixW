#' Taylor Series Log-Determinant for Convex Combination of Weight Matrices
#'
#' Pre-computes trace cross-terms for the Taylor series approximation to
#' \eqn{\log|I - \rho W_c(\gamma)|} where
#' \eqn{W_c = \sum_{m=1}^{M} \gamma_m W_m}.
#'
#' @param Wlist A list of M sparse or dense weight matrices, each \eqn{N \times N}.
#' @param max_order Integer: maximum Taylor order (default 4, supports 2-8).
#'   Higher orders give more accurate log-det approximations at the cost of
#'   more pre-computation time. The number of trace terms at order p is
#'   \eqn{M^p}, so for M=3 and order 8 this is 6561 terms.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{traces}{A named list where \code{traces[[p]]} is a vector of
#'       length \eqn{M^p} containing all \eqn{\text{tr}(W_{i_1} \cdots W_{i_p})}
#'       cross-terms for order p (p = 2, ..., max_order).}
#'     \item{max_order}{The maximum Taylor order stored.}
#'     \item{M}{Number of weight matrices.}
#'     \item{N}{Dimension of each weight matrix.}
#'   }
#'
#' @details
#' The Taylor approximation is:
#' \deqn{\log|I - \rho W_c| \approx -\sum_{p=2}^{\text{order}}
#' \frac{\rho^p}{p} \text{tr}(W_c^p)}
#' (the p=1 term vanishes for zero-diagonal W).
#'
#' Since \eqn{W_c = \sum \gamma_m W_m}, the trace expands via multinomial:
#' \deqn{\text{tr}(W_c^p) = \sum_{i_1, \ldots, i_p}
#' \gamma_{i_1} \cdots \gamma_{i_p} \text{tr}(W_{i_1} \cdots W_{i_p})}
#'
#' These cross-terms are pre-computed once (expensive for high orders) and
#' then rapidly evaluated for any \eqn{(\rho, \gamma)} via Kronecker products.
#'
#' For M=3 weight matrices, the number of trace terms by order:
#' \tabular{rr}{
#'   Order \tab Terms \cr
#'   4 \tab 81 \cr
#'   5 \tab 243 \cr
#'   6 \tab 729 \cr
#'   7 \tab 2187 \cr
#'   8 \tab 6561 \cr
#' }
#'
#' @references
#' Debarsy, N. and LeSage, J.P. (2021). "Bayesian model averaging for spatial
#' autoregressive models based on convex combinations of different types of
#' connectivity matrices." \emph{Journal of Business & Economic Statistics},
#' 40(2), 547-558.
#'
#' @examples
#' \donttest{
#' N <- 20
#' W1 <- matrix(0, N, N); W2 <- matrix(0, N, N)
#' for (i in 1:N) { W1[i, (i %% N) + 1] <- 1; W2[i, ((i-2) %% N) + 1] <- 1 }
#' W1 <- W1 / rowSums(W1); W2 <- W2 / rowSums(W2)
#' traces <- log_det_taylor(list(W1, W2), max_order = 6)
#' eval_taylor_lndet(traces, rho = 0.5, gamma = c(0.7, 0.3))
#' }
#'
#' @export
log_det_taylor <- function(Wlist, max_order = 4L) {

  stopifnot(is.list(Wlist), length(Wlist) >= 1)
  M <- length(Wlist)
  N <- nrow(Wlist[[1]])
  stopifnot(all(sapply(Wlist, function(w) nrow(w) == N && ncol(w) == N)))
  stopifnot(max_order >= 2L, max_order <= 12L)

  Ws <- lapply(Wlist, as.matrix)

  # For each order p = 2, ..., max_order, compute all M^p trace terms.
  # We use an iterative approach: maintain a list of "accumulated products"
  # from the previous order and extend by one more W factor.
  #
  # At order p, we need tr(W_{i1} * W_{i2} * ... * W_{ip}) for all
  # (i1, ..., ip) in {1..M}^p.
  #
  # Strategy: at order p-1 we have M^{p-1} accumulated matrix products.
  # For order p, multiply each by each W_m and take trace against W_m.
  # But storing M^{p-1} full N x N matrices is expensive.
  #
  # More efficient: build incrementally. At each order, we have the
  # products W_{i1}...W_{i_{p-1}} and extend to W_{i1}...W_{i_p}.

  trace_list <- list()

  # Order 2: tr(Wi * Wj) for all (i,j) — M^2 terms
  # Start with "accumulated products" of order 1 = the W matrices themselves
  # Then extend to order 2, 3, ...

  # We store accumulated products as a flat list of M^{p-1} matrices
  # Each matrix is the product W_{i1} * W_{i2} * ... * W_{i_{p-1}}
  # indexed in row-major order of (i1, ..., i_{p-1})

  # Order 1: the M matrices themselves
  accum <- Ws  # length M

  for (p in 2:max_order) {
    n_prev <- length(accum)  # M^{p-1}
    n_new <- n_prev * M      # M^p
    Tp <- numeric(n_new)

    if (p < max_order) {
      # Need to store new accumulated products for next order
      new_accum <- vector("list", n_new)
    }

    cnt <- 1L
    for (a in seq_len(n_prev)) {
      Aprev <- accum[[a]]  # product W_{i1}...W_{i_{p-1}}
      for (m in seq_len(M)) {
        # tr(Aprev * Wm) = sum(Aprev * t(Wm))
        Tp[cnt] <- sum(Aprev * t(Ws[[m]]))
        if (p < max_order) {
          new_accum[[cnt]] <- Aprev %*% Ws[[m]]
        }
        cnt <- cnt + 1L
      }
    }

    trace_list[[as.character(p)]] <- Tp

    if (p < max_order) {
      accum <- new_accum
    }
  }

  # For backward compatibility, also store T2, T3, T4 as named elements
  result <- list(
    traces = trace_list,
    max_order = max_order,
    M = M,
    N = N
  )

  # Backward compat fields
  if ("2" %in% names(trace_list)) result$T2 <- trace_list[["2"]]
  if ("3" %in% names(trace_list)) result$T3 <- trace_list[["3"]]
  if ("4" %in% names(trace_list)) result$T4 <- trace_list[["4"]]

  result
}


#' Evaluate Taylor Series Log-Determinant at Given (rho, gamma)
#'
#' Rapidly evaluates \eqn{\log|I - \rho W_c(\gamma)|} using pre-computed
#' trace terms from \code{\link{log_det_taylor}}.
#'
#' @param traces Output from \code{\link{log_det_taylor}}.
#' @param rho Numeric scalar: spatial autoregressive parameter.
#' @param gamma Numeric vector of length M: convex weights (must sum to 1).
#' @param order Integer: Taylor order to use (default: use all available).
#'   Must be <= \code{traces$max_order}.
#'
#' @return Scalar: approximate \eqn{\log|I - \rho W_c(\gamma)|}.
#'
#' @export
eval_taylor_lndet <- function(traces, rho, gamma, order = NULL) {

  gamma <- as.numeric(gamma)
  M <- traces$M

  if (is.null(order)) order <- traces$max_order
  stopifnot(order >= 2, order <= traces$max_order, length(gamma) == M)

  # Build the Kronecker products of gamma incrementally:
  #   g_p = kron(g_{p-1}, gamma) has length M^p
  # Then tr(Wc^p) = g_p' * Tp

  result <- 0
  g_prev <- gamma  # g_1 = gamma (length M)

  for (p in 2:order) {
    g_p <- as.numeric(kronecker(g_prev, gamma))  # length M^p
    Tp <- traces$traces[[as.character(p)]]
    tr_Wcp <- sum(g_p * Tp)
    result <- result - (rho^p / p) * tr_Wcp
    g_prev <- g_p
  }

  result
}
