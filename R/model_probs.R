#' Posterior Model Probabilities from Log-Marginal Likelihoods
#'
#' Converts a vector of log-marginal likelihoods to posterior model
#' probabilities assuming equal prior model probabilities.
#'
#' @param lmarginal Numeric vector of log-marginal likelihood values,
#'   one per model.
#'
#' @return A numeric vector of posterior model probabilities summing to 1.
#'
#' @details
#' Assumes equal prior probabilities across models. Computes:
#' \deqn{p(M_j | y) = \frac{\exp(\ell_j - \max(\ell))}{\sum_i \exp(\ell_i - \max(\ell))}}
#' where \eqn{\ell_j} is the log-marginal likelihood for model \eqn{j}.
#' The subtraction of the maximum prevents numerical overflow.
#'
#' @references
#' LeSage, J.P. (2014). "What Regional Scientists Need to Know about Spatial
#' Econometrics." \emph{Review of Regional Studies}, 44(1), 13-32.
#'
#' @examples
#' lm <- c(-100, -98, -105)
#' model_probs(lm)
#'
#' @export
model_probs <- function(lmarginal) {

  lmarginal <- as.numeric(lmarginal)
  stopifnot(length(lmarginal) >= 1, all(is.finite(lmarginal)))

  # Scale before exponentiating to avoid overflow
  adj <- max(lmarginal)
  x <- exp(lmarginal - adj)

  x / sum(x)
}
