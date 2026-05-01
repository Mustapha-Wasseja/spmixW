#' Print Method for spmixW_bma Objects
#'
#' @param x An object of class \code{"spmixW_bma"}.
#' @param digits Integer: significant digits.
#' @param ... Further arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.spmixW_bma <- function(x, digits = 4, ...) {

  cat("========================================\n")
  cat("Bayesian Model Averaging Results\n")
  cat(sprintf("Method: %s\n", x$meth))
  cat("========================================\n")
  cat(sprintf("N = %d, T = %d, M = %d weight matrices\n", x$N, x$Time, x$M))
  cat(sprintf("Models evaluated: %d (all non-empty subsets)\n", x$nmodels))
  cat(sprintf("MCMC per model: %d draws, %d burn-in\n", x$ndraw, x$nomit))
  cat(sprintf("Total time: %.1f seconds\n\n", x$time))

  # Model probabilities table
  cat("Model probabilities:\n")
  model_table <- data.frame(
    LogMarginal = round(x$logm, digits),
    Prob = round(x$probs, digits),
    Rho = round(x$rho_means, digits)
  )
  # Add gamma columns
  for (m in seq_len(x$M)) {
    model_table[[paste0("W", m)]] <- round(x$gamma_all[, m], digits)
  }
  rownames(model_table) <- rownames(x$subsets)
  print(model_table)

  # BMA row
  cat(sprintf("\nBMA rho:   %.4f\n", x$rho))
  cat(sprintf("BMA sigma: %.4f\n", x$sige))
  cat("BMA gamma: ", paste(round(x$gamma, digits), collapse = ", "), "\n")

  # BMA coefficients
  cat("\nBMA posterior coefficients:\n")
  bsave <- as.matrix(x$bdraw)
  se <- apply(bsave, 2, sd)
  lower <- apply(bsave, 2, quantile, probs = 0.025)
  upper <- apply(bsave, 2, quantile, probs = 0.975)
  coef_table <- data.frame(
    Mean = round(x$beta, digits),
    StdDev = round(se, digits),
    Lower = round(lower, digits),
    Upper = round(upper, digits)
  )
  rownames(coef_table) <- paste0("X", seq_along(x$beta))
  colnames(coef_table) <- c("Mean", "Std.Dev", "2.5%", "97.5%")
  print(coef_table)

  # BMA effects (if available)
  if (!is.null(x$direct)) {
    cat("\nBMA Direct effects:\n")
    .print_effects_draws(x$direct, digits)
    cat("BMA Indirect effects:\n")
    .print_effects_draws(x$indirect, digits)
    cat("BMA Total effects:\n")
    .print_effects_draws(x$total, digits)
  }

  invisible(x)
}


# Helper to print effects from draws matrix
.print_effects_draws <- function(draws, digits = 4) {
  if (is.null(draws)) return(invisible(NULL))
  p <- ncol(draws)
  out <- data.frame(
    Mean = round(colMeans(draws), digits),
    StdDev = round(apply(draws, 2, sd), digits),
    Lower = round(apply(draws, 2, quantile, 0.025), digits),
    Upper = round(apply(draws, 2, quantile, 0.975), digits)
  )
  rownames(out) <- paste0("X", seq_len(p))
  colnames(out) <- c("Mean", "Std.Dev", "2.5%", "97.5%")
  print(out)
  cat("\n")
}


#' Summary Method for spmixW_bma Objects
#'
#' @param object An object of class \code{"spmixW_bma"}.
#' @param ... Further arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @method summary spmixW_bma
#' @export
summary.spmixW_bma <- function(object, ...) {
  print.spmixW_bma(object, ...)
  invisible(object)
}


#' Coefficient Extractor for spmixW_bma Objects
#'
#' @param object An object of class \code{"spmixW_bma"}.
#' @param ... Further arguments (ignored).
#' @return Named numeric vector of BMA-averaged posterior means.
#' @export
coef.spmixW_bma <- function(object, ...) {
  b <- object$beta
  names(b) <- paste0("X", seq_along(b))
  c(b, rho = object$rho)
}


#' Plot Method for spmixW_bma Objects
#'
#' Produces a multi-panel plot showing model probabilities and BMA
#' posterior densities.
#'
#' @param x An object of class \code{"spmixW_bma"}.
#' @param ... Further arguments (ignored).
#' @return Invisible NULL.
#' @export
plot.spmixW_bma <- function(x, ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  # 1. Model probability bar chart
  barplot(x$probs, names.arg = paste0("M", seq_along(x$probs)),
          main = "Posterior Model Probabilities",
          ylab = "Probability", col = "steelblue", las = 2, cex.names = 0.7)

  # 2. BMA rho posterior density
  rho_draws <- as.numeric(x$pdraw)
  d <- density(rho_draws)
  plot(d, main = expression("BMA posterior: " * rho),
       xlab = expression(rho), ylab = "Density", col = "darkred", lwd = 2)
  abline(v = x$rho, lty = 2, col = "gray40")

  # 3. BMA gamma posterior density (first component)
  if (x$M >= 2) {
    g1_draws <- as.matrix(x$gdraw)[, 1]
    d <- density(g1_draws, from = 0, to = 1)
    plot(d, main = expression("BMA posterior: " * gamma[1]),
         xlab = expression(gamma[1]), ylab = "Density", col = "darkgreen", lwd = 2)
    abline(v = x$gamma[1], lty = 2, col = "gray40")
  } else {
    plot.new()
  }

  # 4. BMA coefficient densities (first variable)
  b1_draws <- as.matrix(x$bdraw)[, 1]
  d <- density(b1_draws)
  plot(d, main = expression("BMA posterior: " * beta[1]),
       xlab = expression(beta[1]), ylab = "Density", col = "navy", lwd = 2)
  abline(v = x$beta[1], lty = 2, col = "gray40")

  invisible(NULL)
}
