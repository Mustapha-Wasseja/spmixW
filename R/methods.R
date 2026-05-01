#' Print Method for spmixW Objects
#'
#' @param x An object of class \code{"spmixW"}.
#' @param digits Integer: number of significant digits to print.
#' @param ... Further arguments (ignored).
#'
#' @return Invisibly returns \code{x}.
#' @export
print.spmixW <- function(x, digits = 4, ...) {

  model_name <- .readable_model_name(x$meth)

  cat("========================================\n")
  cat(model_name, "\n")
  cat("========================================\n")
  cat(sprintf("N = %d, T = %d, NT = %d\n", x$N, x$Time, x$nobs))
  cat(sprintf("MCMC: %d draws, %d burn-in\n", x$ndraw, x$nomit))
  if (!is.null(x$rval)) {
    if (x$rval > 0) {
      cat(sprintf("Heteroscedastic errors (rval = %d)\n", x$rval))
    } else {
      cat("Homoscedastic errors\n")
    }
  }
  cat(sprintf("Elapsed time: %.2f seconds\n\n", x$time))

  # Coefficient table with proper variable names
  nvar <- length(x$beta)
  vnames <- .get_coef_names(x)
  se <- apply(as.matrix(x$bdraw), 2, sd)
  lower <- apply(as.matrix(x$bdraw), 2, quantile, probs = 0.025)
  upper <- apply(as.matrix(x$bdraw), 2, quantile, probs = 0.975)

  cat("Posterior estimates:\n")
  coef_table <- data.frame(
    Mean   = round(x$beta, digits),
    StdDev = round(se, digits),
    tstat  = round(x$tstat[seq_len(nvar)], digits),
    Lower  = round(lower, digits),
    Upper  = round(upper, digits)
  )
  rownames(coef_table) <- vnames
  colnames(coef_table) <- c("Mean", "Std.Dev", "t-stat", "2.5%", "97.5%")
  print(coef_table)

  # Spatial parameter
  if (!is.null(x$rho)) {
    rho_draws <- as.numeric(x$pdraw)
    param_name <- if (.is_sem_type(x$meth)) "lambda" else "rho"
    cat(sprintf("\n%s:  Mean = %.4f, Std.Dev = %.4f, [%.4f, %.4f]\n",
                param_name, x$rho, sd(rho_draws),
                quantile(rho_draws, 0.025), quantile(rho_draws, 0.975)))
  }

  # Gamma weights for convex models
  if (!is.null(x$gamma) && length(x$gamma) > 1) {
    cat("\nGamma weights (convex combination):\n")
    wnames <- .get_w_names(x)
    for (m in seq_along(x$gamma)) {
      cat(sprintf("  %s: %.4f\n", wnames[m], x$gamma[m]))
    }
    if (!is.null(x$rho_acc_rate)) {
      cat(sprintf("\nMH acceptance: rho=%.1f%%, gamma=%.1f%%\n",
                  x$rho_acc_rate * 100, x$gamma_acc_rate * 100))
    }
  }

  cat(sprintf("\nsigma^2: %.4f\n", x$sige))
  cat(sprintf("R-squared: %.4f\n", x$rsqr))
  if (!is.null(x$corr2)) cat(sprintf("Corr-squared: %.4f\n", x$corr2))
  if (!is.null(x$lik)) cat(sprintf("Log-likelihood: %.4f\n", x$lik))

  # Effects with variable names
  if (!is.null(x$direct)) {
    enames <- .get_effect_names(x)
    cat("\nDirect effects:\n")
    de <- round(x$direct, digits); rownames(de) <- enames; print(de)
    cat("\nIndirect effects:\n")
    ie <- round(x$indirect, digits); rownames(ie) <- enames; print(ie)
    cat("\nTotal effects:\n")
    te <- round(x$total, digits); rownames(te) <- enames; print(te)
  }

  invisible(x)
}


#' Summary Method for spmixW Objects
#'
#' @param object An object of class \code{"spmixW"}.
#' @param ... Further arguments (ignored).
#'
#' @return Invisibly returns \code{object} with MCMC diagnostics appended.
#' @method summary spmixW
#' @export
summary.spmixW <- function(object, ...) {

  cat("MCMC Diagnostics:\n")
  cat("-----------------\n")

  ess <- coda::effectiveSize(object$bdraw)
  cat("Effective sample sizes (beta):\n")
  print(round(ess))

  gd <- coda::geweke.diag(object$bdraw)
  cat("\nGeweke z-scores (beta):\n")
  print(round(gd$z, 4))

  if (!is.null(object$pdraw)) {
    ess_rho <- coda::effectiveSize(object$pdraw)
    cat(sprintf("\nEffective sample size (rho): %.1f\n", ess_rho))
  }

  cat("\n")
  print.spmixW(object, ...)

  invisible(object)
}


#' Coefficient Extractor for spmixW Objects
#'
#' @param object An object of class \code{"spmixW"}.
#' @param ... Further arguments (ignored).
#'
#' @return Named numeric vector of posterior means.
#' @export
coef.spmixW <- function(object, ...) {
  b <- object$beta
  names(b) <- .get_coef_names(object)
  if (!is.null(object$rho)) {
    param_name <- if (.is_sem_type(object$meth)) "lambda" else "rho"
    b <- c(b, setNames(object$rho, param_name))
  }
  b
}


#' Log-Likelihood for spmixW Objects
#'
#' @param object An object of class \code{"spmixW"}.
#' @param ... Further arguments (ignored).
#'
#' @return The log-likelihood evaluated at the posterior mean.
#' @export
logLik.spmixW <- function(object, ...) {
  ll <- object$lik
  attr(ll, "df") <- object$nvar + 1
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  ll
}


#' Plot Method for spmixW Objects
#'
#' Produces a 2x2 diagnostic plot layout: trace and density for the spatial
#' parameter, trace for sigma^2, and an effects comparison (or gamma densities
#' for convex combination models).
#'
#' @param x An object of class \code{"spmixW"}.
#' @param ... Further arguments (ignored).
#' @return Invisible NULL.
#' @export
plot.spmixW <- function(x, ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  has_rho <- !is.null(x$pdraw)
  has_gamma <- !is.null(x$gdraw)

  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

  if (has_rho) {
    rho_draws <- as.numeric(x$pdraw)
    pname <- if (.is_sem_type(x$meth)) expression(lambda) else expression(rho)

    plot(rho_draws, type = "l", col = "gray50",
         main = paste("Trace:", if (.is_sem_type(x$meth)) "lambda" else "rho"),
         xlab = "Iteration", ylab = pname)
    abline(h = x$rho, col = "red", lwd = 2)

    d <- density(rho_draws)
    plot(d, main = paste("Density:", if (.is_sem_type(x$meth)) "lambda" else "rho"),
         xlab = pname, col = "darkred", lwd = 2)
    abline(v = x$rho, lty = 2, col = "gray40")
  } else {
    plot.new(); plot.new()
  }

  sige_draws <- as.numeric(x$sdraw)
  plot(sige_draws, type = "l", col = "gray50",
       main = expression("Trace: " * sigma^2), xlab = "Iteration",
       ylab = expression(sigma^2))
  abline(h = x$sige, col = "red", lwd = 2)

  if (has_gamma && ncol(as.matrix(x$gdraw)) >= 2) {
    gdraw_mat <- as.matrix(x$gdraw)
    M <- ncol(gdraw_mat)
    cols <- c("steelblue", "darkgreen", "darkorange", "purple", "brown")[seq_len(M)]
    wnames <- .get_w_names(x)
    d1 <- density(gdraw_mat[, 1], from = 0, to = 1)
    plot(d1, main = "Gamma posteriors",
         xlab = expression(gamma), ylab = "Density",
         col = cols[1], lwd = 2, xlim = c(0, 1),
         ylim = c(0, max(d1$y) * 1.5))
    for (m in 2:M) {
      dm <- density(gdraw_mat[, m], from = 0, to = 1)
      lines(dm, col = cols[m], lwd = 2)
    }
    legend("topright", wnames, col = cols, lwd = 2, cex = 0.7)
  } else if (!is.null(x$direct)) {
    p <- nrow(x$direct)
    if (p > 0) {
      enames <- .get_effect_names(x)
      means <- cbind(x$direct[, 1], x$indirect[, 1], x$total[, 1])
      barplot(t(means), beside = TRUE, names.arg = enames,
              col = c("steelblue", "darkorange", "gray40"),
              main = "Effects estimates", ylab = "Effect")
      legend("topright", c("Direct", "Indirect", "Total"),
             fill = c("steelblue", "darkorange", "gray40"), cex = 0.7)
    }
  } else {
    b1 <- as.matrix(x$bdraw)[, 1]
    d <- density(b1)
    plot(d, main = expression("Density: " * beta[1]),
         xlab = expression(beta[1]), col = "navy", lwd = 2)
    abline(v = x$beta[1], lty = 2)
  }

  invisible(NULL)
}


# ============================================================
# Internal helpers for variable naming
# ============================================================

# Get coefficient names (handles SDM/SDEM augmented X)
.get_coef_names <- function(x) {
  nvar <- length(x$beta)
  pnames <- x$predictor_names

  if (!is.null(pnames)) {
    k_orig <- length(pnames)
    if (nvar == k_orig) {
      return(pnames)
    } else if (nvar > k_orig) {
      # SDM/SDEM/SLX: extra columns are WX terms
      wx_names <- paste0("W-", pnames)
      # Could be multiple W blocks for conv SDM/SDEM
      n_wx <- nvar - k_orig
      if (n_wx == k_orig) {
        return(c(pnames, wx_names))
      } else {
        # Multiple W blocks
        M <- n_wx %/% k_orig
        wx_all <- character(0)
        for (m in seq_len(M)) {
          wn <- .get_w_names(x)
          prefix <- if (!is.null(wn) && m <= length(wn)) wn[m] else paste0("W", m)
          wx_all <- c(wx_all, paste0(prefix, "-", pnames))
        }
        return(c(pnames, wx_all[seq_len(n_wx)]))
      }
    }
  }

  # Fallback
  paste0("var", seq_len(nvar))
}

# Get effect variable names
.get_effect_names <- function(x) {
  p <- if (!is.null(x$p)) x$p else nrow(x$direct)
  # For SDM/SDEM/SLX where p == nvar_orig (original X vars, no augmentation)
  pnames <- x$predictor_names
  if (!is.null(pnames) && length(pnames) >= p) {
    return(pnames[seq_len(p)])
  }
  # For conv SDM/SDEM where p > length(predictor_names),
  # use full coefficient names (includes W-prefixed)
  all_names <- .get_coef_names(x)
  # Skip intercept if present, use first p non-intercept names
  cflag <- !is.null(x$cflag) && x$cflag
  if (cflag && length(all_names) > p) {
    return(all_names[2:(p + 1)])
  }
  if (length(all_names) >= p) return(all_names[seq_len(p)])
  paste0("x", seq_len(p))
}

# Get W matrix names
.get_w_names <- function(x) {
  if (!is.null(x$W_names)) return(x$W_names)
  M <- x$M %||% length(x$gamma)
  if (is.null(M) || M == 0) return(NULL)
  paste0("W", seq_len(M))
}

# Check if model is SEM-type (lambda vs rho)
.is_sem_type <- function(meth) {
  grepl("sem|sdem", meth, ignore.case = TRUE)
}

# Readable model name from internal meth code
.readable_model_name <- function(meth) {
  # Map internal codes to readable names
  base <- if (grepl("ols", meth)) "OLS"
    else if (grepl("sdm", meth)) "SDM"
    else if (grepl("sdem", meth)) "SDEM"
    else if (grepl("sar", meth)) "SAR"
    else if (grepl("sem", meth)) "SEM"
    else if (grepl("slx", meth)) "SLX"
    else meth

  conv <- if (grepl("conv", meth)) " Convex Combination" else ""
  bma <- if (grepl("bma", meth)) " BMA" else ""

  fe <- if (grepl("stfe|twoway", meth)) " (Two-way FE)"
    else if (grepl("sfe", meth) && !grepl("stfe", meth)) " (Region FE)"
    else if (grepl("tfe", meth) && !grepl("stfe", meth)) " (Time FE)"
    else if (grepl("^p", meth)) " (Pooled)"
    else ""

  paste0(base, conv, bma, " Panel", fe)
}
