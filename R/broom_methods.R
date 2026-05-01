#' Tidy Method for spmixW Objects
#'
#' Returns a data frame of estimated parameters with standard errors,
#' test statistics, p-values, and credible intervals.
#'
#' @param x An object of class \code{"spmixW"}.
#' @param effects Logical: if \code{TRUE} (default), include direct, indirect,
#'   and total effects rows in addition to raw coefficients.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with columns: \code{term}, \code{estimate},
#'   \code{std.error}, \code{statistic}, \code{p.value}, \code{conf.low},
#'   \code{conf.high}, \code{type}.
#'
#' @export
tidy.spmixW <- function(x, effects = TRUE, ...) {

  bsave <- as.matrix(x$bdraw)
  nvar <- ncol(bsave)

  # Variable names (uses same helper as print method)
  vnames <- .get_coef_names(x)

  # Coefficients
  rows <- list()
  for (j in seq_len(nvar)) {
    m <- x$beta[j]; s <- sd(bsave[, j])
    tval <- m / s
    pval <- 2 * (1 - pt(abs(tval), df = nrow(bsave) - 1))
    ci <- quantile(bsave[, j], c(0.025, 0.975))
    rows[[length(rows) + 1]] <- data.frame(
      term = vnames[j], estimate = m, std.error = s,
      statistic = tval, p.value = pval,
      conf.low = ci[1], conf.high = ci[2],
      type = "coefficient", stringsAsFactors = FALSE
    )
  }

  # Spatial parameter
  if (!is.null(x$rho)) {
    pd <- as.numeric(x$pdraw)
    m <- x$rho; s <- sd(pd)
    tval <- m / s
    ci <- quantile(pd, c(0.025, 0.975))
    pname <- if (.is_sem_type(x$meth %||% "")) "lambda" else "rho"
    rows[[length(rows) + 1]] <- data.frame(
      term = pname, estimate = m, std.error = s,
      statistic = tval, p.value = 2 * (1 - pt(abs(tval), df = length(pd) - 1)),
      conf.low = ci[1], conf.high = ci[2],
      type = "spatial", stringsAsFactors = FALSE
    )
  }

  # Sigma
  sd_draws <- as.numeric(x$sdraw)
  rows[[length(rows) + 1]] <- data.frame(
    term = "sigma2", estimate = x$sige, std.error = sd(sd_draws),
    statistic = NA_real_, p.value = NA_real_,
    conf.low = quantile(sd_draws, 0.025), conf.high = quantile(sd_draws, 0.975),
    type = "variance", stringsAsFactors = FALSE
  )

  # Gamma (convex weights)
  if (!is.null(x$gamma) && !is.null(x$gdraw)) {
    gd <- as.matrix(x$gdraw)
    for (m_idx in seq_along(x$gamma)) {
      gm <- x$gamma[m_idx]; gs <- sd(gd[, m_idx])
      tval <- gm / gs
      ci <- quantile(gd[, m_idx], c(0.025, 0.975))
      rows[[length(rows) + 1]] <- data.frame(
        term = paste0("gamma", m_idx), estimate = gm, std.error = gs,
        statistic = tval, p.value = 2 * (1 - pt(abs(tval), df = nrow(gd) - 1)),
        conf.low = ci[1], conf.high = ci[2],
        type = "weight", stringsAsFactors = FALSE
      )
    }
  }

  # Effects
  if (effects && !is.null(x$direct)) {
    p <- nrow(x$direct)
    enames <- .get_effect_names(x)

    for (j in seq_len(p)) {
      for (etype in c("direct", "indirect", "total")) {
        draws <- switch(etype,
          "direct" = x$direct_draws[, j],
          "indirect" = x$indirect_draws[, j],
          "total" = x$total_draws[, j])
        eff_mat <- switch(etype, "direct" = x$direct, "indirect" = x$indirect, "total" = x$total)
        m <- eff_mat[j, 1]; s <- sd(draws)
        tval <- m / s
        ci <- quantile(draws, c(0.025, 0.975))
        rows[[length(rows) + 1]] <- data.frame(
          term = enames[j], estimate = m, std.error = s,
          statistic = tval, p.value = 2 * (1 - pt(abs(tval), df = length(draws) - 1)),
          conf.low = ci[1], conf.high = ci[2],
          type = etype, stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Glance Method for spmixW Objects
#'
#' Returns a one-row data frame of model-level summary statistics.
#'
#' @param x An object of class \code{"spmixW"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A one-row data frame.
#'
#' @export
glance.spmixW <- function(x, ...) {

  out <- data.frame(
    r.squared  = x$rsqr,
    sigma2     = x$sige,
    nobs       = x$nobs,
    nregions   = x$N,
    nperiods   = x$Time,
    ndraw      = x$ndraw,
    nomit      = x$nomit,
    stringsAsFactors = FALSE
  )

  if (!is.null(x$rho)) out$rho <- x$rho
  if (!is.null(x$lik)) out$log_likelihood <- x$lik
  if (!is.null(x$logmarginal)) out$log_marginal <- x$logmarginal
  if (!is.null(x$rho_acc_rate)) out$acceptance_rho <- x$rho_acc_rate
  if (!is.null(x$gamma_acc_rate)) out$acceptance_gamma <- x$gamma_acc_rate

  out
}


#' Tidy Method for spmixW_bma Objects
#'
#' @param x An object of class \code{"spmixW_bma"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with BMA-averaged coefficients and model probabilities.
#'
#' @export
tidy.spmixW_bma <- function(x, ...) {

  bsave <- as.matrix(x$bdraw)
  nvar <- ncol(bsave)
  vnames <- x$predictor_names
  if (is.null(vnames)) vnames <- paste0("x", seq_len(nvar))
  if (length(vnames) < nvar) {
    vnames <- c(vnames, paste0("Wx", seq_len(nvar - length(vnames))))
  }

  rows <- list()

  # BMA-averaged coefficients
  for (j in seq_len(nvar)) {
    m <- x$beta[j]; s <- sd(bsave[, j])
    ci <- quantile(bsave[, j], c(0.025, 0.975))
    rows[[length(rows) + 1]] <- data.frame(
      term = vnames[j], estimate = m, std.error = s,
      conf.low = ci[1], conf.high = ci[2],
      type = "bma_coefficient", model_id = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  # Model probabilities
  for (i in seq_len(x$nmodels)) {
    subset_str <- paste(which(x$subsets[i, ] == 1), collapse = "+")
    rows[[length(rows) + 1]] <- data.frame(
      term = paste0("W{", subset_str, "}"),
      estimate = x$probs[i], std.error = NA_real_,
      conf.low = NA_real_, conf.high = NA_real_,
      type = "model_probability", model_id = i,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' Glance Method for spmixW_bma Objects
#'
#' @param x An object of class \code{"spmixW_bma"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A one-row data frame.
#'
#' @export
glance.spmixW_bma <- function(x, ...) {
  data.frame(
    nmodels       = x$nmodels,
    top_prob      = max(x$probs),
    top_model     = which.max(x$probs),
    bma_rho       = x$rho,
    bma_sigma2    = x$sige,
    nobs          = x$nobs,
    M             = x$M,
    ndraw         = x$ndraw,
    nomit         = x$nomit,
    stringsAsFactors = FALSE
  )
}
