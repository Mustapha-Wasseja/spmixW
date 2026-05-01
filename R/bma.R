#' Bayesian Model Averaging for SAR Convex Combination Panel
#'
#' Estimates SAR convex combination models for all \eqn{2^M - 1} non-empty
#' subsets of M weight matrices and averages estimates by posterior model
#' probabilities computed from log-marginal likelihoods.
#'
#' @param y Numeric vector of length NT.
#' @param X Numeric matrix NT x k.
#' @param Wlist A list of M spatial weight matrices (each N x N).
#' @param N Integer: number of cross-sectional units.
#' @param Time Integer: number of time periods.
#' @param ndraw Integer: MCMC draws per model (recommend >= 10000).
#' @param nomit Integer: burn-in draws per model.
#' @param prior List of prior hyperparameters (passed to each model).
#'
#' @return An S3 object of class \code{"spmixW_bma"} containing:
#'   \describe{
#'     \item{beta}{BMA posterior mean of beta.}
#'     \item{rho}{BMA posterior mean of rho.}
#'     \item{gamma}{BMA posterior mean of gamma (M x 1).}
#'     \item{sige}{BMA posterior mean of sigma^2.}
#'     \item{probs}{Vector of posterior model probabilities.}
#'     \item{logm}{Vector of log-marginal likelihoods per model.}
#'     \item{bdraw, pdraw, sdraw, gdraw}{BMA-averaged MCMC draws.}
#'     \item{direct, indirect, total}{BMA-averaged effects draws.}
#'     \item{subsets}{Matrix showing which W's are in each model.}
#'     \item{nmodels}{Number of models evaluated.}
#'     \item{per_model}{List of per-model results summaries.}
#'   }
#'
#' @references
#' Debarsy, N. and LeSage, J.P. (2021). "Bayesian model averaging for spatial
#' autoregressive models based on convex combinations of different types of
#' connectivity matrices." \emph{Journal of Business & Economic Statistics},
#' 40(2), 547-558.
#'
#' @export
sar_conv_bma <- function(y, X, Wlist, N, Time, ndraw = 25000L,
                         nomit = 5000L, prior = list()) {

  .run_bma(y, X, Wlist, N, Time, ndraw, nomit, prior,
           estimator = "sar", meth = "sar_conv_panel_bma_g")
}


#' Bayesian Model Averaging for SDM Convex Combination Panel
#'
#' @inheritParams sar_conv_bma
#' @return An S3 object of class \code{"spmixW_bma"}.
#' @export
sdm_conv_bma <- function(y, X, Wlist, N, Time, ndraw = 25000L,
                         nomit = 5000L, prior = list()) {

  .run_bma(y, X, Wlist, N, Time, ndraw, nomit, prior,
           estimator = "sdm", meth = "sdm_conv_panel_bma_g")
}


#' Bayesian Model Averaging for SDEM Convex Combination Panel
#'
#' @inheritParams sar_conv_bma
#' @return An S3 object of class \code{"spmixW_bma"}.
#' @export
sdem_conv_bma <- function(y, X, Wlist, N, Time, ndraw = 25000L,
                          nomit = 5000L, prior = list()) {

  .run_bma(y, X, Wlist, N, Time, ndraw, nomit, prior,
           estimator = "sdem", meth = "sdem_conv_panel_bma_g")
}


# ============================================================
# Internal: BMA engine shared by all three functions
# ============================================================
.run_bma <- function(y, X, Wlist, N, Time, ndraw, nomit, prior,
                     estimator, meth) {

  start_time <- proc.time()

  y <- as.numeric(y)
  X <- as.matrix(X)
  nobs <- N * Time
  nvar <- ncol(X)
  M <- length(Wlist)

  stopifnot(M >= 2, length(y) == nobs, nrow(X) == nobs)

  # ---- Enumerate all non-empty subsets ----
  # Generate all 2^M - 1 non-empty subsets
  all_subsets <- list()
  for (size in seq_len(M)) {
    combos <- combn(M, size, simplify = FALSE)
    all_subsets <- c(all_subsets, combos)
  }
  nmodels <- length(all_subsets)

  # Build subset indicator matrix for display
  subset_mat <- matrix(0L, nmodels, M)
  for (i in seq_len(nmodels)) {
    subset_mat[i, all_subsets[[i]]] <- 1L
  }
  colnames(subset_mat) <- paste0("W", seq_len(M))
  rownames(subset_mat) <- paste0("Model", seq_len(nmodels))

  # Choose estimator function
  est_fn <- switch(estimator,
    "sar"  = sar_conv_panel,
    "sdm"  = sdm_conv_panel,
    "sdem" = sdem_conv_panel
  )

  # ---- Run each model ----
  results_list <- vector("list", nmodels)
  logm <- numeric(nmodels)

  for (i in seq_len(nmodels)) {
    idx <- all_subsets[[i]]
    Wsub <- Wlist[idx]
    nw <- length(idx)

    message(sprintf("  BMA model %d/%d: W = {%s} (%d matrices)",
                    i, nmodels, paste(idx, collapse = ","), nw))

    if (nw == 1) {
      # Single W: use non-convex estimator for efficiency
      if (estimator == "sar") {
        res_i <- sar_panel(y, X, Wsub[[1]], N, Time, ndraw, nomit, prior)
        # Approximate log-marginal from concentrated likelihood
        res_i$logmarginal <- res_i$lik
      } else if (estimator == "sdm") {
        res_i <- sdm_panel(y, X, Wsub[[1]], N, Time, ndraw, nomit, prior)
        res_i$logmarginal <- res_i$lik
      } else {
        res_i <- sdem_panel(y, X, Wsub[[1]], N, Time, ndraw, nomit, prior)
        res_i$logmarginal <- res_i$lik
      }
      # Wrap gamma as length-1 for consistency
      res_i$gdraw <- coda::mcmc(matrix(1, nrow = nrow(as.matrix(res_i$bdraw)), ncol = 1))
      res_i$gamma <- 1
    } else {
      res_i <- est_fn(y, X, Wsub, N, Time, ndraw, nomit, prior)
    }

    results_list[[i]] <- res_i
    logm[i] <- res_i$logmarginal %||% res_i$lik %||% 0
  }

  # ---- Compute model probabilities ----
  probs <- model_probs(logm)

  # ---- BMA averaging of draws ----
  # For SDM/SDEM, different subsets have different numbers of augmented columns.
  # We average only the first nvar columns (the original X coefficients).
  bdraw_ref <- as.matrix(results_list[[1]]$bdraw)
  n_save <- nrow(bdraw_ref)

  bma_b <- matrix(0, n_save, nvar)
  bma_r <- numeric(n_save)
  bma_s <- numeric(n_save)
  bma_g <- matrix(0, n_save, M)

  # Effects (if available)
  has_effects <- !is.null(results_list[[1]]$direct_draws)
  p <- results_list[[1]]$p %||% nvar
  bma_dir <- bma_ind <- bma_tot <- NULL
  if (has_effects) {
    bma_dir <- matrix(0, n_save, p)
    bma_ind <- matrix(0, n_save, p)
    bma_tot <- matrix(0, n_save, p)
  }

  rho_means <- numeric(nmodels)
  beta_means <- matrix(0, nmodels, nvar)
  gamma_all <- matrix(0, nmodels, M)

  for (i in seq_len(nmodels)) {
    res_i <- results_list[[i]]
    idx <- all_subsets[[i]]
    w_i <- probs[i]

    bi <- as.matrix(res_i$bdraw)
    pi_draws <- as.numeric(res_i$pdraw)
    si_draws <- as.numeric(res_i$sdraw)
    gi_draws <- as.matrix(res_i$gdraw)

    n_i <- min(n_save, nrow(bi))
    # Only average the first nvar columns of beta draws
    ncol_bi <- min(nvar, ncol(bi))
    bma_b[1:n_i, 1:ncol_bi] <- bma_b[1:n_i, 1:ncol_bi] +
      w_i * bi[1:n_i, 1:ncol_bi, drop = FALSE]
    bma_r[1:n_i] <- bma_r[1:n_i] + w_i * pi_draws[1:n_i]
    bma_s[1:n_i] <- bma_s[1:n_i] + w_i * si_draws[1:n_i]

    # Map gamma draws back to full M-length vector
    for (j in seq_along(idx)) {
      if (j <= ncol(gi_draws)) {
        bma_g[1:n_i, idx[j]] <- bma_g[1:n_i, idx[j]] + w_i * gi_draws[1:n_i, j]
      }
    }

    # Effects
    if (has_effects && !is.null(res_i$direct_draws)) {
      di <- res_i$direct_draws
      ii <- res_i$indirect_draws
      ti <- res_i$total_draws
      p_i <- min(p, ncol(di))
      bma_dir[1:n_i, 1:p_i] <- bma_dir[1:n_i, 1:p_i] + w_i * di[1:n_i, 1:p_i, drop = FALSE]
      bma_ind[1:n_i, 1:p_i] <- bma_ind[1:n_i, 1:p_i] + w_i * ii[1:n_i, 1:p_i, drop = FALSE]
      bma_tot[1:n_i, 1:p_i] <- bma_tot[1:n_i, 1:p_i] + w_i * ti[1:n_i, 1:p_i, drop = FALSE]
    }

    rho_means[i] <- res_i$rho
    beta_means[i, ] <- res_i$beta[1:min(nvar, length(res_i$beta))]
    for (j in seq_along(idx)) {
      if (j <= length(res_i$gamma)) {
        gamma_all[i, idx[j]] <- res_i$gamma[j]
      }
    }
  }

  # BMA point estimates
  bma_beta  <- as.numeric(probs %*% beta_means)
  bma_rho   <- sum(probs * rho_means)
  bma_sige  <- mean(bma_s)
  bma_gamma <- as.numeric(probs %*% gamma_all)
  bma_logm  <- sum(probs * logm)

  thin <- prior$thin %||% 1L

  result <- list(
    beta      = bma_beta,
    rho       = bma_rho,
    sige      = bma_sige,
    gamma     = bma_gamma,
    probs     = probs,
    logm      = logm,
    lmarginal = bma_logm,
    bdraw     = coda::mcmc(bma_b, start = nomit + 1, thin = thin),
    pdraw     = coda::mcmc(matrix(bma_r, ncol = 1), start = nomit + 1, thin = thin),
    sdraw     = coda::mcmc(matrix(bma_s, ncol = 1), start = nomit + 1, thin = thin),
    gdraw     = coda::mcmc(bma_g, start = nomit + 1, thin = thin),
    direct    = bma_dir,
    indirect  = bma_ind,
    total     = bma_tot,
    subsets   = subset_mat,
    nmodels   = nmodels,
    rho_means = rho_means,
    beta_means = beta_means,
    gamma_all = gamma_all,
    per_model = results_list,
    nobs      = nobs,
    nvar      = nvar,
    p         = p,
    N         = N,
    Time      = Time,
    M         = M,
    meth      = meth,
    ndraw     = ndraw,
    nomit     = nomit,
    time      = (proc.time() - start_time)[3]
  )

  class(result) <- "spmixW_bma"
  result
}
