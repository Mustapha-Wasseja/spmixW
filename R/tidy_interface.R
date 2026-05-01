#' Formula-Based Entry Point for Spatial Panel Models
#'
#' A user-friendly interface that accepts a formula and data frame, automatically
#' handles sorting, fixed-effects mapping, and dispatches to the appropriate
#' low-level estimation function.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...}. An intercept
#'   is NOT included in the design matrix (fixed effects handle the level).
#' @param data A data frame containing all variables referenced in the formula,
#'   plus the \code{id} and \code{time} columns.
#' @param W A spatial weight matrix (N x N) or a list of M weight matrices for
#'   convex combination models. When a list is provided, the convex combination
#'   variant of the specified model is used automatically.
#' @param model Character string specifying the model type. One of:
#'   \code{"ols"}, \code{"sar"}, \code{"sdm"}, \code{"sem"}, \code{"sdem"},
#'   \code{"slx"}, \code{"sar_bma"}, \code{"sdm_bma"}, \code{"sdem_bma"}.
#' @param id Character: name of the column in \code{data} identifying
#'   cross-sectional units (regions).
#' @param time Character or NULL: name of the column identifying time periods.
#'   If \code{NULL}, a cross-sectional model (T=1) is estimated with
#'   \code{effects = "none"}.
#' @param effects Character: fixed-effects specification. One of \code{"none"},
#'   \code{"region"}, \code{"time"}, \code{"twoway"}. Default \code{"twoway"}.
#' @param heteroscedastic Logical: if \code{TRUE} (default), use heteroscedastic
#'   errors with \code{rval = 5}. If \code{FALSE}, homoscedastic.
#' @param rval Numeric: degrees-of-freedom for heteroscedastic errors. Overrides
#'   \code{heteroscedastic} if provided. Default \code{NULL}.
#' @param ndraw Integer: total MCMC draws. Default 5500.
#' @param nomit Integer: burn-in draws. Default 1500.
#' @param thin Integer: thinning interval. Default 1.
#' @param taylor_order Integer: Taylor series order for convex combination
#'   models. Default 6.
#' @param ... Additional arguments passed to the prior list (e.g.,
#'   \code{beta_mean}, \code{beta_var}, \code{nu}, \code{d0}).
#'
#' @return An S3 object of class \code{"spmixW"} or \code{"spmixW_bma"} with
#'   additional metadata:
#'   \describe{
#'     \item{formula}{The formula used.}
#'     \item{response_name}{Name of the response variable.}
#'     \item{predictor_names}{Names of predictor variables.}
#'     \item{id_var, time_var}{Column names used for panel structure.}
#'   }
#'
#' @details
#' \strong{Auto-sorting:} The data is sorted by time first, then by region
#' within each time period, to match the package's required stacking order
#' (all N regions in period 1, then all N in period 2, etc.). Users do not
#' need to pre-sort their data.
#'
#' \strong{W list detection:} When \code{W} is a list of matrices, the convex
#' combination variant of the model is used automatically (e.g., \code{"sar"}
#' becomes \code{sar_conv_panel()}).
#'
#' \strong{Cross-sectional mode:} When \code{time = NULL}, the function
#' automatically sets T=1 and \code{effects = "none"}.
#'
#' @examples
#' \donttest{
#' coords <- cbind(runif(80), runif(80))
#' W_geo <- make_knw(coords, k = 5)
#' W_trade <- make_knw(coords, k = 10)
#'
#' panel <- simulate_panel(N = 80, T = 8,
#'                          W = list(W_geo, W_trade),
#'                          gamma = c(0.6, 0.4), rho = 0.5,
#'                          beta = c(1.5, -0.8), seed = 123)
#'
#' res <- spmodel(y ~ x1 + x2, data = panel,
#'                W = list(geography = W_geo, trade = W_trade),
#'                model = "sar", id = "region", time = "year",
#'                effects = "twoway", ndraw = 8000, nomit = 2000)
#' print(res)
#' }
#'
#' @export
spmodel <- function(formula, data, W, model = "sar",
                    id, time = NULL,
                    effects = "twoway",
                    heteroscedastic = TRUE,
                    rval = NULL,
                    ndraw = 5500L, nomit = 1500L, thin = 1L,
                    taylor_order = 6L,
                    ...) {

  # ---- 1. Parse formula ----
  mf <- model.frame(formula, data = data, na.action = na.fail)
  y_name <- names(mf)[1]
  x_names <- attr(terms(formula), "term.labels")

  if (length(x_names) == 0) stop("Formula must include at least one predictor.")

  # ---- 2. Determine N and T ----
  if (is.null(time)) {
    message("No 'time' column specified: assuming cross-sectional data (T=1, effects='none').")
    data$`.time_dummy` <- 1L
    time <- ".time_dummy"
    effects <- "none"
  }

  ids <- data[[id]]
  times <- data[[time]]
  N <- length(unique(ids))
  Time <- length(unique(times))

  if (N * Time != nrow(data)) {
    stop(sprintf("Panel is unbalanced: expected N x T = %d x %d = %d rows, got %d",
                 N, Time, N * Time, nrow(data)))
  }

  # ---- 3. Sort data: by time first, then by region within each period ----
  data <- data[order(data[[time]], data[[id]]), ]

  # ---- 4. Extract y and X ----
  y <- data[[y_name]]
  X <- as.matrix(data[, x_names, drop = FALSE])

  # ---- 5. Map effects to model code ----
  model_code <- switch(effects,
    "none" = 0L, "region" = 1L, "time" = 2L, "twoway" = 3L,
    stop(sprintf("Unknown effects: '%s'. Use 'none', 'region', 'time', or 'twoway'.", effects))
  )

  # ---- 6. Build prior list ----
  if (is.null(rval)) {
    rval <- if (heteroscedastic) 5 else 0
  }

  prior <- list(model = model_code, rval = rval, thin = thin,
                taylor_order = taylor_order)
  extra <- list(...)
  prior <- c(prior, extra)

  # ---- 7. Dispatch ----
  is_list_W <- is.list(W) && !inherits(W, "Matrix")
  is_bma <- grepl("_bma$", model)

  result <- switch(model,
    "ols" = {
      if (is_list_W) stop("OLS does not use W. Pass a single W or use a spatial model.")
      ols_panel(y, X, N, Time, ndraw, nomit, prior)
    },
    "sar" = {
      if (is_list_W) sar_conv_panel(y, X, W, N, Time, ndraw, nomit, prior)
      else sar_panel(y, X, W, N, Time, ndraw, nomit, prior)
    },
    "sdm" = {
      if (is_list_W) sdm_conv_panel(y, X, W, N, Time, ndraw, nomit, prior)
      else sdm_panel(y, X, W, N, Time, ndraw, nomit, prior)
    },
    "sem" = {
      if (is_list_W) sem_conv_panel(y, X, W, N, Time, ndraw, nomit, prior)
      else sem_panel(y, X, W, N, Time, ndraw, nomit, prior)
    },
    "sdem" = {
      if (is_list_W) sdem_conv_panel(y, X, W, N, Time, ndraw, nomit, prior)
      else sdem_panel(y, X, W, N, Time, ndraw, nomit, prior)
    },
    "slx" = {
      if (is_list_W) stop("SLX convex combination not supported. Use 'sdem' for spatial lags of X with convex W.")
      slx_panel(y, X, W, N, Time, ndraw, nomit, prior)
    },
    "sar_bma" = {
      if (!is_list_W) stop("BMA requires a list of W matrices.")
      sar_conv_bma(y, X, W, N, Time, ndraw, nomit, prior)
    },
    "sdm_bma" = {
      if (!is_list_W) stop("BMA requires a list of W matrices.")
      sdm_conv_bma(y, X, W, N, Time, ndraw, nomit, prior)
    },
    "sdem_bma" = {
      if (!is_list_W) stop("BMA requires a list of W matrices.")
      sdem_conv_bma(y, X, W, N, Time, ndraw, nomit, prior)
    },
    stop(sprintf("Unknown model: '%s'. Choose from: ols, sar, sdm, sem, sdem, slx, sar_bma, sdm_bma, sdem_bma.", model))
  )

  # ---- 8. Attach metadata ----
  result$formula <- formula
  result$response_name <- y_name
  result$predictor_names <- x_names
  result$id_var <- id
  result$time_var <- if (time == ".time_dummy") NULL else time

  # Store W matrix names for gamma display
  if (is_list_W && !is.null(names(W))) {
    result$W_names <- names(W)
  }

  result
}


#' Compare Spatial Panel Model Specifications
#'
#' Computes log-marginal likelihoods for SLX, SDM, and SDEM specifications
#' and returns posterior model probabilities. A convenience wrapper around
#' \code{\link{lmarginal_panel}}.
#'
#' @inheritParams spmodel
#'
#' @return A data frame with columns: \code{model}, \code{log_marginal},
#'   \code{probability}.
#'
#' @examples
#' \donttest{
#' coords <- cbind(runif(50), runif(50))
#' W <- make_knw(coords, k = 5)
#' panel <- simulate_panel(N = 50, T = 8, W = W, rho = 0.4,
#'                          beta = c(1, -0.5), seed = 42)
#' comp <- compare_models(y ~ x1 + x2, data = panel, W = as.matrix(W),
#'                        id = "region", time = "year", effects = "twoway")
#' print(comp)
#' }
#'
#' @export
compare_models <- function(formula, data, W,
                           id, time = NULL,
                           effects = "twoway", ...) {

  # Parse formula
  mf <- model.frame(formula, data = data, na.action = na.fail)
  y_name <- names(mf)[1]
  x_names <- attr(terms(formula), "term.labels")

  if (is.null(time)) {
    data$`.time_dummy` <- 1L
    time <- ".time_dummy"
    effects <- "none"
  }

  N <- length(unique(data[[id]]))
  Time <- length(unique(data[[time]]))

  if (N * Time != nrow(data)) {
    stop(sprintf("Panel is unbalanced: expected N x T = %d rows, got %d",
                 N * Time, nrow(data)))
  }

  # Sort
  data <- data[order(data[[time]], data[[id]]), ]
  y <- data[[y_name]]
  X <- as.matrix(data[, x_names, drop = FALSE])

  model_code <- switch(effects,
    "none" = 0L, "region" = 1L, "time" = 2L, "twoway" = 3L)

  # Demean
  dm <- demean_panel(y, X, N, Time, model_code)

  # Handle single W or list of W
  extra <- list(...)
  if (is.list(W) && !inherits(W, "Matrix")) {
    # Compare across individual W matrices
    rows <- list()
    for (i in seq_along(W)) {
      lm_i <- lmarginal_panel(dm$ywith, dm$xwith, W[[i]], N, Time, extra)
      for (m in seq_along(lm_i$lmarginal)) {
        rows[[length(rows) + 1]] <- data.frame(
          W = paste0("W", i),
          model = names(lm_i$lmarginal)[m],
          log_marginal = lm_i$lmarginal[m],
          stringsAsFactors = FALSE
        )
      }
    }
    out <- do.call(rbind, rows)
    # Compute global probabilities
    out$probability <- model_probs(out$log_marginal)
  } else {
    lm_res <- lmarginal_panel(dm$ywith, dm$xwith, W, N, Time, extra)
    out <- data.frame(
      model = names(lm_res$lmarginal),
      log_marginal = as.numeric(lm_res$lmarginal),
      probability = as.numeric(lm_res$probs),
      stringsAsFactors = FALSE
    )
  }

  rownames(out) <- NULL
  out
}
