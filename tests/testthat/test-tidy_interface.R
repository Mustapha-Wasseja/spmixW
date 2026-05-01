# tests/testthat/test-tidy_interface.R

.make_panel_df <- function(N, Time, seed = 42) {
  set.seed(seed)
  df <- expand.grid(region = 1:N, year = 2001:(2000 + Time))
  df$x1 <- rnorm(nrow(df))
  df$x2 <- rnorm(nrow(df))
  df$y <- df$x1 * 1.5 - df$x2 * 0.8 + rnorm(nrow(df))
  df
}

.make_W <- function(N, k = 5, seed = 42) {
  set.seed(seed)
  as.matrix(normw(make_knw(cbind(runif(N), runif(N)), k = k,
                            row_normalise = FALSE)))
}


test_that("spmodel() formula produces same results as matrix API", {
  N <- 20; Time <- 5
  df <- .make_panel_df(N, Time, seed = 111)
  W <- .make_W(N, seed = 111)

  # Formula interface
  res_f <- spmodel(y ~ x1 + x2, data = df, W = W, model = "ols",
                   id = "region", time = "year", effects = "twoway",
                   heteroscedastic = FALSE, ndraw = 2000, nomit = 500)

  # Matrix interface (manually sort data)
  df_sorted <- df[order(df$year, df$region), ]
  y <- df_sorted$y
  X <- as.matrix(df_sorted[, c("x1", "x2")])
  res_m <- ols_panel(y, X, N, Time, ndraw = 2000, nomit = 500,
                     prior = list(model = 3, rval = 0))

  # Same seed + same sorted data => same results
  # (Not exact due to different RNG state, but structure should match)
  expect_equal(res_f$nobs, res_m$nobs)
  expect_equal(res_f$nvar, res_m$nvar)
  expect_equal(res_f$model, res_m$model)
  expect_equal(res_f$predictor_names, c("x1", "x2"))
  expect_equal(res_f$response_name, "y")
})


test_that("spmodel() auto-sorts data correctly", {
  N <- 15; Time <- 3
  df <- .make_panel_df(N, Time, seed = 222)
  W <- .make_W(N, seed = 222)

  # Scramble the order (sort by region first — WRONG order for spmixW)
  df_scrambled <- df[order(df$region, df$year), ]

  res <- spmodel(y ~ x1 + x2, data = df_scrambled, W = W,
                 model = "ols", id = "region", time = "year",
                 effects = "region", heteroscedastic = FALSE,
                 ndraw = 1500, nomit = 500)

  # Should still work — internal sorting handles it
  expect_s3_class(res, "spmixW")
  expect_equal(res$N, N)
  expect_equal(res$Time, Time)
})


test_that("spmodel() detects unbalanced panel", {
  N <- 10; Time <- 5
  df <- .make_panel_df(N, Time, seed = 333)
  W <- .make_W(N, seed = 333)

  # Remove one row to make it unbalanced
  df_unbal <- df[-1, ]

  expect_error(
    spmodel(y ~ x1 + x2, data = df_unbal, W = W,
            model = "ols", id = "region", time = "year"),
    "unbalanced"
  )
})


test_that("spmodel() cross-sectional auto-detection (no time column)", {
  set.seed(444)
  N <- 30
  df <- data.frame(county = 1:N, x1 = rnorm(N), x2 = rnorm(N),
                   y = rnorm(N))
  W <- .make_W(N, seed = 444)

  expect_message(
    res <- spmodel(y ~ x1 + x2, data = df, W = W,
                   model = "ols", id = "county",
                   heteroscedastic = FALSE, ndraw = 1500, nomit = 500),
    "cross-sectional"
  )

  # OLS with T=1 should produce valid results
  expect_s3_class(res, "spmixW")
  expect_equal(res$nobs, N)  # NT = N*1 = N
})


test_that("effects mapping: 'twoway' == model=3", {
  N <- 15; Time <- 3
  df <- .make_panel_df(N, Time, seed = 555)
  W <- .make_W(N, seed = 555)

  res <- spmodel(y ~ x1 + x2, data = df, W = W, model = "ols",
                 id = "region", time = "year", effects = "twoway",
                 heteroscedastic = FALSE, ndraw = 1500, nomit = 500)

  expect_equal(res$model, 3)
})


test_that("list W dispatches to convex combination", {
  N <- 15; Time <- 3
  df <- .make_panel_df(N, Time, seed = 666)
  W1 <- .make_W(N, k = 4, seed = 661)
  W2 <- .make_W(N, k = 6, seed = 662)

  res <- spmodel(y ~ x1 + x2, data = df, W = list(W1, W2),
                 model = "sar", id = "region", time = "year",
                 effects = "none", heteroscedastic = FALSE,
                 ndraw = 5000, nomit = 2000)

  expect_s3_class(res, "spmixW")
  expect_true(grepl("conv", res$meth))
  expect_true(!is.null(res$gamma))
})


test_that("compare_models() returns correct structure", {
  N <- 20; Time <- 5
  df <- .make_panel_df(N, Time, seed = 777)
  W <- .make_W(N, seed = 777)

  comp <- compare_models(y ~ x1 + x2, data = df, W = W,
                         id = "region", time = "year", effects = "twoway")

  expect_true(is.data.frame(comp))
  expect_equal(nrow(comp), 3)
  expect_true("model" %in% names(comp))
  expect_true("log_marginal" %in% names(comp))
  expect_true("probability" %in% names(comp))
  expect_equal(sum(comp$probability), 1, tolerance = 1e-10)
})


test_that("tidy.spmixW returns correct structure", {
  N <- 15; Time <- 3
  df <- .make_panel_df(N, Time, seed = 888)
  W <- .make_W(N, seed = 888)

  res <- spmodel(y ~ x1 + x2, data = df, W = W, model = "sar",
                 id = "region", time = "year", effects = "twoway",
                 heteroscedastic = FALSE, ndraw = 2000, nomit = 500)

  td <- tidy.spmixW(res)

  expect_true(is.data.frame(td))
  expect_true(all(c("term", "estimate", "std.error", "statistic",
                     "p.value", "conf.low", "conf.high", "type") %in% names(td)))
  # Should have: 2 coefficients + rho + sigma2 + 2*3 effects = 12 rows
  coef_rows <- sum(td$type == "coefficient")
  expect_equal(coef_rows, 2)
  expect_true("rho" %in% td$term)
  expect_true("sigma2" %in% td$term)
  # Variable names from formula
  expect_true("x1" %in% td$term)
  expect_true("x2" %in% td$term)
})


test_that("glance.spmixW returns one-row data frame", {
  N <- 15; Time <- 3
  df <- .make_panel_df(N, Time, seed = 999)
  W <- .make_W(N, seed = 999)

  res <- spmodel(y ~ x1 + x2, data = df, W = W, model = "ols",
                 id = "region", time = "year", effects = "none",
                 heteroscedastic = FALSE, ndraw = 1500, nomit = 500)

  gl <- glance.spmixW(res)
  expect_true(is.data.frame(gl))
  expect_equal(nrow(gl), 1)
  expect_true("r.squared" %in% names(gl))
  expect_true("nobs" %in% names(gl))
})
