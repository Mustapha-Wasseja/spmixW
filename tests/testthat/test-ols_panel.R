# tests/testthat/test-ols_panel.R

test_that("ols_panel recovers true beta from known DGP (homoscedastic, no FE)", {
  set.seed(12345)
  N <- 50; Time <- 10; k <- 2
  beta_true <- c(1.0, 2.0)
  sige_true <- 1.0

  X <- matrix(rnorm(N * Time * k), ncol = k)
  y <- as.numeric(X %*% beta_true + rnorm(N * Time, sd = sqrt(sige_true)))

  res <- ols_panel(y, X, N, Time, ndraw = 6000, nomit = 2000,
                   prior = list(model = 0, rval = 0))

  # Posterior mean should recover truth within +/- 0.15
  expect_equal(res$beta[1], beta_true[1], tolerance = 0.15)
  expect_equal(res$beta[2], beta_true[2], tolerance = 0.15)

  # Sigma should be in the right ballpark
  expect_equal(res$sige, sige_true, tolerance = 0.5)
})


test_that("ols_panel recovers true beta with region FE (model=1)", {
  set.seed(54321)
  N <- 50; Time <- 10; k <- 2
  beta_true <- c(1.0, 2.0)

  # Generate region fixed effects
  mu <- rnorm(N, sd = 2)
  X <- matrix(rnorm(N * Time * k), ncol = k)
  y <- as.numeric(X %*% beta_true) +
    rep(mu, times = Time) +
    rnorm(N * Time, sd = 1)

  res <- ols_panel(y, X, N, Time, ndraw = 6000, nomit = 2000,
                   prior = list(model = 1, rval = 0))

  expect_equal(res$beta[1], beta_true[1], tolerance = 0.15)
  expect_equal(res$beta[2], beta_true[2], tolerance = 0.15)
})


test_that("ols_panel recovers true beta with two-way FE (model=3)", {
  set.seed(11111)
  N <- 50; Time <- 10; k <- 2
  beta_true <- c(1.0, 2.0)

  # Generate both region and time fixed effects
  mu <- rnorm(N, sd = 2)    # region FE
  nu <- rnorm(Time, sd = 1) # time FE
  X <- matrix(rnorm(N * Time * k), ncol = k)
  y <- as.numeric(X %*% beta_true) +
    rep(mu, times = Time) +
    rep(nu, each = N) +
    rnorm(N * Time, sd = 1)

  res <- ols_panel(y, X, N, Time, ndraw = 6000, nomit = 2000,
                   prior = list(model = 3, rval = 0))

  expect_equal(res$beta[1], beta_true[1], tolerance = 0.15)
  expect_equal(res$beta[2], beta_true[2], tolerance = 0.15)
})


test_that("ols_panel heteroscedastic model still recovers beta", {
  set.seed(22222)
  N <- 50; Time <- 10; k <- 2
  beta_true <- c(1.0, 2.0)

  X <- matrix(rnorm(N * Time * k), ncol = k)
  y <- as.numeric(X %*% beta_true + rnorm(N * Time, sd = 1))

  res <- ols_panel(y, X, N, Time, ndraw = 6000, nomit = 2000,
                   prior = list(model = 0, rval = 4))

  expect_equal(res$beta[1], beta_true[1], tolerance = 0.15)
  expect_equal(res$beta[2], beta_true[2], tolerance = 0.15)
})


test_that("ols_panel returns correct S3 class and structure", {
  set.seed(33333)
  N <- 20; Time <- 5
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- ols_panel(y, X, N, Time, ndraw = 1000, nomit = 500,
                   prior = list(model = 1, rval = 0))

  expect_s3_class(res, "spmixW")
  expect_true(inherits(res$bdraw, "mcmc"))
  expect_true(inherits(res$sdraw, "mcmc"))
  expect_equal(nrow(as.matrix(res$bdraw)), 500)  # ndraw - nomit
  expect_equal(ncol(as.matrix(res$bdraw)), 2)
  expect_length(res$beta, 2)
  expect_length(res$yhat, N * Time)
  expect_length(res$resid, N * Time)
  expect_true(is.numeric(res$rsqr))
  expect_true(is.numeric(res$corr2))
  expect_true(is.numeric(res$lik))
  expect_length(res$sfe, N)
  expect_length(res$tfe, 0)  # model=1 has no time FE (empty numeric)
})


test_that("ols_panel with model=3 returns both sfe and tfe", {
  set.seed(44444)
  N <- 20; Time <- 5
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- ols_panel(y, X, N, Time, ndraw = 1000, nomit = 500,
                   prior = list(model = 3, rval = 0))

  expect_length(res$sfe, N)
  expect_length(res$tfe, Time)
  expect_true(!is.null(res$intercept))
})


test_that("ols_panel thinning works correctly", {
  set.seed(55555)
  N <- 20; Time <- 5
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- ols_panel(y, X, N, Time, ndraw = 2000, nomit = 500,
                   prior = list(model = 0, rval = 0, thin = 3))

  # (2000 - 500) / 3 = 500 saved draws
  expect_equal(nrow(as.matrix(res$bdraw)), 500)
})


test_that("ols_panel rejects invalid inputs", {
  N <- 10; Time <- 5
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y_wrong <- rnorm(N * Time + 1)

  expect_error(ols_panel(y_wrong, X, N, Time))
  expect_error(ols_panel(rnorm(N * Time), X, N, Time, ndraw = 100, nomit = 200))
})
