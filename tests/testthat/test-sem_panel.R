# tests/testthat/test-sem_panel.R

.make_test_W <- function(N) {
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, ((i - 2) %% N) + 1] <- 0.5
    W[i, (i %% N) + 1] <- 0.5
  }
  W
}

test_that("sem_panel recovers lambda and beta from known DGP", {
  set.seed(55556)
  N <- 50; Time <- 10
  lambda_true <- 0.5; beta_true <- c(1.0, 2.0)

  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  B <- diag(N * Time) - lambda_true * Wbig
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  u <- as.numeric(solve(B) %*% rnorm(N * Time))
  y <- as.numeric(X %*% beta_true + u)

  res <- sem_panel(y, X, W, N, Time, ndraw = 6000, nomit = 2000,
                   prior = list(model = 0, rval = 0))

  expect_equal(res$rho, lambda_true, tolerance = 0.1)
  expect_equal(res$beta[1], beta_true[1], tolerance = 0.15)
  expect_equal(res$beta[2], beta_true[2], tolerance = 0.15)
})


test_that("sem_panel returns correct S3 class", {
  set.seed(66667)
  N <- 20; Time <- 5
  W <- .make_test_W(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- sem_panel(y, X, W, N, Time, ndraw = 1500, nomit = 500,
                   prior = list(model = 0, rval = 0))

  expect_s3_class(res, "spmixW")
  expect_true(grepl("sem", res$meth))
  expect_true(!is.null(res$rho))
  expect_true(!is.null(res$pdraw))
})


test_that("sem_panel works with region FE", {
  set.seed(77778)
  N <- 30; Time <- 5
  lambda_true <- 0.4; beta_true <- c(1, -0.5)

  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  mu <- rnorm(N, sd = 2)  # region FE
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  u <- as.numeric(solve(diag(N*Time) - lambda_true * Wbig) %*% rnorm(N*Time))
  y <- as.numeric(X %*% beta_true) + rep(mu, times = Time) + u

  res <- sem_panel(y, X, W, N, Time, ndraw = 5000, nomit = 2000,
                   prior = list(model = 1, rval = 0))

  expect_equal(res$rho, lambda_true, tolerance = 0.2)
  expect_equal(res$beta[1], beta_true[1], tolerance = 0.2)
})
