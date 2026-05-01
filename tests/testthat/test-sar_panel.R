# tests/testthat/test-sar_panel.R

.make_test_W <- function(N) {
  # Simple ring lattice W (row-normalised)
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, ((i - 2) %% N) + 1] <- 0.5
    W[i, (i %% N) + 1] <- 0.5
  }
  W
}

test_that("sar_panel recovers rho and beta from known DGP", {
  set.seed(77777)
  N <- 50; Time <- 10
  rho_true <- 0.5; beta_true <- c(1.0, 2.0)

  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  A <- diag(N * Time) - rho_true * Wbig
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- as.numeric(solve(A) %*% (X %*% beta_true + rnorm(N * Time)))

  res <- sar_panel(y, X, W, N, Time, ndraw = 6000, nomit = 2000,
                   prior = list(model = 0, rval = 0))

  expect_equal(res$rho, rho_true, tolerance = 0.1)
  expect_equal(res$beta[1], beta_true[1], tolerance = 0.15)
  expect_equal(res$beta[2], beta_true[2], tolerance = 0.15)
})


test_that("sar_panel effects: direct + indirect = total", {
  set.seed(88888)
  N <- 30; Time <- 5
  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- as.numeric(solve(diag(N*Time) - 0.4 * Wbig) %*%
                    (X %*% c(1, -0.5) + rnorm(N*Time)))

  res <- sar_panel(y, X, W, N, Time, ndraw = 4000, nomit = 1500,
                   prior = list(model = 0, rval = 0))

  # For each variable, direct + indirect should equal total
  for (j in seq_len(res$p)) {
    d <- res$direct_draws[, j]
    ind <- res$indirect_draws[, j]
    tot <- res$total_draws[, j]
    expect_equal(d + ind, tot, tolerance = 1e-10)
  }
})


test_that("sar_panel returns correct S3 class and structure", {
  set.seed(99999)
  N <- 20; Time <- 5
  W <- .make_test_W(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- sar_panel(y, X, W, N, Time, ndraw = 1500, nomit = 500,
                   prior = list(model = 1, rval = 0))

  expect_s3_class(res, "spmixW")
  expect_true(inherits(res$bdraw, "mcmc"))
  expect_true(inherits(res$pdraw, "mcmc"))
  expect_true(!is.null(res$rho))
  expect_true(!is.null(res$direct))
  expect_true(!is.null(res$indirect))
  expect_true(!is.null(res$total))
  expect_true(!is.null(res$lndet))
  expect_equal(nrow(res$direct), res$p)
})


test_that("sar_panel works with heteroscedastic errors", {
  set.seed(42042)
  N <- 50; Time <- 10
  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- as.numeric(solve(diag(N*Time) - 0.3 * Wbig) %*%
                    (X %*% c(1, -0.5) + rnorm(N*Time)))

  res <- sar_panel(y, X, W, N, Time, ndraw = 6000, nomit = 2000,
                   prior = list(model = 0, rval = 4))

  expect_true(res$rval == 4)
  expect_equal(res$rho, 0.3, tolerance = 0.2)
})
