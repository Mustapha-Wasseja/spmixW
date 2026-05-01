# tests/testthat/test-sdm_panel.R

.make_test_W <- function(N) {
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, ((i - 2) %% N) + 1] <- 0.5
    W[i, (i %% N) + 1] <- 0.5
  }
  W
}

test_that("sdm_panel recovers rho from known DGP", {
  set.seed(22223)
  N <- 50; Time <- 10
  rho_true <- 0.5; beta_true <- c(1, -0.5); theta_true <- c(0.3, 0.2)

  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  A <- diag(N * Time) - rho_true * Wbig
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  WX <- Wbig %*% X
  y <- as.numeric(solve(A) %*% (X %*% beta_true + WX %*% theta_true + rnorm(N * Time)))

  res <- sdm_panel(y, X, W, N, Time, ndraw = 6000, nomit = 2000,
                   prior = list(model = 0, rval = 0))

  expect_equal(res$rho, rho_true, tolerance = 0.12)
})


test_that("sdm_panel effects: direct + indirect = total", {
  set.seed(33334)
  N <- 30; Time <- 5
  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  WX <- Wbig %*% X
  y <- as.numeric(solve(diag(N*Time) - 0.4 * Wbig) %*%
                    (X %*% c(1, -0.5) + WX %*% c(0.2, 0.1) + rnorm(N*Time)))

  res <- sdm_panel(y, X, W, N, Time, ndraw = 4000, nomit = 1500,
                   prior = list(model = 0, rval = 0))

  for (j in seq_len(res$p)) {
    expect_equal(res$direct_draws[, j] + res$indirect_draws[, j],
                 res$total_draws[, j], tolerance = 1e-10)
  }
})


test_that("sdm_panel returns correct structure", {
  set.seed(44445)
  N <- 20; Time <- 5
  W <- .make_test_W(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- sdm_panel(y, X, W, N, Time, ndraw = 1500, nomit = 500,
                   prior = list(model = 0, rval = 0))

  expect_s3_class(res, "spmixW")
  expect_true(grepl("sdm", res$meth))
  expect_equal(res$p, 2)  # 2 non-intercept X variables
  expect_equal(nrow(res$direct), 2)
})
