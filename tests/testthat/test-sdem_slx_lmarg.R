# tests/testthat/test-sdem_slx_lmarg.R

.make_test_W <- function(N) {
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, ((i - 2) %% N) + 1] <- 0.5
    W[i, (i %% N) + 1] <- 0.5
  }
  W
}

# ---- sdem_panel tests ----

test_that("sdem_panel recovers lambda from SDEM DGP", {
  set.seed(11113)
  N <- 50; Time <- 10
  lambda_true <- 0.5; beta_true <- c(1, -0.5); theta_true <- c(0.3, 0.2)

  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  WX <- Wbig %*% X
  u <- as.numeric(solve(diag(N*Time) - lambda_true * Wbig) %*% rnorm(N*Time))
  y <- as.numeric(X %*% beta_true + WX %*% theta_true + u)

  res <- sdem_panel(y, X, W, N, Time, ndraw = 6000, nomit = 2000,
                    prior = list(model = 0, rval = 0))

  expect_equal(res$rho, lambda_true, tolerance = 0.12)
})


test_that("sdem_panel effects: direct + indirect = total", {
  set.seed(22224)
  N <- 30; Time <- 5
  W <- .make_test_W(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- sdem_panel(y, X, W, N, Time, ndraw = 3000, nomit = 1000,
                    prior = list(model = 0, rval = 0))

  for (j in seq_len(res$p)) {
    expect_equal(res$direct_draws[, j] + res$indirect_draws[, j],
                 res$total_draws[, j], tolerance = 1e-10)
  }
})


# ---- slx_panel tests ----

test_that("slx_panel recovers beta and theta from known DGP", {
  set.seed(33335)
  N <- 50; Time <- 10
  beta_true <- c(1, -0.5); theta_true <- c(0.3, 0.2)

  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  WX <- Wbig %*% X
  y <- as.numeric(X %*% beta_true + WX %*% theta_true + rnorm(N * Time, sd = 0.5))

  res <- slx_panel(y, X, W, N, Time, ndraw = 5000, nomit = 2000,
                   prior = list(model = 0, rval = 0))

  # Direct effects should recover beta
  expect_equal(unname(res$direct[1, "Mean"]), beta_true[1], tolerance = 0.15)
  expect_equal(unname(res$direct[2, "Mean"]), beta_true[2], tolerance = 0.15)
  # Indirect effects should recover theta
  expect_equal(unname(res$indirect[1, "Mean"]), theta_true[1], tolerance = 0.15)
  expect_equal(unname(res$indirect[2, "Mean"]), theta_true[2], tolerance = 0.2)
})


test_that("slx_panel effects: direct + indirect = total", {
  set.seed(44446)
  N <- 20; Time <- 5
  W <- .make_test_W(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- slx_panel(y, X, W, N, Time, ndraw = 2000, nomit = 500,
                   prior = list(model = 0, rval = 0))

  for (j in seq_len(res$p)) {
    expect_equal(res$direct_draws[, j] + res$indirect_draws[, j],
                 res$total_draws[, j], tolerance = 1e-10)
  }
})


# ---- lmarginal_panel tests ----

test_that("lmarginal_panel returns valid probabilities summing to 1", {
  set.seed(55557)
  N <- 30; Time <- 5
  W <- .make_test_W(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- lmarginal_panel(y, X, W, N, Time)

  expect_length(res$lmarginal, 3)
  expect_named(res$lmarginal, c("SLX", "SDM", "SDEM"))
  expect_equal(sum(res$probs), 1, tolerance = 1e-10)
  expect_true(all(res$probs >= 0))
  expect_true(all(is.finite(res$lmarginal)))
})


test_that("lmarginal_panel favours SDM for SAR-generated data", {
  set.seed(66668)
  N <- 40; Time <- 10
  rho_true <- 0.5
  W <- .make_test_W(N)
  Wbig <- kronecker(diag(Time), W)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- as.numeric(solve(diag(N*Time) - rho_true * Wbig) %*%
                    (X %*% c(1, -0.5) + rnorm(N * Time)))

  res <- lmarginal_panel(y, X, W, N, Time)

  # SDM should have higher marginal than SLX for spatially autocorrelated data
  expect_true(res$logm_sdm > res$logm_slx)
})
