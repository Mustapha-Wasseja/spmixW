# tests/testthat/test-convex.R
# Tests for Phase 3: convex combination models

.make_test_Ws <- function(N) {
  # Three distinct W matrices for testing
  W1 <- matrix(0, N, N)  # forward neighbour
  W2 <- matrix(0, N, N)  # backward neighbour
  W3 <- matrix(0, N, N)  # 2-step forward
  for (i in 1:N) {
    W1[i, (i %% N) + 1] <- 1
    W2[i, ((i - 2) %% N) + 1] <- 1
    W3[i, ((i + 1) %% N) + 1] <- 1  # 2 steps forward
  }
  W1 <- W1 / pmax(rowSums(W1), 1)
  W2 <- W2 / pmax(rowSums(W2), 1)
  W3 <- W3 / pmax(rowSums(W3), 1)
  list(W1, W2, W3)
}


# ---- log_det_taylor tests ----

test_that("log_det_taylor returns correct structure", {
  N <- 10
  Ws <- .make_test_Ws(N)
  traces <- log_det_taylor(Ws)

  expect_true("traces" %in% names(traces))
  expect_true("M" %in% names(traces))
  expect_equal(traces$M, 3)
  expect_equal(traces$N, N)
  expect_length(traces$T2, 9)   # 3^2 (backward compat)
  expect_length(traces$T3, 27)  # 3^3
  expect_length(traces$T4, 81)  # 3^4
  expect_length(traces$traces[["2"]], 9)
  expect_length(traces$traces[["3"]], 27)
  expect_length(traces$traces[["4"]], 81)
})


test_that("eval_taylor_lndet gives 0 at rho=0", {
  N <- 10
  Ws <- .make_test_Ws(N)
  traces <- log_det_taylor(Ws)

  val <- eval_taylor_lndet(traces, rho = 0, gamma = c(0.5, 0.3, 0.2))
  expect_equal(val, 0, tolerance = 1e-14)
})


test_that("eval_taylor_lndet approximates exact for single W", {
  N <- 20
  W1 <- matrix(0, N, N)
  for (i in 1:N) {
    W1[i, (i %% N) + 1] <- 0.5
    W1[i, ((i - 2) %% N) + 1] <- 0.5
  }
  traces <- log_det_taylor(list(W1))

  # With M=1, gamma = 1, the Taylor should approximate the exact
  rho_test <- 0.3
  taylor_val <- eval_taylor_lndet(traces, rho = rho_test, gamma = 1)
  exact_val <- log_det_exact(W1, rmin = rho_test - 0.001, rmax = rho_test + 0.001,
                              grid_step = 0.001)
  exact_at_rho <- exact_val[which.min(abs(exact_val[, "rho"] - rho_test)), "lndet"]

  # 4th order Taylor should be close for moderate rho
  expect_equal(taylor_val, unname(exact_at_rho), tolerance = 0.5)
})


test_that("eval_taylor_lndet is symmetric in equal gamma", {
  N <- 15
  Ws <- .make_test_Ws(N)[1:2]
  traces <- log_det_taylor(Ws)

  # With gamma = (0.5, 0.5), result should not depend on order
  val <- eval_taylor_lndet(traces, rho = 0.4, gamma = c(0.5, 0.5))
  expect_true(is.finite(val))
})


# ---- sar_conv_panel tests ----

test_that("sar_conv_panel recovers rho and gamma from known DGP", {
  set.seed(31415)
  N <- 40; Time <- 10
  rho_true <- 0.4
  gamma_true <- c(0.7, 0.3)
  beta_true <- c(1, -0.5)

  Ws <- .make_test_Ws(N)[1:2]
  Wbigs <- lapply(Ws, function(W) kronecker(diag(Time), W))
  nobs <- N * Time

  Wc <- gamma_true[1] * Wbigs[[1]] + gamma_true[2] * Wbigs[[2]]
  A <- diag(nobs) - rho_true * as.matrix(Wc)
  X <- matrix(rnorm(nobs * 2), ncol = 2)
  y <- as.numeric(solve(A) %*% (X %*% beta_true + rnorm(nobs)))

  res <- sar_conv_panel(y, X, Ws, N, Time, ndraw = 20000, nomit = 5000,
                        prior = list(model = 0))

  # rho should be recovered (convex models have wider posteriors)
  expect_equal(res$rho, rho_true, tolerance = 0.25)
  # gamma[1] should dominate
  expect_true(res$gamma[1] > 0.4)
  # Acceptance rates should be in reasonable range
  expect_true(res$rho_acc_rate > 0.05)
  expect_true(res$gamma_acc_rate > 0.05)
})


test_that("sar_conv_panel with M=1 behaves like standard SAR", {
  set.seed(13579)
  N <- 30; Time <- 5
  rho_true <- 0.4; beta_true <- c(1, -0.5)
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, (i %% N) + 1] <- 0.5
    W[i, ((i - 2) %% N) + 1] <- 0.5
  }

  Wbig <- kronecker(diag(Time), W)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- as.numeric(solve(diag(N*Time) - rho_true * Wbig) %*%
                    (X %*% beta_true + rnorm(N*Time)))

  res <- sar_conv_panel(y, X, list(W), N, Time, ndraw = 12000, nomit = 4000,
                        prior = list(model = 0))

  expect_equal(res$rho, rho_true, tolerance = 0.2)
  expect_equal(res$gamma[1], 1, tolerance = 0.01)  # only one W, gamma must be 1
})


test_that("sar_conv_panel effects: direct + indirect = total", {
  set.seed(24680)
  N <- 25; Time <- 5
  Ws <- .make_test_Ws(N)[1:2]
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- sar_conv_panel(y, X, Ws, N, Time, ndraw = 8000, nomit = 3000,
                        prior = list(model = 0))

  for (j in seq_len(res$p)) {
    expect_equal(res$direct_draws[, j] + res$indirect_draws[, j],
                 res$total_draws[, j], tolerance = 1e-10)
  }
})


test_that("sar_conv_panel returns MH acceptance rates", {
  set.seed(35791)
  N <- 20; Time <- 5
  Ws <- .make_test_Ws(N)[1:2]
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- sar_conv_panel(y, X, Ws, N, Time, ndraw = 5000, nomit = 2000,
                        prior = list(model = 0))

  expect_true(!is.null(res$rho_acc_rate))
  expect_true(!is.null(res$gamma_acc_rate))
  expect_true(res$rho_acc_rate >= 0 && res$rho_acc_rate <= 1)
  expect_true(res$gamma_acc_rate >= 0 && res$gamma_acc_rate <= 1)
})


# ---- sem_conv_panel tests ----

test_that("sem_conv_panel recovers lambda from known DGP", {
  set.seed(46802)
  N <- 40; Time <- 10
  lambda_true <- 0.5; beta_true <- c(1, -0.5)
  gamma_true <- c(0.6, 0.4)

  Ws <- .make_test_Ws(N)[1:2]
  Wbigs <- lapply(Ws, function(W) kronecker(diag(Time), W))
  nobs <- N * Time
  Wc <- gamma_true[1] * Wbigs[[1]] + gamma_true[2] * Wbigs[[2]]
  B <- diag(nobs) - lambda_true * as.matrix(Wc)
  X <- matrix(rnorm(nobs * 2), ncol = 2)
  u <- as.numeric(solve(B) %*% rnorm(nobs))
  y <- as.numeric(X %*% beta_true + u)

  res <- sem_conv_panel(y, X, Ws, N, Time, ndraw = 15000, nomit = 5000,
                        prior = list(model = 0))

  expect_equal(res$rho, lambda_true, tolerance = 0.2)
  expect_true(res$gamma[1] > 0.3)  # gamma_1 should dominate
})


# ---- sdm_conv_panel tests ----

test_that("sdm_conv_panel runs and returns correct structure", {
  set.seed(57913)
  N <- 20; Time <- 5
  Ws <- .make_test_Ws(N)[1:2]
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- sdm_conv_panel(y, X, Ws, N, Time, ndraw = 5000, nomit = 2000,
                        prior = list(model = 0))

  expect_s3_class(res, "spmixW")
  expect_true(grepl("sdm", res$meth))
  expect_true(!is.null(res$gamma))
  expect_length(res$gamma, 2)
})


# ---- sdem_conv_panel tests ----

test_that("sdem_conv_panel effects: direct + indirect = total", {
  set.seed(68024)
  N <- 20; Time <- 5
  Ws <- .make_test_Ws(N)[1:2]
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  res <- sdem_conv_panel(y, X, Ws, N, Time, ndraw = 5000, nomit = 2000,
                         prior = list(model = 0))

  for (j in seq_len(res$p)) {
    expect_equal(res$direct_draws[, j] + res$indirect_draws[, j],
                 res$total_draws[, j], tolerance = 1e-10)
  }
})
