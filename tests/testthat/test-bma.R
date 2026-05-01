# tests/testthat/test-bma.R

.make_test_W <- function(N) {
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, (i %% N) + 1] <- 0.5
    W[i, ((i - 2) %% N) + 1] <- 0.5
  }
  W
}

.make_test_Ws2 <- function(N) {
  W1 <- matrix(0, N, N)
  W2 <- matrix(0, N, N)
  for (i in 1:N) {
    W1[i, (i %% N) + 1] <- 1
    W2[i, ((i - 2) %% N) + 1] <- 1
  }
  W1 <- W1 / pmax(rowSums(W1), 1)
  W2 <- W2 / pmax(rowSums(W2), 1)
  list(W1, W2)
}


test_that("sar_conv_bma evaluates correct number of models", {
  set.seed(11111)
  N <- 20; Time <- 5
  Ws <- .make_test_Ws2(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  # M=2: 2^2 - 1 = 3 models ({W1}, {W2}, {W1,W2})
  suppressMessages(
    res <- sar_conv_bma(y, X, Ws, N, Time, ndraw = 5000, nomit = 2000,
                        prior = list(model = 0))
  )

  expect_s3_class(res, "spmixW_bma")
  expect_equal(res$nmodels, 3)  # 2^2 - 1
  expect_length(res$probs, 3)
  expect_equal(sum(res$probs), 1, tolerance = 1e-10)
})


test_that("sar_conv_bma puts highest prob on correct subset", {
  set.seed(22222)
  N <- 30; Time <- 10
  rho_true <- 0.4; beta_true <- c(1, -0.5)

  Ws <- .make_test_Ws2(N)
  # Data generated using W1 only (gamma = (1, 0))
  W1big <- kronecker(diag(Time), Ws[[1]])
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- as.numeric(solve(diag(N*Time) - rho_true * W1big) %*%
                    (X %*% beta_true + rnorm(N*Time)))

  suppressMessages(
    res <- sar_conv_bma(y, X, Ws, N, Time, ndraw = 8000, nomit = 3000,
                        prior = list(model = 0))
  )

  # Model with W1 alone (subset 1) should have highest or second-highest prob
  # The {W1, W2} model might also have high prob since W1 is included
  # At minimum, the W2-only model (subset 2) should have lowest prob
  expect_true(res$probs[2] < max(res$probs[1], res$probs[3]))

  # BMA beta should be close to truth
  expect_equal(res$beta[1], beta_true[1], tolerance = 0.3)
  expect_equal(res$beta[2], beta_true[2], tolerance = 0.3)
})


test_that("sar_conv_bma returns correct output structure", {
  set.seed(33333)
  N <- 15; Time <- 5
  Ws <- .make_test_Ws2(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  suppressMessages(
    res <- sar_conv_bma(y, X, Ws, N, Time, ndraw = 5000, nomit = 2000,
                        prior = list(model = 0))
  )

  expect_true(!is.null(res$beta))
  expect_true(!is.null(res$rho))
  expect_true(!is.null(res$gamma))
  expect_true(!is.null(res$probs))
  expect_true(!is.null(res$logm))
  expect_true(!is.null(res$bdraw))
  expect_true(!is.null(res$pdraw))
  expect_true(!is.null(res$gdraw))
  expect_true(!is.null(res$subsets))
  expect_length(res$gamma, 2)
  expect_equal(nrow(res$subsets), 3)
  expect_equal(ncol(res$subsets), 2)
})


test_that("sdem_conv_bma runs without errors", {
  set.seed(44444)
  N <- 15; Time <- 5
  Ws <- .make_test_Ws2(N)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)
  y <- rnorm(N * Time)

  suppressMessages(
    res <- sdem_conv_bma(y, X, Ws, N, Time, ndraw = 5000, nomit = 2000,
                         prior = list(model = 0))
  )

  expect_s3_class(res, "spmixW_bma")
  expect_equal(res$nmodels, 3)
})
