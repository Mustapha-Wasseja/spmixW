# tests/testthat/test-log_det.R

test_that("log_det_exact returns correct format", {
  W <- matrix(c(0, 0.5, 0.5, 0,
                0.5, 0, 0, 0.5,
                0.5, 0, 0, 0.5,
                0, 0.5, 0.5, 0), 4, 4)

  detval <- log_det_exact(W)

  expect_true(is.matrix(detval))
  expect_equal(ncol(detval), 2)
  expect_equal(colnames(detval), c("rho", "lndet"))

  # Check grid endpoints
  expect_equal(unname(detval[1, "rho"]), -1)
  expect_equal(unname(detval[nrow(detval), "rho"]), 1)
})


test_that("log_det_exact gives log|I - rho*W| = 0 at rho = 0", {
  W <- matrix(c(0, 0.5, 0.5, 0,
                0.5, 0, 0, 0.5,
                0.5, 0, 0, 0.5,
                0, 0.5, 0.5, 0), 4, 4)

  detval <- log_det_exact(W, rmin = -0.01, rmax = 0.01, grid_step = 0.001)

  # At rho = 0: log|I| = 0
  idx0 <- which(abs(detval[, "rho"]) < 1e-10)
  expect_equal(unname(detval[idx0, "lndet"]), 0, tolerance = 1e-12)
})


test_that("log_det_exact matches manual computation for small W", {
  # 2x2 matrix: W = [0, 1; 1, 0] (rook on 2 units)
  # Eigenvalues: +1, -1
  # log|I - rho*W| = log(1 - rho) + log(1 + rho) = log(1 - rho^2)
  W <- matrix(c(0, 1, 1, 0), 2, 2)

  detval <- log_det_exact(W, rmin = -0.99, rmax = 0.99, grid_step = 0.01)

  rho_grid <- detval[, "rho"]
  expected <- log(1 - rho_grid^2)

  expect_equal(detval[, "lndet"], expected, tolerance = 1e-12)
})


test_that("log_det_mc approximates exact for moderate N", {
  # Build a row-normalised 20x20 ring lattice
  N <- 20
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, ((i - 2) %% N) + 1] <- 0.5
    W[i, (i %% N) + 1] <- 0.5
  }

  set.seed(42)
  exact <- log_det_exact(W, rmin = -0.9, rmax = 0.9, grid_step = 0.1)
  mc    <- log_det_mc(W, rmin = -0.9, rmax = 0.9, grid_step = 0.1,
                      order = 50, iter = 30)

  # MC should be within 5% of exact for moderate rho values
  # (accuracy degrades near |rho| = 1)
  mid <- abs(exact[, "rho"]) < 0.7
  max_err <- max(abs(exact[mid, "lndet"] - mc[mid, "lndet"]))
  expect_true(max_err < 0.5,
              info = sprintf("Max MC error = %.4f (should be < 0.5)", max_err))
})


test_that("log_det_mc output is compatible with approx() lookup", {
  N <- 10
  W <- matrix(0, N, N)
  for (i in 1:N) {
    W[i, (i %% N) + 1] <- 1
  }
  W <- W / rowSums(W)

  set.seed(123)
  detval <- log_det_mc(W)

  # Should work with approx() for arbitrary rho values
  rho_test <- 0.35
  result <- approx(detval[, "rho"], detval[, "lndet"], xout = rho_test)
  expect_true(is.finite(result$y))
  expect_true(result$y < 0)  # log-det should be negative for positive rho
})
