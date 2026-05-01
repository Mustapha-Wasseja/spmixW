# tests/testthat/test-spatial_utils.R

test_that("normw produces row-stochastic matrix", {
  W <- matrix(c(0, 3, 1, 0,
                2, 0, 0, 4,
                1, 0, 0, 1,
                0, 5, 3, 0), 4, 4, byrow = TRUE)

  Wn <- normw(W)
  rs <- Matrix::rowSums(Wn)

  expect_equal(as.numeric(rs), rep(1, 4), tolerance = 1e-14)
})


test_that("normw handles zero rows (isolates) without NaN", {
  W <- matrix(c(0, 1, 0,
                1, 0, 0,
                0, 0, 0), 3, 3, byrow = TRUE)

  Wn <- normw(W)

  expect_true(all(is.finite(as.matrix(Wn))))
  expect_equal(as.numeric(Matrix::rowSums(Wn)), c(1, 1, 0))
})


test_that("normw returns sparse matrix", {
  W <- matrix(c(0, 1, 1, 0,
                1, 0, 0, 1,
                1, 0, 0, 1,
                0, 1, 1, 0), 4, 4)

  Wn <- normw(W)
  expect_s4_class(Wn, "dgCMatrix")
})


test_that("make_knw returns correct dimensions and row sums", {
  set.seed(1)
  coords <- cbind(runif(30), runif(30))
  k <- 5

  W <- make_knw(coords, k)

  expect_equal(dim(W), c(30, 30))
  # Row-normalised: all rows sum to 1
  expect_equal(as.numeric(Matrix::rowSums(W)), rep(1, 30), tolerance = 1e-14)
})


test_that("make_knw has exactly k non-zero entries per row before normalisation", {
  set.seed(2)
  coords <- cbind(runif(20), runif(20))
  k <- 4

  W <- make_knw(coords, k, row_normalise = FALSE)

  # Each row should have exactly k = 4 nonzero entries
  nnz_per_row <- Matrix::rowSums(W != 0)
  expect_equal(as.numeric(nnz_per_row), rep(k, 20))
})


test_that("make_knw diagonal is zero (no self-neighbours)", {
  set.seed(3)
  coords <- cbind(runif(15), runif(15))
  W <- make_knw(coords, k = 3)

  expect_equal(as.numeric(Matrix::diag(W)), rep(0, 15))
})
