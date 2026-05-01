# tests/testthat/test-model_probs.R

test_that("model_probs sums to 1", {
  lm <- c(-100, -98, -105, -102)
  probs <- model_probs(lm)

  expect_length(probs, 4)
  expect_equal(sum(probs), 1, tolerance = 1e-14)
})


test_that("model_probs assigns highest probability to highest log-marginal", {
  lm <- c(-100, -95, -110)
  probs <- model_probs(lm)

  expect_equal(which.max(probs), 2)
})


test_that("model_probs handles equal log-marginals", {
  lm <- c(-50, -50, -50)
  probs <- model_probs(lm)

  expect_equal(probs, rep(1/3, 3), tolerance = 1e-14)
})


test_that("model_probs handles single model", {
  probs <- model_probs(-100)
  expect_equal(probs, 1)
})


test_that("model_probs handles very different magnitudes without overflow", {
  # The scaling trick (subtract max) should prevent overflow
  lm <- c(-1000, -500, -10)
  probs <- model_probs(lm)

  expect_equal(sum(probs), 1, tolerance = 1e-14)
  # Model 3 dominates overwhelmingly

  expect_true(probs[3] > 0.999)
})
