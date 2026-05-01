# tests/testthat/test-demean_panel.R
# Hand-computable example: N=4 regions, T=3 time periods

test_that("demean_panel model=0 returns unchanged data", {
  N <- 4; Time <- 3; k <- 2
  y <- 1:(N * Time)
  X <- matrix(1:(N * Time * k), ncol = k)

  dm <- demean_panel(y, X, N, Time, model = 0)

  expect_equal(dm$ywith, as.numeric(y))
  expect_equal(dm$xwith, X)
  expect_equal(dm$meanny, rep(0, N))
  expect_equal(dm$meanty, rep(0, Time))
})


test_that("demean_panel model=1 (region FE) subtracts region means correctly", {
  N <- 4; Time <- 3

  # Stacking: [r1t1, r2t1, r3t1, r4t1, r1t2, r2t2, r3t2, r4t2, r1t3, ...]
  # Region 1 obs: positions 1, 5, 9 => values 1, 5, 9 => mean = 5
  # Region 2 obs: positions 2, 6, 10 => values 2, 6, 10 => mean = 6
  # Region 3 obs: positions 3, 7, 11 => values 3, 7, 11 => mean = 7
  # Region 4 obs: positions 4, 8, 12 => values 4, 8, 12 => mean = 8
  y <- as.numeric(1:12)
  X <- matrix(1:24, ncol = 2)

  dm <- demean_panel(y, X, N, Time, model = 1)

  expect_equal(dm$meanny, c(5, 6, 7, 8))

  # ywith for region 1: 1-5=-4, 5-5=0, 9-5=4
  # ywith for region 2: 2-6=-4, 6-6=0, 10-6=4
  expected_ywith <- y - rep(c(5, 6, 7, 8), times = 3)
  expect_equal(dm$ywith, expected_ywith)

  # Demeaned region means should be zero
  ymat <- matrix(dm$ywith, nrow = N, ncol = Time)
  expect_equal(rowMeans(ymat), rep(0, N), tolerance = 1e-14)
})


test_that("demean_panel model=2 (time FE) subtracts time means correctly", {
  N <- 4; Time <- 3

  y <- as.numeric(1:12)
  X <- matrix(1:24, ncol = 2)

  dm <- demean_panel(y, X, N, Time, model = 2)

  # Time period 1: positions 1:4 => values 1,2,3,4 => mean = 2.5
  # Time period 2: positions 5:8 => values 5,6,7,8 => mean = 6.5
  # Time period 3: positions 9:12 => values 9,10,11,12 => mean = 10.5
  expect_equal(dm$meanty, c(2.5, 6.5, 10.5))

  expected_ywith <- y - rep(c(2.5, 6.5, 10.5), each = N)
  expect_equal(dm$ywith, expected_ywith)

  # Demeaned time means should be zero
  ymat <- matrix(dm$ywith, nrow = N, ncol = Time)
  expect_equal(colMeans(ymat), rep(0, Time), tolerance = 1e-14)
})


test_that("demean_panel model=3 (two-way FE) applies correct formula", {
  N <- 4; Time <- 3

  y <- as.numeric(1:12)
  X <- matrix(1:24, ncol = 2)

  dm <- demean_panel(y, X, N, Time, model = 3)

  # Hand computation:
  # Grand mean = mean(1:12) = 6.5
  # Region means: 5, 6, 7, 8
  # Time means: 2.5, 6.5, 10.5
  #
  # ywith[1] = y[1] - meanny[1] - meanty[1] + grand = 1 - 5 - 2.5 + 6.5 = 0
  # ywith[2] = y[2] - meanny[2] - meanty[1] + grand = 2 - 6 - 2.5 + 6.5 = 0
  # ywith[5] = y[5] - meanny[1] - meanty[2] + grand = 5 - 5 - 6.5 + 6.5 = 0
  # In fact, for y = 1:12 with this structure, all ywith should be 0
  # because y_it = i + (t-1)*N is perfectly explained by region + time FE

  grand_mean <- 6.5
  region_means <- c(5, 6, 7, 8)
  time_means <- c(2.5, 6.5, 10.5)

  expected_ywith <- y -
    rep(region_means, times = 3) -
    rep(time_means, each = 4) +
    grand_mean

  expect_equal(dm$ywith, expected_ywith, tolerance = 1e-14)
  expect_equal(dm$meanny, region_means)
  expect_equal(dm$meanty, time_means)

  # For this linear sequence, two-way demeaning should give all zeros
  expect_true(all(abs(dm$ywith) < 1e-12))

  # Both region and time means of demeaned data should be (approximately) zero
  ymat <- matrix(dm$ywith, nrow = N, ncol = Time)
  expect_equal(rowMeans(ymat), rep(0, N), tolerance = 1e-14)
  expect_equal(colMeans(ymat), rep(0, Time), tolerance = 1e-14)
})


test_that("demean_panel returns all six required output components", {
  N <- 4; Time <- 3
  y <- rnorm(N * Time)
  X <- matrix(rnorm(N * Time * 2), ncol = 2)

  for (m in 0:3) {
    dm <- demean_panel(y, X, N, Time, model = m)
    expect_named(dm, c("ywith", "xwith", "meanny", "meannx", "meanty", "meantx"))
    expect_length(dm$ywith, N * Time)
    expect_equal(dim(dm$xwith), c(N * Time, 2))
    expect_length(dm$meanny, N)
    expect_equal(dim(dm$meannx), c(N, 2))
    expect_length(dm$meanty, Time)
    expect_equal(dim(dm$meantx), c(Time, 2))
  }
})


test_that("demean_panel X demeaning is consistent with y demeaning", {
  # If we put y as a column of X, the demeaned versions should match
  N <- 4; Time <- 3
  set.seed(99)
  y <- rnorm(N * Time)
  X <- matrix(c(y, rnorm(N * Time)), ncol = 2)

  for (m in 1:3) {
    dm <- demean_panel(y, X, N, Time, model = m)
    expect_equal(dm$ywith, dm$xwith[, 1], tolerance = 1e-14)
  }
})


test_that("demean_panel rejects invalid inputs", {
  expect_error(demean_panel(1:10, matrix(1:20, ncol = 2), N = 3, Time = 3, model = 1))
  expect_error(demean_panel(1:12, matrix(1:24, ncol = 2), N = 4, Time = 3, model = 5))
})
