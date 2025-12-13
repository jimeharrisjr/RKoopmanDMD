# Tests for Hankel-DMD (Krylov subspace method)

test_that("hankel_dmd works with scalar time series", {
  # Simple oscillator
  t <- seq(0, 10, by = 0.1)
  y <- cos(2 * pi * 0.5 * t)

  model <- hankel_dmd(y, delays = 5)

  expect_s3_class(model, "hankel_dmd")
  expect_s3_class(model, "dmd")
  expect_equal(model$n_obs, 1)
  expect_equal(model$delays, 5)
  expect_equal(nrow(model$hankel), 6)  # delays + 1
})

test_that("hankel_dmd works with multivariate time series", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  model <- hankel_dmd(X, delays = 3)

  expect_s3_class(model, "hankel_dmd")
  expect_equal(model$n_obs, 2)
  expect_equal(model$delays, 3)
  expect_equal(nrow(model$hankel), 8)  # (delays + 1) * n_obs
})

test_that("hankel_dmd detects oscillation frequencies", {
  dt <- 0.1
  t <- seq(0, 20, by = dt)
  freq <- 0.5  # Hz
  y <- cos(2 * pi * freq * t)

  model <- hankel_dmd(y, delays = 10, dt = dt)
  spectrum <- dmd_spectrum(model, dt = dt)

  # Should find eigenvalues near exp(±i * 2*pi*freq*dt)
  expected_phase <- 2 * pi * freq * dt
  phases <- abs(spectrum$phase)

  expect_true(any(abs(phases - expected_phase) < 0.1))
})

test_that("hankel_dmd predictions work", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t)

  model <- hankel_dmd(y, delays = 5)
  pred <- predict(model, n_ahead = 10)

  expect_equal(nrow(pred), 1)  # 1 observable

  expect_equal(ncol(pred), 10)
  expect_true(all(is.finite(pred)))
})

test_that("hankel_dmd predictions with return_full works", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t)

  model <- hankel_dmd(y, delays = 5)
  pred_full <- predict(model, n_ahead = 10, return_full = TRUE)

  expect_equal(nrow(pred_full), 6)  # delays + 1
  expect_equal(ncol(pred_full), 10)
})

test_that("build_hankel_matrix creates correct structure", {
  y <- matrix(1:10, nrow = 1)
  H <- RKoopmanDMD:::build_hankel_matrix(y, delays = 2)

  expect_equal(nrow(H), 3)  # delays + 1
  expect_equal(ncol(H), 8)  # n_time - delays

  # Check Hankel structure - each column is a delay-embedded snapshot
  # Row 1 = y(t), Row 2 = y(t+1), Row 3 = y(t+2)
  # Column 1: starts at t=1, so y[1], y[2], y[3]
  # Using unname to ignore row names
  expect_equal(unname(H[1, 1]), 1)   # y[1]
  expect_equal(unname(H[2, 1]), 2)   # y[2]
  expect_equal(unname(H[3, 1]), 3)   # y[3]
  expect_equal(unname(H[1, 2]), 2)   # y[2] - next column shifts by 1
})

test_that("hankel_reconstruct works", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t)

  model <- hankel_dmd(y, delays = 5)
  recon <- hankel_reconstruct(model)

  expect_equal(nrow(recon), 1)
  expect_equal(ncol(recon), ncol(model$hankel))
  expect_true(all(is.finite(recon)))
})

test_that("hankel_dmd handles automatic delay selection", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t)

  model <- hankel_dmd(y)  # delays = NULL

  expect_true(model$delays >= 2)
  expect_true(model$delays < length(y) / 2)
})

test_that("hankel_dmd residual is computed", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t) + 0.1 * rnorm(length(t))

  model <- hankel_dmd(y, delays = 5)

  expect_true(!is.null(model$residual))
  expect_true(is.numeric(model$residual))
  expect_true(model$residual >= 0)
})

test_that("hankel_dmd print method works", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t)

  model <- hankel_dmd(y, delays = 5)

  expect_output(print(model), "Hankel-DMD")
  expect_output(print(model), "Delays used")
})

test_that("hankel_dmd works with dmd_spectrum", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t)

  model <- hankel_dmd(y, delays = 5)
  spectrum <- dmd_spectrum(model)

  expect_s3_class(spectrum, "data.frame")
  expect_true("magnitude" %in% names(spectrum))
  expect_true("frequency" %in% names(spectrum))
})

test_that("hankel_dmd works with dmd_stability", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t)

  model <- hankel_dmd(y, delays = 5)
  stability <- dmd_stability(model)

  expect_true(is.list(stability))
  expect_true("is_stable" %in% names(stability))
})

test_that("hankel_dmd validates inputs", {
  expect_error(hankel_dmd("not a vector"), "numeric")
  expect_error(hankel_dmd(c(1, NA, 3)), "non-finite")
  expect_error(hankel_dmd(1:5, delays = 10), "less than")
})
