# Tests for predict.dmd() and related prediction functions

test_that("predict.dmd() returns correct dimensions", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  n_ahead <- 20
  preds <- predict(model, n_ahead = n_ahead)

  expect_equal(nrow(preds), nrow(X))
  expect_equal(ncol(preds), n_ahead)
})

test_that("predict.dmd() works with both methods", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  pred_modes <- predict(model, n_ahead = 10, method = "modes")
  pred_matrix <- predict(model, n_ahead = 10, method = "matrix")

  expect_equal(dim(pred_modes), dim(pred_matrix))
  # Results should be similar (not identical due to numerical differences)
  expect_lt(max(abs(pred_modes - pred_matrix)), 0.1)
})

test_that("predict.dmd() respects custom initial condition", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  x0 <- c(0.5, 0.5)
  preds <- predict(model, n_ahead = 10, x0 = x0)

  expect_equal(nrow(preds), 2)
  expect_equal(ncol(preds), 10)
})

test_that("predict.dmd() validates x0 dimension", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))  # 2 state variables
  model <- dmd(X)

  # Wrong dimension for x0
  expect_error(predict(model, n_ahead = 10, x0 = c(1, 2, 3)),
               "length 2")
})

test_that("predict.dmd() validates n_ahead", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  expect_error(predict(model, n_ahead = 0), "positive")
  expect_error(predict(model, n_ahead = -5), "positive")
})

test_that("predict.dmd() continues oscillatory trajectory", {
  # Pure oscillator should continue oscillating
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  preds <- predict(model, n_ahead = 50)

  # Check that predictions maintain similar amplitude
  orig_amp <- sqrt(mean(X[1,]^2 + X[2,]^2))
  pred_amp <- sqrt(mean(preds[1,]^2 + preds[2,]^2))

  expect_lt(abs(orig_amp - pred_amp) / orig_amp, 0.1)
})

test_that("predict.dmd() decays for damped system", {
  # Damped oscillator predictions should decay
  t <- seq(0, 10, by = 0.1)
  X <- rbind(
    exp(-0.2 * t) * cos(t),
    exp(-0.2 * t) * sin(t)
  )
  model <- dmd(X)

  preds <- predict(model, n_ahead = 50)

  # Last predictions should be smaller than first
  early_amp <- sqrt(preds[1, 1]^2 + preds[2, 1]^2)
  late_amp <- sqrt(preds[1, 50]^2 + preds[2, 50]^2)

  expect_lt(late_amp, early_amp)
})

test_that("predict.dmd() adds column names", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  preds <- predict(model, n_ahead = 5)

  expect_equal(colnames(preds), c("t+1", "t+2", "t+3", "t+4", "t+5"))
})

test_that("dmd_reconstruct() matches original data", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  recon <- dmd_reconstruct(model)

  # Reconstruction should be close to original
  rel_error <- sqrt(sum((X - recon)^2)) / sqrt(sum(X^2))
  expect_lt(rel_error, 0.01)
})

test_that("dmd_reconstruct() works with mode subset", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  # Use only first mode
  recon_1 <- dmd_reconstruct(model, modes = 1)

  expect_equal(nrow(recon_1), nrow(X))
  expect_equal(ncol(recon_1), ncol(X))
})

test_that("dmd_error() computes correct metrics", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  errors <- dmd_error(model, X)

  expect_true("rmse" %in% names(errors))
  expect_true("mae" %in% names(errors))
  expect_true("relative_error" %in% names(errors))
  expect_true("per_variable" %in% names(errors))

  # Error should be small for clean oscillator
  expect_lt(errors$relative_error, 0.01)
})

test_that("dmd_forecast() returns confidence bounds", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  model <- dmd(X)

  forecast <- dmd_forecast(model, n_ahead = 10)

  expect_true("fit" %in% names(forecast))
  expect_true("lower" %in% names(forecast))
  expect_true("upper" %in% names(forecast))
  expect_true("level" %in% names(forecast))

  # Upper should be greater than lower
  expect_true(all(forecast$upper >= forecast$lower))
})
