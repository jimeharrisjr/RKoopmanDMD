# Tests for dmd() function

test_that("dmd() works with basic input", {
  # Simple oscillator
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X)

  expect_s3_class(model, "dmd")
  expect_true(!is.null(model$A))
  expect_true(!is.null(model$modes))
  expect_true(!is.null(model$eigenvalues))
  expect_true(!is.null(model$amplitudes))
})

test_that("dmd() returns correct dimensions", {
  n_vars <- 3
  n_time <- 50
  X <- matrix(rnorm(n_vars * n_time), nrow = n_vars, ncol = n_time)

  model <- dmd(X)

  expect_equal(nrow(model$A), n_vars)
  expect_equal(ncol(model$A), n_vars)
  expect_equal(nrow(model$modes), n_vars)
  expect_equal(length(model$eigenvalues), model$rank)
  expect_equal(length(model$amplitudes), model$rank)
  expect_equal(model$data_dim, c(n_vars, n_time))
})

test_that("dmd() respects rank parameter", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t), cos(2*t))

  model_r2 <- dmd(X, rank = 2)
  model_auto <- dmd(X)

  expect_equal(model_r2$rank, 2)
  expect_equal(length(model_r2$eigenvalues), 2)
  expect_gte(model_auto$rank, 1)
})

test_that("dmd() handles centering", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t) + 5, sin(t) + 3)  # Non-zero mean

  model_centered <- dmd(X, center = TRUE)
  model_uncentered <- dmd(X, center = FALSE)

  expect_true(model_centered$center)
  expect_false(model_uncentered$center)
  expect_false(is.null(model_centered$X_mean))
  expect_true(is.null(model_uncentered$X_mean))
})

test_that("dmd() validates input correctly", {
  # Not a matrix

  expect_error(dmd(1:10), "must be a matrix")

  # Too few columns (need at least 3)
  expect_error(dmd(matrix(1:4, nrow = 2)), "at least 3 columns")

  # Non-numeric
  expect_error(dmd(matrix(letters[1:10], nrow = 2)), "numeric")

  # Contains NA
  X <- matrix(1:20, nrow = 2)
  X[1, 1] <- NA
  expect_error(dmd(X), "non-finite")

  # Contains Inf
  X <- matrix(1:20, nrow = 2)
  X[1, 1] <- Inf
  expect_error(dmd(X), "non-finite")
})

test_that("dmd() works with data.frame input", {
  df <- data.frame(
    t1 = c(1, 2, 3),
    t2 = c(2, 3, 4),
    t3 = c(3, 4, 5),
    t4 = c(4, 5, 6)
  )

  expect_no_error(model <- dmd(as.matrix(t(df))))
  expect_s3_class(model, "dmd")
})

test_that("dmd() captures oscillatory dynamics", {
  # Create pure oscillator (eigenvalues should be on unit circle)
  t <- seq(0, 20, by = 0.1)
  freq <- 0.5
  X <- rbind(cos(2 * pi * freq * t), sin(2 * pi * freq * t))

  model <- dmd(X)

  # Eigenvalue magnitudes should be close to 1
  mags <- Mod(model$eigenvalues)
  expect_true(all(abs(mags - 1) < 0.01))
})

test_that("dmd() captures decaying dynamics", {
  # Damped oscillator
  t <- seq(0, 10, by = 0.1)
  decay <- 0.1
  X <- rbind(
    exp(-decay * t) * cos(t),
    exp(-decay * t) * sin(t)
  )

  model <- dmd(X)

  # Eigenvalue magnitudes should be less than 1
  mags <- Mod(model$eigenvalues)
  expect_true(all(mags < 1))
})

test_that("dmd() stores call correctly", {
  X <- matrix(rnorm(30), nrow = 3)
  model <- dmd(X, rank = 2)

  expect_true(!is.null(model$call))
  expect_equal(model$call[[1]], as.name("dmd"))
})
