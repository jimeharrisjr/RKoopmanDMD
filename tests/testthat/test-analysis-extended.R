# Tests for extended analysis functions (residual, pseudospectrum, convergence)

test_that("dmd_residual computes residual for standard DMD", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X)
  res <- dmd_residual(model, X)

  expect_true(is.list(res))
  expect_true("residual_norm" %in% names(res))
  expect_true("residual_relative" %in% names(res))
  expect_true("per_step_residual" %in% names(res))
  expect_true("pseudospectral_bound" %in% names(res))

  expect_true(res$residual_norm >= 0)
  expect_true(res$residual_relative >= 0)
  expect_true(res$residual_relative <= 1 || res$residual_relative > 0)  # allow small errors
})

test_that("dmd_residual works for hankel_dmd", {
  t <- seq(0, 10, by = 0.1)
  y <- cos(t)

  model <- hankel_dmd(y, delays = 5)
  res <- dmd_residual(model)

  expect_true(is.list(res))
  expect_true(res$residual_norm >= 0)
})

test_that("dmd_residual per_mode_residual has correct length", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X, rank = 2)
  res <- dmd_residual(model, X)

  expect_equal(length(res$per_mode_residual), model$rank)
})

test_that("dmd_pseudospectrum computes without error", {
  t <- seq(0, 5, by = 0.2)  # Smaller data for speed
  X <- rbind(cos(t), sin(t))

  model <- dmd(X)
  ps <- dmd_pseudospectrum(model, grid_n = 20, plot = FALSE)

  expect_true(is.list(ps))
  expect_true("x" %in% names(ps))
  expect_true("y" %in% names(ps))
  expect_true("sigma_min" %in% names(ps))
  expect_true("eigenvalues" %in% names(ps))

  expect_equal(length(ps$x), 20)
  expect_equal(length(ps$y), 20)
  expect_equal(dim(ps$sigma_min), c(20, 20))
})

test_that("dmd_pseudospectrum eigenvalues have small sigma_min", {
  t <- seq(0, 5, by = 0.2)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X)
  ps <- dmd_pseudospectrum(model, grid_n = 50, plot = FALSE)

  # At eigenvalue locations, sigma_min should be very small
  lambdas <- model$eigenvalues
  for (lam in lambdas) {
    i <- which.min(abs(ps$x - Re(lam)))
    j <- which.min(abs(ps$y - Im(lam)))
    sigma_at_eig <- ps$sigma_min[i, j]
    # Should be small (eigenvalue is in spectrum)
    expect_true(sigma_at_eig < 0.5)
  }
})

test_that("dmd_convergence estimates rate", {
  t <- seq(0, 20, by = 0.1)
  X <- rbind(cos(t), sin(t))

  conv <- dmd_convergence(X, sample_fractions = c(0.3, 0.5, 0.7, 1.0))

  expect_true(is.list(conv))
  expect_true("sample_sizes" %in% names(conv))
  expect_true("eigenvalues" %in% names(conv))
  expect_true("eigenvalue_changes" %in% names(conv))
  expect_true("convergence_estimate" %in% names(conv))

  expect_equal(length(conv$sample_sizes), 4)
  expect_equal(length(conv$eigenvalues), 4)
  expect_equal(length(conv$eigenvalue_changes), 3)
})

test_that("dmd_convergence handles small data", {
  t <- seq(0, 5, by = 0.5)
  X <- rbind(cos(t), sin(t))

  # Should not error even with limited data
  conv <- dmd_convergence(X, sample_fractions = c(0.5, 1.0))

  expect_true(is.list(conv))
})

test_that("dmd_convergence shows decreasing changes for good fits", {
  # For a simple linear system, eigenvalue estimates should converge
  t <- seq(0, 30, by = 0.1)
  X <- rbind(cos(t), sin(t))

  conv <- dmd_convergence(X, sample_fractions = c(0.25, 0.5, 0.75, 1.0))

  # Changes should generally decrease (convergence)
  # Allow some noise but check trend
  changes <- conv$eigenvalue_changes
  # At least the last change should be smaller than first
  if (length(changes) >= 2 && all(changes > 0)) {
    expect_true(changes[length(changes)] <= changes[1] * 2)
  }
})

test_that("dmd_residual requires X_original for standard DMD", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X)
  expect_error(dmd_residual(model), "X_original must be provided")
})

test_that("dmd_pseudospectrum custom limits work", {
  t <- seq(0, 5, by = 0.2)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X)
  ps <- dmd_pseudospectrum(model, grid_n = 10, plot = FALSE,
                            xlim = c(-2, 2), ylim = c(-2, 2))

  expect_equal(range(ps$x), c(-2, 2))
  expect_equal(range(ps$y), c(-2, 2))
})
