# Tests for Generalized Laplace Analysis (GLA)

test_that("gla works with scalar time series", {
  t <- seq(0, 20, by = 0.1)
  y <- cos(2 * pi * 0.1 * t)

  result <- gla(y, n_eigenvalues = 4)

  expect_s3_class(result, "gla")
  expect_equal(result$n_obs, 1)
  # n_eigenvalues is a max; actual number depends on preliminary DMD
  expect_true(length(result$eigenvalues) >= 1)
  expect_true(length(result$eigenvalues) <= 4)
  expect_equal(nrow(result$modes), 1)
})

test_that("gla works with multivariate time series", {
  t <- seq(0, 20, by = 0.1)
  X <- rbind(cos(t), sin(t))

  result <- gla(X, n_eigenvalues = 3)

  expect_s3_class(result, "gla")
  expect_equal(result$n_obs, 2)
  expect_equal(nrow(result$modes), 2)
  expect_true(length(result$eigenvalues) >= 1)
})

test_that("gla works with provided eigenvalues", {
  dt <- 0.1
  t <- seq(0, 20, by = dt)
  freq <- 0.1
  y <- cos(2 * pi * freq * t)

  # Provide known eigenvalues
  omega <- 2 * pi * freq * dt
  known_eigs <- c(exp(1i * omega), exp(-1i * omega))

  result <- gla(y, eigenvalues = known_eigs)

  expect_equal(length(result$eigenvalues), 2)
})

test_that("gla eigenfunctions have eigenvalue relationship", {
  t <- seq(0, 30, by = 0.1)
  y <- cos(2 * pi * 0.05 * t)  # Low frequency for better convergence

  result <- gla(y, n_eigenvalues = 2)

  # eigenfunction_errors should exist and be finite

  expect_true(length(result$eigenfunction_errors) >= 1)
  expect_true(all(is.finite(result$eigenfunction_errors)))

  # For converged modes, eigenfunction error should be small
  for (k in seq_along(result$eigenvalues)) {
    if (result$convergence[k]) {
      expect_true(result$eigenfunction_errors[k] < 0.5)
    }
  }
})

test_that("gla print method works", {
  t <- seq(0, 20, by = 0.1)
  y <- cos(t)

  result <- gla(y, n_eigenvalues = 3)

  expect_output(print(result), "Generalized Laplace Analysis")
  expect_output(print(result), "Eigenvalues")
})

test_that("gla predict works", {
  t <- seq(0, 20, by = 0.1)
  y <- cos(t)

  result <- gla(y, n_eigenvalues = 4)
  pred <- predict(result, n_ahead = 10)

  expect_equal(nrow(pred), 1)
  expect_equal(ncol(pred), 10)
  expect_true(all(is.finite(pred)))
})

test_that("gla_reconstruct works", {
  t <- seq(0, 20, by = 0.1)
  y <- cos(t)

  result <- gla(y, n_eigenvalues = 4)
  recon <- gla_reconstruct(result)

  expect_equal(nrow(recon), 1)
  expect_equal(ncol(recon), result$n_time)
  expect_true(all(is.finite(recon)))
})

test_that("gla_reconstruct with subset of modes works", {
  t <- seq(0, 20, by = 0.1)
  y <- cos(t)

  result <- gla(y, n_eigenvalues = 4)
  recon <- gla_reconstruct(result, modes_to_use = 1:2)

  expect_equal(nrow(recon), 1)
  expect_equal(ncol(recon), result$n_time)
})

test_that("gla validates inputs", {
  expect_error(gla("not numeric"), "numeric")
  expect_error(gla(c(1, NA, 3)), "non-finite")
})

test_that("gla handles max_iter parameter", {
  t <- seq(0, 20, by = 0.1)
  y <- cos(t)

  result <- gla(y, n_eigenvalues = 2, max_iter = 50)

  expect_equal(result$n_iter, 50)
})

test_that("gla residuals are computed", {
  t <- seq(0, 20, by = 0.1)
  y <- cos(t)

  result <- gla(y, n_eigenvalues = 3)

  # residuals should match number of eigenvalues found
  expect_equal(length(result$residuals), length(result$eigenvalues))
  expect_true(all(is.finite(result$residuals)))
})

test_that("gla convergence is tracked", {
  t <- seq(0, 20, by = 0.1)
  y <- cos(t)

  result <- gla(y, n_eigenvalues = 3)

  # convergence should match number of eigenvalues found
  expect_equal(length(result$convergence), length(result$eigenvalues))
  expect_true(all(is.logical(result$convergence)))
})
