# Tests for lifting functions

test_that("polynomial lifting produces correct dimensions", {
  X <- matrix(1:20, nrow = 2, ncol = 10)

  # Degree 2: original + squared = 2 + 2 = 4 rows
  lifted <- RKoopmanDMD:::lift_polynomial(X, degree = 2)
  expect_equal(nrow(lifted), 4)
  expect_equal(ncol(lifted), 10)

  # Degree 3: 2 + 2 + 2 = 6 rows
  lifted3 <- RKoopmanDMD:::lift_polynomial(X, degree = 3)
  expect_equal(nrow(lifted3), 6)

  # First rows should be original data (ignore row names)
  expect_equal(unname(lifted[1:2, ]), unname(X))

  # Second set should be X^2
  expect_equal(unname(lifted[3:4, ]), unname(X^2))
})


test_that("polynomial cross-term lifting works", {
  X <- matrix(1:6, nrow = 2, ncol = 3)

  # Degree 2 for 2 vars: x1, x2, x1^2, x1*x2, x2^2 = 5 terms
  lifted <- RKoopmanDMD:::lift_poly_cross(X, degree = 2)
  expect_equal(nrow(lifted), 5)
  expect_equal(ncol(lifted), 3)

  # Check row names are set
  expect_true(!is.null(rownames(lifted)))
})


test_that("trigonometric lifting works", {
  X <- matrix(seq(0, 2*pi, length.out = 10), nrow = 1)

  # 1 harmonic: x, sin(x), cos(x) = 3 rows
  lifted <- RKoopmanDMD:::lift_trigonometric(X, harmonics = 1)
  expect_equal(nrow(lifted), 3)

  # 2 harmonics: x, sin(x), cos(x), sin(2x), cos(2x) = 5 rows
  lifted2 <- RKoopmanDMD:::lift_trigonometric(X, harmonics = 2)
  expect_equal(nrow(lifted2), 5)

  # Check sin/cos values
  expect_equal(lifted[2, ], sin(X[1, ]))
  expect_equal(lifted[3, ], cos(X[1, ]))
})


test_that("time-delay embedding works", {
  X <- matrix(1:10, nrow = 1, ncol = 10)

  # 2 delays: [x(t), x(t-1), x(t-2)] with 8 valid time points
  lifted <- RKoopmanDMD:::lift_delay(X, delays = 2)
  expect_equal(nrow(lifted), 3)  # 1 * (2 + 1) = 3

  expect_equal(ncol(lifted), 8)  # 10 - 2 = 8

  # Check delay structure (use unname to ignore row names)
  expect_equal(unname(lifted[1, 1]), unname(X[1, 3]))  # First valid time point
  expect_equal(unname(lifted[2, 1]), unname(X[1, 2]))  # t-1
  expect_equal(unname(lifted[3, 1]), unname(X[1, 1]))  # t-2
})


test_that("RBF lifting works", {
  X <- matrix(rnorm(20), nrow = 2, ncol = 10)
  centers <- matrix(rnorm(6), nrow = 2, ncol = 3)

  lifted <- RKoopmanDMD:::lift_rbf(X, centers, sigma = 1)
  expect_equal(nrow(lifted), 5)  # 2 original + 3 RBF
  expect_equal(ncol(lifted), 10)

  # RBF values should be between 0 and 1
  expect_true(all(lifted[3:5, ] >= 0 & lifted[3:5, ] <= 1))
})


test_that("make_lifting_fn handles character shortcuts", {
  X <- matrix(1:20, nrow = 2, ncol = 10)

  # Test poly2
  fn_poly2 <- RKoopmanDMD:::make_lifting_fn("poly2", X)
  expect_true(is.function(fn_poly2))
  expect_equal(nrow(fn_poly2(X)), 4)

  # Test trig
  fn_trig <- RKoopmanDMD:::make_lifting_fn("trig", X)
  expect_true(is.function(fn_trig))
  expect_equal(nrow(fn_trig(X)), 6)  # 2 + 2*2 = 6
})


test_that("make_lifting_fn handles list specification", {
  X <- matrix(1:20, nrow = 2, ncol = 10)

  # Polynomial with custom degree
  fn_poly4 <- RKoopmanDMD:::make_lifting_fn(list(type = "poly", degree = 4), X)
  expect_equal(nrow(fn_poly4(X)), 8)  # 2 * 4 = 8

  # Delay with custom delays
  fn_delay <- RKoopmanDMD:::make_lifting_fn(list(type = "delay", delays = 3), X)
  result <- fn_delay(X)
  expect_equal(nrow(result), 8)  # 2 * (3 + 1) = 8
  expect_equal(ncol(result), 7)  # 10 - 3 = 7
})


test_that("make_lifting_fn handles custom functions", {
  X <- matrix(1:20, nrow = 2, ncol = 10)

  custom_fn <- function(X) rbind(X, X^2, X^3)
  fn <- RKoopmanDMD:::make_lifting_fn(custom_fn, X)
  expect_true(is.function(fn))
  expect_equal(nrow(fn(X)), 6)
})


test_that("dmd works with lifting parameter", {
  # Create simple oscillatory data
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  # Standard DMD
  model_std <- dmd(X)
  expect_null(model_std$lifting)
  expect_equal(model_std$n_vars_original, 2)
  expect_equal(model_std$n_vars_lifted, 2)

  # DMD with polynomial lifting
  model_poly <- dmd(X, lifting = "poly2")
  expect_equal(model_poly$lifting, "poly2")
  expect_equal(model_poly$n_vars_original, 2)
  expect_equal(model_poly$n_vars_lifted, 4)
  expect_equal(model_poly$observables, 1:2)
})


test_that("predict.dmd works with lifted models", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X, lifting = "poly2")

  # Predict returns original dimension by default
  pred <- predict(model, n_ahead = 5)
  expect_equal(nrow(pred), 2)
  expect_equal(ncol(pred), 5)

  # Can return lifted predictions
  pred_lifted <- predict(model, n_ahead = 5, return_lifted = TRUE)
  expect_equal(nrow(pred_lifted), 4)

  # Custom x0 in original space
  x0 <- c(1, 0)
  pred_x0 <- predict(model, n_ahead = 5, x0 = x0)
  expect_equal(nrow(pred_x0), 2)
})


test_that("dmd_lift function works", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X, lifting = "poly2")

  # Get lifting info
  info <- dmd_lift(model)
  expect_equal(info$original_dim, 2)
  expect_equal(info$lifted_dim, 4)
  expect_equal(info$observables, 1:2)

  # Apply lifting to a state
  x0 <- c(1, 0)
  lifted_x0 <- dmd_lift(model, x0)
  expect_equal(length(lifted_x0), 4)
  expect_equal(lifted_x0[1:2], x0)
  expect_equal(lifted_x0[3:4], x0^2)
})


test_that("list_lifting_functions returns expected info", {
  info <- list_lifting_functions()
  expect_true(is.data.frame(info))
  expect_true("name" %in% names(info))
  expect_true("description" %in% names(info))
  expect_true("poly2" %in% info$name)
  expect_true("trig" %in% info$name)
})


test_that("print and summary work with lifted models", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  model <- dmd(X, lifting = "poly2")

  # Print should not error and should mention lifting
  output <- capture.output(print(model))
  expect_true(any(grepl("Lifted", output)))
  expect_true(any(grepl("poly2", output)))

  # Summary should include lifting info
  summ <- summary(model)
  expect_equal(summ$lifting, "poly2")
  expect_equal(summ$n_vars_original, 2)
  expect_equal(summ$n_vars_lifted, 4)
})


test_that("delay embedding reduces time points appropriately", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))
  n_time <- ncol(X)

  # With 3 delays, we lose 3 time points
  model <- dmd(X, lifting = "delay3")
  expect_equal(model$data_dim[2], n_time - 3)
})


test_that("lifting validation catches errors", {
  X <- matrix(1:20, nrow = 2, ncol = 10)

  # Invalid string

  expect_error(RKoopmanDMD:::make_lifting_fn("invalid_type", X))

  # List without type
  expect_error(RKoopmanDMD:::make_lifting_fn(list(degree = 2), X))

  # RBF without centers
  expect_error(RKoopmanDMD:::make_lifting_fn(list(type = "rbf"), X))
})


test_that("lifting works with custom observables", {
  t <- seq(0, 10, by = 0.1)
  X <- rbind(cos(t), sin(t))

  # Custom lifting that reorders
  custom_lift <- function(X) {
    rbind(X^2, X)  # Squared first, then original
  }

  model <- dmd(X, lifting = custom_lift, observables = 3:4)
  expect_equal(model$observables, 3:4)

  # Predictions should return the correct observables
  pred <- predict(model, n_ahead = 5)
  expect_equal(nrow(pred), 2)
})
