# Koopman Dynamic Mode Decomposition (DMD) in R
# This implementation learns dynamics from data and makes predictions

# Function to perform DMD
koopman_dmd <- function(X, r = NULL) {
  # X: data matrix where columns are snapshots over time
  # r: rank for truncation (if NULL, uses full rank)

  # Split data into X (current) and Y (next state)
  X1 <- X[, 1:(ncol(X) - 1)]
  X2 <- X[, 2:ncol(X)]

  # Perform SVD on X1
  svd_result <- svd(X1)
  U <- svd_result$u
  S <- svd_result$d
  V <- svd_result$v

  # Determine rank
  if (is.null(r)) {
    r <- length(S)
  } else {
    r <- min(r, length(S))
  }

  # Truncate SVD
  U_r <- U[, 1:r]
  S_r <- S[1:r]
  V_r <- V[, 1:r]

  # Compute DMD matrix (reduced coordinates)
  A_tilde <- t(U_r) %*% X2 %*% V_r %*% diag(1/S_r, nrow = r)

  # Eigendecomposition of A_tilde
  eigen_result <- eigen(A_tilde)
  eigenvalues <- eigen_result$values
  W <- eigen_result$vectors

  # Compute DMD modes
  modes <- X2 %*% V_r %*% diag(1/S_r, nrow = r) %*% W

  # Compute initial amplitudes
  b <- solve(modes, X[, 1])

  return(list(
    modes = modes,
    eigenvalues = eigenvalues,
    amplitudes = b,
    rank = r
  ))
}

# Function to predict future states using DMD
predict_dmd <- function(dmd_result, n_steps, x0 = NULL) {
  # dmd_result: output from koopman_dmd()
  # n_steps: number of time steps to predict
  # x0: initial condition - the state from which to start prediction
  #     - If NULL: uses stored amplitudes (reconstructs training trajectory)
  #     - If provided: projects x0 onto DMD modes to get new amplitudes

  modes <- dmd_result$modes
  lambdas <- dmd_result$eigenvalues

  if (is.null(x0)) {
    # Use amplitudes from training (b = Modes^(-1) * X[:,1])
    # This reconstructs the learned trajectory
    b <- dmd_result$amplitudes
  } else {
    # Project new initial condition onto DMD modes
    # b = Modes^(-1) * x0
    # This allows prediction from ANY starting point
    b <- solve(modes, x0)
  }

  # Predict future states
  predictions <- matrix(0, nrow = nrow(modes), ncol = n_steps)

  for (k in 1:n_steps) {
    # x(t) = Σ b_i * λ_i^t * φ_i
    # Each mode evolves independently with its eigenvalue
    time_evolution <- lambdas^(k - 1)
    predictions[, k] <- Re(modes %*% (b * time_evolution))
  }

  return(predictions)
}

# Example: Generate synthetic data from a simple dynamical system
set.seed(42)
n <- 100  # number of time steps
dt <- 0.1

# Create a spiral/oscillatory system
t <- seq(0, (n-1)*dt, by = dt)
x1 <- exp(-0.05*t) * cos(2*pi*t)
x2 <- exp(-0.05*t) * sin(2*pi*t)
X <- rbind(x1, x2)

# Learn DMD model
dmd_model <- koopman_dmd(X, r = 2)

cat("DMD Analysis Results:\n")
cat("Rank used:", dmd_model$rank, "\n")
cat("Eigenvalues (first 5):\n")
print(head(dmd_model$eigenvalues, 5))

# Predict next 50 steps
n_predict <- 50
predictions <- predict_dmd(dmd_model, n_predict, x0 = X[, ncol(X)])

# Visualize results
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

# Plot original data and predictions
plot(t, x1, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "State Variable 1",
     main = "DMD: Original Data and Predictions")
t_pred <- seq(max(t) + dt, max(t) + n_predict*dt, by = dt)
lines(t_pred, predictions[1, ], col = "red", lwd = 2, lty = 2)
legend("topright", c("Original", "Predicted"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# Phase portrait
plot(x1, x2, type = "l", col = "blue", lwd = 2,
     xlab = "State 1", ylab = "State 2",
     main = "Phase Portrait")
lines(predictions[1, ], predictions[2, ], col = "red", lwd = 2, lty = 2)
legend("topright", c("Original", "Predicted"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

# Print reconstruction error
X_test <- X[, 1:min(50, ncol(X))]
X_recon <- predict_dmd(dmd_model, ncol(X_test), x0 = X[, 1])
error <- norm(X_test - X_recon, "F") / norm(X_test, "F")
cat("\nRelative reconstruction error:", round(error, 6), "\n")

# =============================================================================
# PRACTICAL EXAMPLES OF x0 USAGE
# =============================================================================

cat("\n=== EXAMPLE 1: Reconstructing Training Trajectory ===\n")
# When x0 = NULL, uses stored amplitudes from training data
recon_from_training <- predict_dmd(dmd_model, 30, x0 = NULL)
cat("Prediction starts from first training snapshot\n")
cat("Shape:", dim(recon_from_training), "\n")

cat("\n=== EXAMPLE 2: Forecasting from Last Observed State ===\n")
# Predict future from the last observed data point
last_state <- X[, ncol(X)]
forecast <- predict_dmd(dmd_model, 50, x0 = last_state)
cat("Starting from last observed state:", round(last_state, 3), "\n")
cat("Forecasting", ncol(forecast), "steps into the future\n")

cat("\n=== EXAMPLE 3: Prediction from Arbitrary State ===\n")
# What if system starts from a different condition?
new_initial <- c(0.5, 0.3)  # Custom initial state
custom_pred <- predict_dmd(dmd_model, 40, x0 = new_initial)
cat("Starting from custom state:", new_initial, "\n")
cat("DMD projects this onto learned modes and evolves forward\n")

cat("\n=== EXAMPLE 4: Interpolation - Predict from Middle of Dataset ===\n")
# Start from middle of training data and see if it matches
mid_point <- X[, 50]
interp_pred <- predict_dmd(dmd_model, 20, x0 = mid_point)
actual_continuation <- X[, 51:70]
interp_error <- norm(interp_pred - actual_continuation, "F") / norm(actual_continuation, "F")
cat("Interpolation error:", round(interp_error, 6), "\n")
cat("This tests how well DMD captures the dynamics\n")

# Visualization of x0 usage
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Plot 1: Reconstruction (x0 = NULL)
plot(1:ncol(X), X[1,], type = "l", col = "blue", lwd = 2,
     main = "x0 = NULL: Reconstruct Training",
     xlab = "Time Step", ylab = "State 1")
lines(1:ncol(recon_from_training), recon_from_training[1,],
      col = "red", lwd = 2, lty = 2)
legend("topright", c("Original", "DMD Reconstruction"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, cex = 0.8)

# Plot 2: Forecasting (x0 = last state)
plot((ncol(X)-10):ncol(X), X[1, (ncol(X)-10):ncol(X)],
     type = "l", col = "blue", lwd = 2,
     xlim = c(ncol(X)-10, ncol(X) + ncol(forecast)),
     ylim = range(c(X[1,], forecast[1,])),
     main = "x0 = Last State: Forecast",
     xlab = "Time Step", ylab = "State 1")
lines((ncol(X)+1):(ncol(X)+ncol(forecast)), forecast[1,],
      col = "red", lwd = 2, lty = 2)
abline(v = ncol(X), lty = 3, col = "gray")
legend("topright", c("Historical", "Forecast"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, cex = 0.8)

# Plot 3: Custom initial condition
plot(1:ncol(custom_pred), custom_pred[1,], type = "l", col = "red", lwd = 2,
     main = "x0 = Custom State [0.5, 0.3]",
     xlab = "Time Step", ylab = "State 1")
points(1, new_initial[1], pch = 19, col = "darkred", cex = 2)
text(5, new_initial[1], "Starting point", pos = 3)

# Plot 4: Interpolation test
plot(50:70, X[1, 50:70], type = "l", col = "blue", lwd = 2,
     main = "x0 = Mid-point: Interpolation Test",
     xlab = "Time Step", ylab = "State 1")
lines(51:70, interp_pred[1,], col = "red", lwd = 2, lty = 2)
points(50, mid_point[1], pch = 19, col = "darkgreen", cex = 2)
legend("topright", c("Actual", "DMD Predicted"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, cex = 0.8)
