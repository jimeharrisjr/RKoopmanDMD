# RKoopmanDMD

Koopman Operator and Dynamic Mode Decomposition for Dynamical Systems Analysis in R

## Overview

RKoopmanDMD implements Dynamic Mode Decomposition (DMD), a data-driven method for analyzing dynamical systems. DMD extracts spatiotemporal coherent structures from time-series data and approximates the Koopman operator - an infinite-dimensional linear operator that captures the evolution of observables for nonlinear systems.

## Installation

```r
# Install from GitHub (development version)
# install.packages("devtools")
devtools::install_github("jimeharrisjr/RKoopmanDMD")
```

## Quick Start

```r
library(RKoopmanDMD)

# Create sample data: damped oscillator
t <- seq(0, 10, by = 0.1)
x <- exp(-0.1 * t) * cos(2 * pi * 0.5 * t)
y <- exp(-0.1 * t) * sin(2 * pi * 0.5 * t)
X <- rbind(x, y)

# Fit DMD model
model <- dmd(X)

# Examine the model
print(model)
summary(model)

# Visualize eigenvalue spectrum
plot(model)

# Predict future states
future <- predict(model, n_ahead = 50)

# Check reconstruction accuracy
errors <- dmd_error(model, X)
cat("Relative error:", errors$relative_error, "\n")
```

## Key Features

- **DMD Analysis**: Extract modes, eigenvalues, and the Koopman matrix from time-series data
- **Lifting Functions**: Transform state space to capture nonlinear dynamics with extended DMD
- **Flexible Prediction**: Forecast from any initial condition using mode-based or matrix methods
- **Stability Analysis**: Determine system stability from eigenvalue spectrum
- **Rich Diagnostics**: Print, summary, and plot methods for model inspection
- **Reconstruction**: Assess model fit by reconstructing training data

## Core Functions

| Function | Description |
|----------|-------------|
| `dmd()` | Perform Dynamic Mode Decomposition |
| `predict()` | Forecast future states |
| `dmd_spectrum()` | Analyze eigenvalue properties |
| `dmd_reconstruct()` | Reconstruct data from modes |
| `dmd_stability()` | Check system stability |
| `dmd_error()` | Compute reconstruction error |
| `dmd_lift()` | Inspect lifting transformation |
| `list_lifting_functions()` | Show available lifting options |

## Advanced Methods

| Function | Description |
|----------|-------------|
| `hankel_dmd()` | Hankel-DMD using Krylov subspace/time-delay embedding |
| `hankel_reconstruct()` | Reconstruct time series from Hankel-DMD model |
| `gla()` | Generalized Laplace Analysis for computing Koopman eigenfunctions |
| `gla_reconstruct()` | Reconstruct signal using GLA model |
| `dmd_dominant_modes()` | Extract most important DMD modes by amplitude/energy |
| `dmd_residual()` | Compute DMD residual for error assessment |
| `dmd_pseudospectrum()` | Pseudospectrum analysis for eigenvalue reliability |
| `dmd_convergence()` | Estimate convergence rate of DMD approximation |

## Understanding the Output

A `dmd` object contains:

- **`A`**: The DMD matrix (Koopman operator approximation). Multiply by current state to get next state.
- **`modes`**: Spatial patterns that evolve with single frequencies
- **`eigenvalues`**: Complex values encoding growth rates and frequencies
  - `|eigenvalue| < 1`: Decaying mode (stable)
  - `|eigenvalue| = 1`: Persistent mode (marginally stable)
  - `|eigenvalue| > 1`: Growing mode (unstable)
- **`amplitudes`**: Initial contribution of each mode

## Lifting Functions for Nonlinear Dynamics

Standard DMD finds a linear operator to approximate dynamics. For systems with nonlinear behavior, **lifting functions** transform the state space into a higher-dimensional feature space where the dynamics become approximately linear. This is the key idea behind Extended DMD (EDMD) and Koopman operator theory.

### When to Use Lifting

- Your system has nonlinear dynamics (e.g., `x' = x²`, oscillators with harmonics)
- Standard DMD gives poor predictions
- You want to capture higher-order interactions between state variables

### Built-in Lifting Functions

| Option | Description | Example |
|--------|-------------|---------|
| `"poly2"`, `"poly3"` | Polynomial powers | `x, x², x³` |
| `"poly_cross2"` | Polynomials with cross-terms | `x, y, x², xy, y²` |
| `"trig"` | Trigonometric | `x, sin(x), cos(x)` |
| `"delay2"`, `"delay3"` | Time-delay embedding | `x(t), x(t-1), x(t-2)` |

### Example: Nonlinear Oscillator

```r
library(RKoopmanDMD)

# System with nonlinear component
t <- seq(0, 10, by = 0.1)
x1 <- cos(t) + 0.3 * cos(t)^2  # Nonlinear term
x2 <- sin(t)
X <- rbind(x1, x2)

# Standard DMD
model_std <- dmd(X)

# DMD with polynomial lifting (captures x² terms)
model_lift <- dmd(X, lifting = "poly2")
print(model_lift)
#> Dynamic Mode Decomposition Model
#> ================================
#>   Original vars:    2
#>   Lifted vars:      4
#>   Lifting:          poly2
#>   ...

# Predictions return original 2D space automatically
pred <- predict(model_lift, n_ahead = 20)
dim(pred)  # 2 x 20

# Inspect the lifting transformation
dmd_lift(model_lift)
#> $original_dim
#> [1] 2
#> $lifted_dim
#> [1] 4
#> $observables
#> [1] 1 2
```

### Custom Lifting Functions

You can define your own lifting function for domain-specific transformations:

```r
# Custom lifting: include specific nonlinear terms
my_lift <- function(X) {
  x <- X[1, ]
  y <- X[2, ]
  rbind(
    X,           # Original states
    x^2,         # x squared
    x * y,       # Interaction term
    sin(x)       # Trigonometric term
  )
}

model <- dmd(X, lifting = my_lift, observables = 1:2)
```

### Parametric Lifting

For fine-grained control, use list specifications:
```r
# Polynomial of degree 5
model <- dmd(X, lifting = list(type = "poly", degree = 5))

# Time-delay embedding with 4 delays
model <- dmd(X, lifting = list(type = "delay", delays = 4))

# RBF lifting with custom centers
centers <- matrix(rnorm(10), nrow = 2, ncol = 5)
model <- dmd(X, lifting = list(type = "rbf", centers = centers, sigma = 0.5))
```

## Example: Satellite Trajectory Prediction

The package includes real satellite position data from the IDAO 2020 Competition. Here's how to predict satellite positions using DMD:

```r
library(RKoopmanDMD)

# Load satellite data
data(sat0)

# Prepare state matrix (position + velocity)
positions <- as.matrix(sat0[, c("x", "y", "z")])
velocities <- as.matrix(sat0[, c("Vx", "Vy", "Vz")])
state_matrix <- t(cbind(positions, velocities))

# Train on first 600 points
X_train <- state_matrix[, 1:600]
model <- dmd(X_train)

# Predict next 100 positions
last_state <- X_train[, ncol(X_train)]
predictions <- predict(model, n_ahead = 100, x0 = last_state)

# 3D visualization with plotly
library(plotly)
plot_ly() %>%
  add_trace(x = X_train[1,], y = X_train[2,], z = X_train[3,],
            type = "scatter3d", mode = "lines",
            line = list(color = "blue"), name = "Training") %>%
  add_trace(x = predictions[1,], y = predictions[2,], z = predictions[3,],
            type = "scatter3d", mode = "lines",
            line = list(color = "red"), name = "Predicted")
```

See `vignette("satellite-prediction")` for a complete tutorial with validation and error analysis.

## Example: Analyzing a Spiral System

```r
# Generate spiral trajectory
set.seed(42)
n <- 100
dt <- 0.1
t <- seq(0, (n-1)*dt, by = dt)
x1 <- exp(-0.05*t) * cos(2*pi*t)
x2 <- exp(-0.05*t) * sin(2*pi*t)
X <- rbind(x1, x2)

# Fit and analyze
model <- dmd(X)

# Check stability
stability <- dmd_stability(model)
cat("System is stable:", stability$is_stable, "\n")
cat("Spectral radius:", stability$spectral_radius, "\n")

# Get detailed spectrum analysis
spectrum <- dmd_spectrum(model, dt = dt)
print(spectrum[, c("mode", "frequency", "growth_rate", "stability")])

# Forecast beyond training data
last_state <- X[, ncol(X)]
forecast <- predict(model, n_ahead = 50, x0 = last_state)

# Plot results
par(mfrow = c(1, 2))
plot(t, X[1,], type = "l", col = "blue",
     xlab = "Time", ylab = "State 1", main = "Training + Forecast")
t_forecast <- max(t) + seq_len(50) * dt
lines(t_forecast, forecast[1,], col = "red", lty = 2)
legend("topright", c("Training", "Forecast"), col = c("blue", "red"), lty = 1:2)

# Phase portrait
plot(X[1,], X[2,], type = "l", col = "blue",
     xlab = "State 1", ylab = "State 2", main = "Phase Portrait")
lines(forecast[1,], forecast[2,], col = "red", lty = 2)
```

## Vignettes

The package includes detailed tutorials on various topics:

| Vignette | Description |
|----------|-------------|
| `vignette("introduction")` | Introduction to RKoopmanDMD |
| `vignette("satellite-prediction")` | Predicting Satellite Positions with DMD |
| `vignette("lifting-functions")` | Lifting Functions for Improved Model Predictions |
| `vignette("hankel-dmd")` | Hankel-DMD: Time-Delay Embedding for Scalar Time Series |
| `vignette("generalized-laplace-analysis")` | Generalized Laplace Analysis: Computing Koopman Eigenfunctions |
| `vignette("diagnostics")` | DMD Diagnostics: Residuals, Pseudospectra, and Convergence |

## References

1. Schmid, P. J. (2010). Dynamic mode decomposition of numerical and experimental data. *Journal of Fluid Mechanics*, 656, 5-28.

2. Kutz, J. N., Brunton, S. L., Brunton, B. W., & Proctor, J. L. (2016). *Dynamic Mode Decomposition: Data-Driven Modeling of Complex Systems*. SIAM.

3. Williams, M. O., Kevrekidis, I. G., & Rowley, C. W. (2015). A data-driven approximation of the Koopman operator. *Journal of Nonlinear Science*, 25(6), 1307-1346.

4. Mezic, I. (2020). On Numerical Approximations of the Koopman Operator. arXiv:2009.05883.

## License

MIT License
