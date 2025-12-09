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

## Understanding the Output

A `dmd` object contains:

- **`A`**: The DMD matrix (Koopman operator approximation). Multiply by current state to get next state.
- **`modes`**: Spatial patterns that evolve with single frequencies
- **`eigenvalues`**: Complex values encoding growth rates and frequencies
  - `|eigenvalue| < 1`: Decaying mode (stable)
  - `|eigenvalue| = 1`: Persistent mode (marginally stable)
  - `|eigenvalue| > 1`: Growing mode (unstable)
- **`amplitudes`**: Initial contribution of each mode

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

## References

1. Schmid, P. J. (2010). Dynamic mode decomposition of numerical and experimental data. *Journal of Fluid Mechanics*, 656, 5-28.

2. Kutz, J. N., Brunton, S. L., Brunton, B. W., & Proctor, J. L. (2016). *Dynamic Mode Decomposition: Data-Driven Modeling of Complex Systems*. SIAM.

3. Williams, M. O., Kevrekidis, I. G., & Rowley, C. W. (2015). A data-driven approximation of the Koopman operator. *Journal of Nonlinear Science*, 25(6), 1307-1346.

## License

MIT License
