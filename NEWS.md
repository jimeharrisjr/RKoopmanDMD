# RKoopmanDMD 0.1.0

## Initial Release

This is the first release of RKoopmanDMD, providing tools for analyzing dynamical systems using Dynamic Mode Decomposition (DMD) and the Koopman operator framework.

### Core Features

* `dmd()`: Main function to perform Dynamic Mode Decomposition on time-series data
  - Extracts DMD modes, eigenvalues, and amplitudes
  - Computes the Koopman operator approximation (A matrix)
  - Supports automatic or user-specified rank truncation
  - Optional data centering

* `predict.dmd()`: Forecast future states using a fitted DMD model
  - Two prediction methods: mode-based and matrix multiplication
  - Supports prediction from arbitrary initial conditions
  - Returns properly labeled prediction matrix

### Analysis Functions

* `dmd_spectrum()`: Detailed eigenvalue spectrum analysis including frequencies, growth rates, and stability classification
* `dmd_reconstruct()`: Reconstruct training data from DMD modes (supports mode subsetting)
* `dmd_error()`: Compute reconstruction error metrics (RMSE, MAE, relative error)
* `dmd_stability()`: Check system stability based on eigenvalue magnitudes
* `dmd_dominant_modes()`: Extract most important modes by various criteria
* `dmd_forecast()`: Experimental forecasting with confidence bounds

### S3 Methods

* `print.dmd()`: Concise model summary
* `summary.dmd()`: Detailed model statistics and eigenvalue table
* `plot.dmd()`: Diagnostic plots (eigenvalue spectrum, magnitudes, amplitudes, singular values)
* `coef.dmd()`: Extract eigenvalues, amplitudes, or modes

### Documentation

* Comprehensive function documentation with examples
* Introduction vignette demonstrating typical workflows
* README with quick-start guide

### Deprecated Functions

The following legacy functions are included for backward compatibility but are deprecated:
* `getAMatrix()`: Use `dmd()` instead
* `predict.kdmd()`: Use `predict.dmd()` instead
* `koopman_dmd()`: Use `dmd()` instead
* `predict_dmd()`: Use `predict.dmd()` instead
