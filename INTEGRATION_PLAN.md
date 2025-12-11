# RKoopmanDMD Integration Plan

## Overview

This plan outlines the steps to integrate the Koopman Dynamic Mode Decomposition (DMD) functions from `claudecode_sugested.R` into the existing package structure, creating a CRAN-worthy R package for analyzing dynamical systems.

---

## Current State Analysis

### Existing Package Structure
```
RKoopmanDMD/
├── DESCRIPTION          # Placeholder metadata
├── NAMESPACE            # Generic export pattern
├── R/
│   ├── main.R           # Basic kdmd class + getAMatrix + predict.kdmd
│   └── claudecode_sugested.R  # koopman_dmd + predict_dmd + examples
├── man/
│   └── hello.Rd         # Placeholder documentation
└── RKoopmanDMD.Rproj
```

### Code to Integrate

**From `main.R`** (existing):
- `kdmd()` - Constructor for kdmd class
- `getAMatrix()` - Computes Koopman matrix using DMD
- `predict.kdmd()` - S3 predict method for kdmd objects

**From `claudecode_sugested.R`** (to integrate):
- `koopman_dmd()` - Alternative DMD implementation returning modes, eigenvalues, amplitudes
- `predict_dmd()` - Prediction function supporting arbitrary initial conditions

### Issues in Current Code

1. **main.R line 55**: `file_connection` variable is undefined (debugging artifact)
2. **main.R line 17**: Typo in example `m[3,]<-3,12` should be `m[3,]<-3:12`
3. **Inconsistent interfaces**: Two different DMD implementations with different return types
4. **Missing documentation**: Roxygen headers incomplete or missing
5. **No input validation**: `claudecode_sugested.R` lacks defensive checks
6. **Dependency issues**: `pracma::pinv` import needs proper handling

---

## Integration Strategy

### Design Decision: Unified API

Rather than keeping two parallel implementations, create a unified interface that:
1. Uses a single `dmd()` constructor function
2. Returns a rich S3 object containing modes, eigenvalues, amplitudes, AND the A matrix
3. Provides S3 methods for `predict()`, `print()`, `summary()`, `plot()`
4. Supports both use cases: matrix multiplication prediction AND mode-based forecasting

---

## Implementation Plan

### Phase 1: Package Infrastructure

#### 1.1 Update DESCRIPTION
- Set proper Title, Author, Maintainer, Description
- Add License (suggest MIT or GPL-3)
- Declare Imports: `pracma` (for pinv), `graphics`, `stats`
- Suggests: `testthat`, `knitr`, `rmarkdown`
- Set minimum R version (>= 3.5.0)
- Add URL and BugReports fields

#### 1.2 Update NAMESPACE
- Replace `exportPattern` with explicit `export()` statements
- Add `importFrom()` for external dependencies
- Register S3 methods with `S3method()`

#### 1.3 Create Package-Level Documentation
- Create `R/RKoopmanDMD-package.R` with package documentation
- Include `@docType package` roxygen block

---

### Phase 2: Core Functions

#### 2.1 Create `R/dmd.R` - Main DMD Functions

```r
# Contents:
# - dmd() - Main constructor (combining best of both implementations)
# - validate_data() - Internal helper for input validation
```

**Function: `dmd()`**
- Parameters: `X` (data matrix), `rank` (truncation rank, default NULL = auto), `method` (default "exact")
- Returns: S3 object of class "dmd" containing:
  - `A` - The Koopman/DMD matrix
  - `modes` - DMD modes (eigenvectors in physical space)
  - `eigenvalues` - Complex eigenvalues
  - `amplitudes` - Initial mode amplitudes
  - `rank` - Rank used
  - `U`, `V`, `S` - SVD components (for advanced users)
  - `call` - The matched call
  - `data_dim` - Original data dimensions

#### 2.2 Create `R/predict.R` - Prediction Methods

```r
# Contents:
# - predict.dmd() - S3 predict method
```

**Function: `predict.dmd()`**
- Parameters:
  - `object` - dmd object
  - `n_ahead` - Number of steps to forecast
  - `x0` - Initial state (default: last column of training data)
  - `method` - "matrix" (A^n approach) or "modes" (eigenvalue evolution)
  - `...` - For S3 compatibility
- Returns: Matrix of predictions

#### 2.3 Create `R/methods.R` - S3 Methods

```r
# Contents:
# - print.dmd()
# - summary.dmd()
# - plot.dmd()
# - coef.dmd() - Extract eigenvalues/modes
```

#### 2.4 Create `R/analysis.R` - Analysis Functions

```r
# Contents:
# - dmd_spectrum() - Analyze eigenvalue spectrum
# - dmd_stability() - Check system stability (|eigenvalues| relative to 1)
# - dmd_reconstruct() - Reconstruct training data
# - dmd_error() - Compute reconstruction error
```

#### 2.5 Create `R/utils.R` - Utility Functions

```r
# Contents:
# - validate_matrix() - Input validation
# - determine_rank() - Auto rank selection
# - safe_pinv() - Wrapper for pseudo-inverse with fallback
```

---

### Phase 3: Documentation

#### 3.1 Function Documentation
All exported functions get complete roxygen2 documentation:
- `@title` - Brief title
- `@description` - Detailed description
- `@param` - All parameters documented
- `@return` - Return value description
- `@examples` - Runnable examples
- `@references` - Key DMD/Koopman papers
- `@seealso` - Related functions
- `@export` / `@keywords internal`

#### 3.2 Create Vignettes

**Vignette 1: `vignettes/introduction.Rmd`**
- What is Koopman/DMD?
- Basic usage workflow
- Simple examples

**Vignette 2: `vignettes/dynamical-systems.Rmd`**
- Understanding eigenvalues and stability
- Phase portraits
- Real-world dynamical systems examples

**Vignette 3: `vignettes/advanced-usage.Rmd`**
- Rank selection strategies
- Handling noisy data
- Comparison with other methods

#### 3.3 README.md
- Package overview
- Installation instructions
- Quick start example
- Badge placeholders for CRAN, CI, coverage

---

### Phase 4: Testing

#### 4.1 Unit Tests (`tests/testthat/`)

```
tests/
├── testthat.R
└── testthat/
    ├── test-dmd.R           # Core DMD function tests
    ├── test-predict.R       # Prediction tests
    ├── test-validation.R    # Input validation tests
    ├── test-analysis.R      # Analysis function tests
    └── test-edge-cases.R    # Edge cases and error handling
```

**Test categories:**
- Input validation (wrong types, dimensions, NA handling)
- Numerical accuracy (known analytical systems)
- Edge cases (single mode, full rank, near-singular)
- Prediction accuracy (reconstruction error bounds)
- S3 method dispatch

#### 4.2 Example Data
Create `data/` directory with:
- `lorenz` - Lorenz attractor sample data
- `oscillator` - Damped harmonic oscillator
- `vanderpol` - Van der Pol oscillator

---

### Phase 5: CRAN Compliance

#### 5.1 R CMD check Requirements
- No ERRORs, WARNINGs, or NOTEs (ideally)
- Examples run in < 5 seconds each
- Vignettes build without error
- All documentation complete

#### 5.2 Files to Add
- `NEWS.md` - Changelog
- `cran-comments.md` - CRAN submission notes
- `.Rbuildignore` - Update to exclude dev files
- `LICENSE` / `LICENSE.md` - License file
- `inst/CITATION` - Citation information

#### 5.3 Continuous Integration
- GitHub Actions workflow for R CMD check
- Code coverage with codecov
- pkgdown site generation

---

## Detailed Function Specifications

### `dmd()` - Main Constructor

```r
dmd <- function(X, rank = NULL, center = FALSE) {

  # Validate input
  X <- validate_matrix(X)

  # Optionally center data
  if (center) {
    X_mean <- rowMeans(X)
    X <- X - X_mean
  }

  # Split into current/next state pairs
  X1 <- X[, -ncol(X), drop = FALSE]
  X2 <- X[, -1, drop = FALSE]

  # SVD of X1
  svd_result <- svd(X1)

  # Determine rank (auto or specified)
  r <- determine_rank(svd_result$d, rank)

  # Truncate
  U_r <- svd_result$u[, 1:r, drop = FALSE]
  S_r <- svd_result$d[1:r]
  V_r <- svd_result$v[, 1:r, drop = FALSE]

  # Compute reduced DMD matrix
  A_tilde <- t(U_r) %*% X2 %*% V_r %*% diag(1/S_r, nrow = r, ncol = r)

  # Eigendecomposition
  eig <- eigen(A_tilde)

  # DMD modes in physical space
  Phi <- X2 %*% V_r %*% diag(1/S_r, nrow = r, ncol = r) %*% eig$vectors

  # Full A matrix (optional, for matrix-based prediction)
  A <- Phi %*% diag(eig$values) %*% pracma::pinv(Phi)

  # Initial amplitudes
  b <- solve(Phi, X[, 1])

  structure(
    list(
      A = Re(A),
      modes = Phi,
      eigenvalues = eig$values,
      amplitudes = b,
      rank = r,
      svd = list(U = U_r, S = S_r, V = V_r),
      A_tilde = A_tilde,
      X_last = X[, ncol(X)],
      data_dim = dim(X),
      center = center,
      X_mean = if (center) X_mean else NULL,
      call = match.call()
    ),
    class = "dmd"
  )
}
```

### `predict.dmd()` - Prediction Method

```r
predict.dmd <- function(object, n_ahead = 1, x0 = NULL,
                        method = c("modes", "matrix"), ...) {
  method <- match.arg(method)

  if (is.null(x0)) {
    x0 <- object$X_last
  }

  if (method == "matrix") {
    # Direct matrix multiplication
    predictions <- matrix(0, nrow = length(x0), ncol = n_ahead)
    current <- x0
    for (i in seq_len(n_ahead)) {
      current <- Re(object$A %*% current)
      predictions[, i] <- current
    }
  } else {
    # Mode-based evolution
    b <- solve(object$modes, x0)
    predictions <- matrix(0, nrow = length(x0), ncol = n_ahead)
    for (k in seq_len(n_ahead)) {
      evolution <- object$eigenvalues^k
      predictions[, k] <- Re(object$modes %*% (b * evolution))
    }
  }

  predictions
}
```

---

## Migration Path

### Step 1: Deprecate Old API
- Keep `getAMatrix()` and `predict.kdmd()` temporarily
- Add `.Deprecated()` calls pointing to new functions
- Document migration in vignette

### Step 2: Provide Wrappers
```r
getAMatrix <- function(...) {
  .Deprecated("dmd")
  result <- dmd(...)
  result$A
}
```

### Step 3: Remove in v1.0.0
- Plan removal for next major version
- Document in NEWS.md

---

## File Structure (Final)

```
RKoopmanDMD/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE.md
├── NEWS.md
├── README.md
├── .Rbuildignore
├── .gitignore
├── R/
│   ├── RKoopmanDMD-package.R    # Package documentation
│   ├── dmd.R                     # Main dmd() function
│   ├── predict.R                 # predict.dmd()
│   ├── methods.R                 # print, summary, plot, coef
│   ├── analysis.R                # spectrum, stability, reconstruct
│   ├── utils.R                   # Internal helpers
│   └── deprecated.R              # Old API with deprecation warnings
├── man/                          # Generated by roxygen2
├── data/
│   ├── lorenz.rda
│   └── oscillator.rda
├── inst/
│   └── CITATION
├── tests/
│   ├── testthat.R
│   └── testthat/
│       ├── test-dmd.R
│       ├── test-predict.R
│       └── ...
├── vignettes/
│   ├── introduction.Rmd
│   └── dynamical-systems.Rmd
└── .github/
    └── workflows/
        └── R-CMD-check.yaml
```

---

## Priority Order

1. **High Priority** (Core functionality)
   - Fix bugs in main.R (undefined variable, typo)
   - Implement unified `dmd()` function
   - Implement `predict.dmd()`
   - Update DESCRIPTION properly
   - Generate NAMESPACE with roxygen2

2. **Medium Priority** (User experience)
   - Add print/summary/plot methods
   - Write introduction vignette
   - Create test suite
   - Add example datasets

3. **Lower Priority** (Polish)
   - Advanced analysis functions
   - Additional vignettes
   - GitHub Actions CI
   - pkgdown site

---

## Key References to Include

1. Schmid, P. J. (2010). Dynamic mode decomposition of numerical and experimental data. Journal of Fluid Mechanics, 656, 5-28.

2. Kutz, J. N., Brunton, S. L., Brunton, B. W., & Proctor, J. L. (2016). Dynamic Mode Decomposition: Data-Driven Modeling of Complex Systems. SIAM.

3. Williams, M. O., Kevrekidis, I. G., & Rowley, C. W. (2015). A data-driven approximation of the Koopman operator: Extending dynamic mode decomposition. Journal of Nonlinear Science, 25(6), 1307-1346.

---

## Success Criteria

The package will be considered CRAN-ready when:

- [ ] `R CMD check` passes with 0 errors, 0 warnings, 0 notes
- [ ] All exported functions have complete documentation with examples
- [ ] Test coverage > 80%
- [ ] At least one vignette demonstrating typical workflow
- [ ] README with installation and quick-start guide
- [ ] Properly declared dependencies
- [ ] License file present
- [ ] Examples run in reasonable time (< 5s each)
- [ ] Package can be installed from source without errors

---

## Phase 6: Lifting Functions (Extended Koopman DMD)

### Overview

Lifting functions are a key feature for Koopman-based analysis. They transform the original state space into a higher-dimensional feature space where nonlinear dynamics become (approximately) linear. This enables DMD to capture nonlinear behavior that standard DMD cannot.

The mathematical idea:
- Original state: `x ∈ ℝⁿ`
- Lifted state: `g(x) ∈ ℝᵐ` where `m ≥ n`
- DMD operates on `g(x)` to find linear dynamics in lifted space
- Predictions are made in lifted space, then projected back

### 6.1 Design Goals

1. **User-friendly API**: Users can easily specify lifting functions
2. **Flexibility**: Support custom user-defined lifting functions
3. **Built-in options**: Provide common lifting functions out of the box
4. **Seamless integration**: Works transparently with existing `dmd()` and `predict()` API
5. **Inverse mapping**: Handle projection back to original state space

### 6.2 API Design

#### New Parameter in `dmd()`

```r
dmd <- function(X, rank = NULL, center = FALSE,
                lifting = NULL,          # NEW: lifting function specification
                observables = 1:nrow(X)) # NEW: which lifted dimensions are observables
```

**`lifting` parameter options:**

1. **NULL** (default): No lifting, standard DMD
2. **Character string**: Use a built-in lifting function
   - `"poly2"` - Polynomial degree 2 (includes x, x²)
   - `"poly3"` - Polynomial degree 3 (includes x, x², x³)
   - `"poly_cross2"` - Polynomial with cross-terms up to degree 2
   - `"trig"` - Trigonometric (includes x, sin(x), cos(x))
   - `"rbf"` - Radial basis functions (requires centers)
   - `"delay"` - Time-delay embedding
3. **Function**: Custom lifting function `function(X) -> X_lifted`
4. **List**: Specification with parameters

```r
# Examples of lifting parameter usage:

# Built-in polynomial lifting
model <- dmd(X, lifting = "poly2")

# Built-in with parameters
model <- dmd(X, lifting = list(type = "poly", degree = 4))
model <- dmd(X, lifting = list(type = "delay", delays = 3))
model <- dmd(X, lifting = list(type = "rbf", centers = centers_matrix, sigma = 0.5))

# Custom lifting function
my_lift <- function(X) {
  rbind(X, X^2, sin(X))
}
model <- dmd(X, lifting = my_lift)
```

**`observables` parameter:**

Specifies which rows of the lifted state correspond to the original observables. This is crucial for projecting predictions back to the original state space.

- Default: `1:nrow(X)` (assumes first `n` rows are original states)
- For custom lifting: user specifies which rows map back

### 6.3 Implementation Plan

#### 6.3.1 Create `R/lifting.R` - Lifting Function Infrastructure

```r
# Contents:
# - lift_data() - Apply lifting transformation to data matrix
# - unlift_data() - Project lifted predictions back to observable space
# - make_lifting_fn() - Create lifting function from specification
# - Built-in lifting functions:
#   - lift_polynomial()
#   - lift_poly_cross()
#   - lift_trigonometric()
#   - lift_delay()
#   - lift_rbf()
# - validate_lifting() - Validate lifting function specification
```

#### 6.3.2 Built-in Lifting Functions

**Polynomial Lifting (`lift_polynomial`)**
```r
lift_polynomial <- function(X, degree = 2, include_original = TRUE) {
  # X: n x T matrix
  # Returns: (n * degree) x T matrix
  # Contains: X, X^2, ..., X^degree

  lifted <- X
  for (d in 2:degree) {
    lifted <- rbind(lifted, X^d)
  }
  lifted
}
```

**Polynomial with Cross-terms (`lift_poly_cross`)**
```r
lift_poly_cross <- function(X, degree = 2) {
  # For n=2, degree=2: returns [x1, x2, x1^2, x1*x2, x2^2]
  # Uses polynomial expansion with all cross-terms
}
```

**Trigonometric Lifting (`lift_trigonometric`)**
```r
lift_trigonometric <- function(X, harmonics = 1) {
  # Returns: X, sin(X), cos(X), sin(2X), cos(2X), ...
  lifted <- X
  for (h in seq_len(harmonics)) {
    lifted <- rbind(lifted, sin(h * X), cos(h * X))
  }
  lifted
}
```

**Time-delay Embedding (`lift_delay`)**
```r
lift_delay <- function(X, delays = 2) {
  # Creates time-delay coordinates
  # Returns: [x(t), x(t-1), x(t-2), ...] for each time point
  # Note: This reduces the number of usable time points
}
```
**Radial Basis Functions (`lift_rbf`)**
```r
lift_rbf <- function(X, centers, sigma = 1) {
  # Gaussian RBF: phi_i(x) = exp(-||x - c_i||^2 / (2*sigma^2))
  # centers: matrix where each column is a center point
  # Returns: original X plus RBF features
}
```

#### 6.3.3 Modify `dmd()` in `R/dmd.R`

Add lifting transformation after input validation but before data splitting:

```r
dmd <- function(X, rank = NULL, center = FALSE,
                lifting = NULL, observables = NULL) {

  cl <- match.call()
  X <- validate_matrix(X, min_rows = 1, min_cols = 3)

  n_vars_original <- nrow(X)
  n_time <- ncol(X)

  # Store original data info
  X_original_first <- X[, 1]
  X_original_last <- X[, n_time]

  # Apply lifting transformation if specified
  lifting_fn <- NULL
  if (!is.null(lifting)) {
    lifting_fn <- make_lifting_fn(lifting, X)
    X <- lift_data(X, lifting_fn)

    # Update observables default if not specified
    if (is.null(observables)) {
      observables <- 1:n_vars_original
    }
  } else {
    observables <- 1:n_vars_original
  }

  n_vars_lifted <- nrow(X)

  # ... rest of DMD algorithm on lifted X ...

  # Store lifting info in result
  result <- structure(
    list(
      # ... existing fields ...
      lifting = lifting,
      lifting_fn = lifting_fn,
      observables = observables,
      n_vars_original = n_vars_original,
      n_vars_lifted = n_vars_lifted,
      X_original_first = X_original_first,
      X_original_last = X_original_last
    ),
    class = "dmd"
  )
}
```

#### 6.3.4 Modify `predict.dmd()` in `R/predict.R`

Handle lifted state in predictions:

```r
predict.dmd <- function(object, n_ahead = 10, x0 = NULL,
                        method = c("modes", "matrix"),
                        return_lifted = FALSE, ...) {

  # If x0 provided in original space, lift it
  if (!is.null(x0) && !is.null(object$lifting_fn)) {
    if (length(x0) == object$n_vars_original) {
      x0_matrix <- matrix(x0, ncol = 1)
      x0 <- as.vector(lift_data(x0_matrix, object$lifting_fn))
    }
  }

  # ... perform prediction in lifted space ...

  # Project back to observable space unless return_lifted = TRUE
  if (!return_lifted && !is.null(object$observables)) {
    predictions <- predictions[object$observables, , drop = FALSE]
  }

  predictions
}
```

### 6.4 Convenience Functions

#### `dmd_lift()` - Inspect Lifting Transformation

```r
#' Examine the lifting transformation
#' @param object A dmd object
#' @param x0 A state vector to transform (optional)
#' @return Information about the lifting or the lifted state
#' @export
dmd_lift <- function(object, x0 = NULL) {
  if (is.null(object$lifting)) {
    message("No lifting function was used in this model")
    return(invisible(NULL))
  }

  if (!is.null(x0)) {
    # Apply lifting to provided state
    x0_matrix <- matrix(x0, ncol = 1)
    return(as.vector(lift_data(x0_matrix, object$lifting_fn)))
  }

  # Return lifting info
  list(
    type = object$lifting,
    original_dim = object$n_vars_original,
    lifted_dim = object$n_vars_lifted,
    observables = object$observables
  )
}
```

### 6.5 Examples for Documentation

```r
# Example 1: Polynomial lifting for nonlinear oscillator
t <- seq(0, 20, by = 0.1)
x1 <- cos(t) + 0.5 * cos(t)^2  # Nonlinear component
x2 <- sin(t)
X <- rbind(x1, x2)

# Standard DMD (may struggle with nonlinearity)
model_standard <- dmd(X)

# DMD with polynomial lifting
model_poly <- dmd(X, lifting = "poly2")

# Compare predictions
pred_standard <- predict(model_standard, n_ahead = 50)
pred_poly <- predict(model_poly, n_ahead = 50)

# Example 2: Custom lifting function
my_lifting <- function(X) {
  # State: [x, y]
  # Lifted: [x, y, x^2, xy, y^2, x^3]
  x <- X[1, ]
  y <- X[2, ]
  rbind(X, x^2, x*y, y^2, x^3)
}

model_custom <- dmd(X, lifting = my_lifting, observables = 1:2)

# Example 3: Time-delay embedding for scalar time series
y <- sin(seq(0, 10, by = 0.1)) + 0.3 * sin(3 * seq(0, 10, by = 0.1))
Y <- matrix(y, nrow = 1)

model_delay <- dmd(Y, lifting = list(type = "delay", delays = 5))

# Example 4: RBF lifting for complex dynamics
centers <- matrix(rnorm(20), nrow = 2, ncol = 10)  # 10 RBF centers
model_rbf <- dmd(X, lifting = list(type = "rbf", centers = centers, sigma = 0.5))
```

### 6.6 Testing Plan

Add `tests/testthat/test-lifting.R`:

```r
# Test categories:
# 1. Built-in lifting functions produce correct dimensions
# 2. Custom lifting functions are validated properly
# 3. Predictions are correctly projected back to original space
# 4. Lifting improves accuracy on known nonlinear systems
# 5. Edge cases: lifting to very high dimensions, singular cases
# 6. observables parameter correctly identifies original states
```

### 6.7 Documentation Updates

1. **Update `R/dmd.R`** - Add `@param` for `lifting` and `observables`
2. **Create `R/lifting.R`** - Full roxygen documentation for all lifting functions
3. **Update vignettes** - Add section on lifting functions
4. **Create new vignette** - `vignettes/lifting-functions.Rmd` covering:
   - When to use lifting functions
   - Choosing the right lifting function
   - Custom lifting functions
   - Performance considerations

### 6.8 Priority and Dependencies

**Dependencies:**
- Phase 1-4 must be complete (core functionality working)
- Good test coverage for base `dmd()` before modifying

**Priority within Phase 6:**
1. Core `lifting.R` infrastructure and `make_lifting_fn()`
2. Polynomial lifting (most common use case)
3. Time-delay embedding (popular for scalar time series)
4. Modify `dmd()` to accept lifting parameter
5. Modify `predict.dmd()` to handle lifting
6. RBF and trigonometric lifting
7. Documentation and vignettes
8. Comprehensive tests

### 6.9 Future Extensions

Potential future enhancements (post-v1.0):

1. **Automatic lifting selection** - Cross-validation to choose best lifting
2. **Sparse lifting** - L1-regularized selection of lifting functions
3. **Neural network lifting** - Learned lifting via autoencoder (requires torch)
4. **Kernel DMD** - Implicit infinite-dimensional lifting via kernel trick
5. **Dictionary learning** - Adaptive basis selection from data
