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
