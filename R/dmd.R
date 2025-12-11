#' Dynamic Mode Decomposition
#'
#' Performs Dynamic Mode Decomposition (DMD) on time-series data to extract
#' spatiotemporal modes, eigenvalues, and a linear operator approximating
#' the system dynamics.
#'
#' @param X A numeric matrix where rows represent state variables and columns
#'   represent time snapshots. Must have at least 2 rows and 3 columns.
#' @param rank Integer specifying the rank for SVD truncation. If `NULL`
#'   (default), rank is automatically selected to capture 99% of variance.
#' @param center Logical; if `TRUE`, the data is centered by subtracting the
#'   row means before analysis. Default is `FALSE`.
#' @param lifting Lifting function specification for extended DMD. Can be:
#'   \describe{
#'     \item{NULL}{(default) No lifting, standard DMD}
#'     \item{character}{Built-in lifting: `"poly2"`, `"poly3"`, `"poly4"`,
#'       `"poly_cross2"`, `"poly_cross3"`, `"trig"`, `"trig2"`,
#'       `"delay2"`, `"delay3"`, `"delay5"`}
#'     \item{function}{Custom lifting function taking matrix X and returning
#'       lifted matrix with same number of columns}
#'     \item{list}{Parametric specification, e.g.,
#'       `list(type = "poly", degree = 4)` or
#'       `list(type = "rbf", centers = centers_matrix, sigma = 0.5)`}
#'   }
#'   See [list_lifting_functions()] for available options.
#' @param observables Integer vector specifying which rows of the lifted state
#'   correspond to the original observable states. Default is `1:nrow(X)`,
#'   assuming the first n rows are the original states. Used for projecting
#'   predictions back to observable space.
#'
#' @return An object of class `"dmd"` containing:
#' \describe{
#'   \item{A}{The DMD matrix (Koopman operator approximation) in the original
#'     state space. Multiplying a state by A gives the predicted next state.}
#'   \item{modes}{Complex matrix of DMD modes (columns). Each mode represents
#'     a spatial pattern that evolves with a single frequency and growth rate.}
#'   \item{eigenvalues}{Complex vector of DMD eigenvalues. The magnitude
#'     indicates growth/decay rate; the angle indicates oscillation frequency.}
#'   \item{amplitudes}{Complex vector of initial mode amplitudes, computed
#'     by projecting the first snapshot onto the modes.}
#'   \item{rank}{Integer; the rank used for truncation.}
#'   \item{svd}{List containing truncated SVD components (U, S, V).}
#'   \item{A_tilde}{The reduced DMD matrix in the SVD coordinate system.}
#'   \item{X_first}{First column of input data (initial condition).}
#'   \item{X_last}{Last column of input data (for prediction initialization).}
#'   \item{data_dim}{Dimensions of the input matrix (rows, columns).}
#'   \item{center}{Logical; whether centering was applied.}
#'   \item{X_mean}{Row means if centering was applied, otherwise NULL.}
#'   \item{dt}{Time step (currently set to 1; for future extension).}
#'   \item{call}{The matched function call.}
#'   \item{lifting}{The lifting specification used (NULL if none).}
#'   \item{lifting_fn}{The actual lifting function (NULL if none).}
#'   \item{observables}{Indices of observable states in lifted space.}
#'   \item{n_vars_original}{Number of original state variables.}
#'   \item{n_vars_lifted}{Number of lifted state variables.}
#' }
#'
#' @details
#' DMD approximates the dynamics of a system by finding the best-fit linear
#' operator \eqn{A} such that \eqn{X_2 \approx A X_1}, where \eqn{X_1} and
#' \eqn{X_2} are time-shifted versions of the data.
#'
#' The algorithm:
#' 1. Splits data into \eqn{X_1 = X[, 1:(n-1)]} and \eqn{X_2 = X[, 2:n]}
#' 2. Computes truncated SVD of \eqn{X_1}: \eqn{X_1 = U \Sigma V^T}
#' 3. Projects onto reduced coordinates: \eqn{\tilde{A} = U^T X_2 V \Sigma^{-1}}
#' 4. Computes eigendecomposition of \eqn{\tilde{A}}
#' 5. Recovers full-space modes: \eqn{\Phi = X_2 V \Sigma^{-1} W}
#'
#' @section Interpreting Results:
#' - **Eigenvalue magnitude < 1**: Mode is decaying (stable)
#' - **Eigenvalue magnitude = 1**: Mode is neutral (marginally stable)
#' - **Eigenvalue magnitude > 1**: Mode is growing (unstable)
#' - **Eigenvalue phase**: Determines oscillation frequency as \eqn{\omega = \arg(\lambda) / \Delta t}
#'
#' @examples
#' # Create a simple oscillatory system
#' t <- seq(0, 10, by = 0.1)
#' x1 <- cos(2 * pi * 0.5 * t)
#' x2 <- sin(2 * pi * 0.5 * t)
#' X <- rbind(x1, x2)
#'
#' # Fit DMD model
#' model <- dmd(X)
#'
#' # Examine the model
#' print(model)
#' summary(model)
#'
#' # Check eigenvalue magnitudes (should be ~1 for oscillatory system)
#' Mod(model$eigenvalues)
#'
#' # Example with lifting for nonlinear dynamics
#' t <- seq(0, 10, by = 0.1)
#' x1 <- cos(t) + 0.3 * cos(t)^2  # Nonlinear component
#' x2 <- sin(t)
#' X_nonlin <- rbind(x1, x2)
#'
#' # Standard DMD
#' model_std <- dmd(X_nonlin)
#'
#' # DMD with polynomial lifting
#' model_lift <- dmd(X_nonlin, lifting = "poly2")
#'
#' # Check lifted dimensions
#' model_lift$n_vars_original  # 2
#' model_lift$n_vars_lifted    # 4 (x1, x2, x1^2, x2^2)
#'
#' @references
#' Schmid, P. J. (2010). Dynamic mode decomposition of numerical and
#' experimental data. Journal of Fluid Mechanics, 656, 5-28.
#'
#' @seealso [predict.dmd()] for forecasting, [dmd_spectrum()] for eigenvalue
#'   analysis, [dmd_reconstruct()] for data reconstruction,
#'   [dmd_lift()] for inspecting lifting transformations,
#'   [list_lifting_functions()] for available lifting options.
#'
#' @export
dmd <- function(X, rank = NULL, center = FALSE, lifting = NULL, observables = NULL) {

  # Capture the call
  cl <- match.call()

  # Validate input
  X <- validate_matrix(X, min_rows = 1, min_cols = 3)

  n_vars_original <- nrow(X)
  n_time_original <- ncol(X)

  # Store original endpoints before any transformation
  X_original_first <- X[, 1]
  X_original_last <- X[, n_time_original]

  # Apply lifting transformation if specified
  lifting_fn <- NULL
  if (!is.null(lifting)) {
    lifting_fn <- make_lifting_fn(lifting, X)
    # Allow column reduction for delay-based lifting
    allow_col_reduction <- is_delay_lifting(lifting)
    X <- lift_data(X, lifting_fn, allow_col_reduction = allow_col_reduction)

    # Set default observables if not specified
    if (is.null(observables)) {
      observables <- seq_len(n_vars_original)
    }
  } else {
    # No lifting - observables are all rows
    observables <- seq_len(n_vars_original)
  }

  n_vars <- nrow(X)
  n_time <- ncol(X)

  # Optionally center the data (after lifting)
  X_mean <- NULL
  if (center) {
    X_mean <- rowMeans(X)
    X <- X - X_mean
  }

  # Split into current and next state
  X1 <- X[, -n_time, drop = FALSE]
  X2 <- X[, -1, drop = FALSE]

  # SVD of X1
  svd_result <- svd(X1)

  # Determine rank (auto or specified)
  r <- determine_rank(svd_result$d, rank)

  # Truncate SVD
  U_r <- svd_result$u[, 1:r, drop = FALSE]
  S_r <- svd_result$d[1:r]
  V_r <- svd_result$v[, 1:r, drop = FALSE]

  # Compute reduced DMD matrix: A_tilde = U^T X2 V S^{-1}
  S_inv <- diag(1 / S_r, nrow = r, ncol = r)
  A_tilde <- t(U_r) %*% X2 %*% V_r %*% S_inv

  # Eigendecomposition of reduced matrix
  eig <- eigen(A_tilde)
  eigenvalues <- eig$values
  W <- eig$vectors

  # DMD modes in physical space: Phi = X2 V S^{-1} W
  Phi <- X2 %*% V_r %*% S_inv %*% W

  # Full A matrix in original coordinates
  # A = Phi * diag(lambda) * Phi^{-1}
  A <- Phi %*% diag(eigenvalues, nrow = r, ncol = r) %*% safe_solve(Phi, diag(n_vars))

  # Initial amplitudes: project first snapshot onto modes
  b <- safe_solve(Phi, X[, 1])

  # Construct result object
  result <- structure(
    list(
      A = Re(A),
      modes = Phi,
      eigenvalues = eigenvalues,
      amplitudes = b,
      rank = r,
      svd = list(U = U_r, S = S_r, V = V_r),
      A_tilde = A_tilde,
      X_first = X[, 1] + if (center) X_mean else 0,
      X_last = X[, n_time] + if (center) X_mean else 0,
      data_dim = c(n_vars, n_time),
      center = center,
      X_mean = X_mean,
      dt = 1,
      call = cl,
      # Lifting-related fields
      lifting = lifting,
      lifting_fn = lifting_fn,
      observables = observables,
      n_vars_original = n_vars_original,
      n_vars_lifted = n_vars,
      X_original_first = X_original_first,
      X_original_last = X_original_last
    ),
    class = "dmd"
  )

  result
}
