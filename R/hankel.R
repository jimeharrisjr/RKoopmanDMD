#' Hankel-DMD (Krylov Subspace Method)
#'
#' Performs Dynamic Mode Decomposition using time-delayed observables
#' forming a Krylov subspace. This method avoids the curse of dimensionality
#' that affects standard EDMD with polynomial/RBF basis functions.
#'
#' @param y Numeric vector or matrix. For vector: a scalar time series.
#'   For matrix: rows are observables, columns are time points.
#' @param delays Integer; number of time delays to use (embedding dimension - 1).
#'   If NULL (default), automatically selected based on data length.
#' @param rank Integer; rank for SVD truncation. If NULL (default),
#'   automatically selected to capture 99% of variance.
#' @param dt Numeric; time step between observations. Default is 1.
#'
#' @return An object of class `c("hankel_dmd", "dmd")` containing:
#' \describe{
#'   \item{A}{The DMD matrix (Koopman operator approximation)}
#'   \item{modes}{Complex matrix of DMD modes}
#'   \item{eigenvalues}{Complex vector of DMD eigenvalues}
#'   \item{amplitudes}{Complex vector of initial mode amplitudes}
#'   \item{rank}{Integer; the rank used for truncation}
#'   \item{delays}{Number of delays used in Hankel matrix}
#'   \item{hankel}{The Hankel-Takens matrix}
#'   \item{companion}{The companion matrix (if computed)}
#'   \item{residual}{Residual norm for error assessment}
#'   \item{svd}{List containing truncated SVD components}
#'   \item{y_original}{Original time series data}
#'   \item{dt}{Time step}
#' }
#'
#' @details
#' The Hankel-DMD method constructs a Hankel (Takens) matrix from time-delayed
#' copies of the observable:
#' \deqn{H = \begin{bmatrix} y_1 & y_2 & \cdots & y_n \\
#'                           y_2 & y_3 & \cdots & y_{n+1} \\
#'                           \vdots & & & \vdots \\
#'                           y_m & y_{m+1} & \cdots & y_{m+n-1}
#'       \end{bmatrix}}
#'
#' This is equivalent to using the Krylov sequence \eqn{(f, Uf, U^2f, \ldots)}
#' as the basis for the finite section approximation.
#'
#' @section Advantages over Standard DMD:
#' \itemize{
#'   \item Avoids curse of dimensionality: basis size equals number of delays,
#'     not exponential in state dimension
#'   \item Natural choice of observables selected by the dynamics
#'   \item Proven pseudospectral convergence for systems with pure point spectrum
#'   \item Connected to Takens' embedding theorem for attractor reconstruction
#' }
#'
#' @section Theoretical Background:
#' This implementation is based on Section 5 of Mezic (2020), which shows that
#' Krylov subspace methods provide pseudospectral convergence without requiring
#' an exponentially large basis. The companion matrix representation (eq. 83)
#' captures the dynamics in the time-delay coordinate system.
#'
#' @examples
#' # Scalar time series from oscillator
#' t <- seq(0, 20, by = 0.1)
#' y <- cos(2 * pi * 0.1 * t) + 0.5 * sin(2 * pi * 0.25 * t)
#'
#' # Fit Hankel-DMD
#' model <- hankel_dmd(y, delays = 10)
#' print(model)
#'
#' # Check eigenvalues (should reveal frequencies)
#' spectrum <- dmd_spectrum(model, dt = 0.1)
#' print(spectrum)
#'
#' # Predict future values
#' future <- predict(model, n_ahead = 50)
#'
#' # Multivariate time series
#' X <- rbind(cos(t), sin(t))
#' model_mv <- hankel_dmd(X, delays = 5)
#'
#' @references
#' Mezic, I. (2020). On Numerical Approximations of the Koopman Operator.
#' arXiv:2009.05883, Section 5.
#'
#' Arbabi, H. and Mezic, I. (2017). Ergodic theory, dynamic mode decomposition,
#' and computation of spectral properties of the Koopman operator.
#' SIAM Journal on Applied Dynamical Systems, 16(4):2096-2126.
#'
#' Takens, F. (1981). Detecting strange attractors in turbulence.
#' In Dynamical Systems and Turbulence, pages 366-381. Springer.
#'
#' @seealso [dmd()] for standard DMD, [dmd_spectrum()] for eigenvalue analysis,
#'   [predict.dmd()] for forecasting.
#'
#' @export
hankel_dmd <- function(y, delays = NULL, rank = NULL, dt = 1) {


  cl <- match.call()


  # Handle input: convert vector to 1-row matrix

  if (is.vector(y)) {
    y <- matrix(y, nrow = 1)
  }

  if (is.data.frame(y)) {
    y <- as.matrix(y)
  }

  if (!is.matrix(y) || !is.numeric(y)) {
    stop("y must be a numeric vector or matrix", call. = FALSE)
  }

  if (any(!is.finite(y))) {
    stop("y contains non-finite values (NA, NaN, or Inf)", call. = FALSE)
  }

  n_obs <- nrow(y)  # number of observables

  n_time <- ncol(y)  # number of time points

  # Determine number of delays
  if (is.null(delays)) {
    # Heuristic: use roughly sqrt(n_time) delays, capped
    delays <- min(floor(sqrt(n_time)), floor(n_time / 3), 50)
    delays <- max(delays, 2)  # at least 2 delays
  }

  delays <- as.integer(delays)
  if (delays < 1) {
    stop("delays must be at least 1", call. = FALSE)
  }

  if (delays >= n_time - 1) {
    stop(sprintf("delays (%d) must be less than n_time - 1 (%d)",
                 delays, n_time - 1), call. = FALSE)
  }

  # Build Hankel-Takens matrix (eq. 118 from paper)
  # H[i,j] = y[:, i+j-1] for appropriate ranges
  # Rows: embedding dimension (delays + 1) * n_obs
  # Cols: n_time - delays
  H <- build_hankel_matrix(y, delays)

  n_rows_H <- nrow(H)
  n_cols_H <- ncol(H)

  # Split into H1 (all but last column) and H2 (all but first column)
  H1 <- H[, -n_cols_H, drop = FALSE]
  H2 <- H[, -1, drop = FALSE]

  # SVD of H1
  svd_result <- svd(H1)

  # Determine rank
  r <- determine_rank(svd_result$d, rank)

  # Truncate SVD
  U_r <- svd_result$u[, 1:r, drop = FALSE]
  S_r <- svd_result$d[1:r]
  V_r <- svd_result$v[, 1:r, drop = FALSE]

  # Compute reduced DMD matrix: A_tilde = U^T H2 V S^{-1}
  S_inv <- diag(1 / S_r, nrow = r, ncol = r)
  A_tilde <- t(U_r) %*% H2 %*% V_r %*% S_inv

  # Eigendecomposition of reduced matrix
  eig <- eigen(A_tilde)
  eigenvalues <- eig$values
  W <- eig$vectors

  # DMD modes in Hankel space: Phi = H2 V S^{-1} W
  Phi <- H2 %*% V_r %*% S_inv %*% W

  # Full A matrix in Hankel coordinates
  A <- Phi %*% diag(eigenvalues, nrow = r, ncol = r) %*% safe_solve(Phi, diag(n_rows_H))

  # Initial amplitudes: project first Hankel column onto modes
  b <- safe_solve(Phi, H[, 1])

  # Compute residual for error assessment (eq. 101)
  # residual = H2[:, end] - A %*% H1[:, end]
  H_last <- H[, n_cols_H]
  H_pred <- Re(A %*% H[, n_cols_H - 1])
  residual_vec <- H_last - H_pred
  residual_norm <- sqrt(sum(residual_vec^2))

  # Also compute companion matrix coefficients if rank allows
  companion <- NULL
  if (r >= n_rows_H) {
    # Can compute full companion matrix
    companion <- compute_companion_matrix(H)
  }

  # Store first/last states for prediction
  H_first <- H[, 1]
  H_last_state <- H[, n_cols_H]

  # Construct result object
  result <- structure(
    list(
      A = Re(A),
      modes = Phi,
      eigenvalues = eigenvalues,
      amplitudes = b,
      rank = r,
      delays = delays,
      hankel = H,
      companion = companion,
      residual = residual_norm,
      residual_vector = residual_vec,
      svd = list(U = U_r, S = S_r, V = V_r),
      A_tilde = A_tilde,
      X_first = H_first,
      X_last = H_last_state,
      data_dim = c(n_rows_H, n_cols_H),
      center = FALSE,
      X_mean = NULL,
      dt = dt,
      call = cl,
      # Original data info
      y_original = y,
      n_obs = n_obs,
      n_time = n_time,
      # For compatibility with dmd methods
      lifting = NULL,
      lifting_fn = NULL,
      observables = seq_len(n_obs),
      n_vars_original = n_obs,
      n_vars_lifted = n_rows_H
    ),
    class = c("hankel_dmd", "dmd")
  )

  result
}


#' Build Hankel-Takens Matrix
#'
#' Constructs the Hankel (Takens) matrix from time series data for
#' delay embedding.
#'
#' @param y Numeric matrix (n_obs x n_time).
#' @param delays Number of delay steps.
#'
#' @return Hankel matrix of dimension ((delays + 1) * n_obs) x (n_time - delays).
#'
#' @details
#' The Hankel-Takens matrix has the structure (eq. 118 from Mezic 2020):
#' \deqn{H_{i,j} = y(:, i+j-1)}
#'
#' For a scalar time series y = [y1, y2, ..., yn], with d delays:
#' \deqn{H = \begin{bmatrix}
#'   y_1 & y_2 & \cdots & y_{n-d} \\
#'   y_2 & y_3 & \cdots & y_{n-d+1} \\
#'   \vdots & & & \vdots \\
#'   y_{d+1} & y_{d+2} & \cdots & y_n
#' \end{bmatrix}}
#'
#' @keywords internal
build_hankel_matrix <- function(y, delays) {

  n_obs <- nrow(y)
  n_time <- ncol(y)

  n_rows <- (delays + 1) * n_obs
  n_cols <- n_time - delays

  H <- matrix(0, nrow = n_rows, ncol = n_cols)

  for (d in 0:delays) {
    row_start <- d * n_obs + 1
    row_end <- (d + 1) * n_obs
    col_indices <- (d + 1):(d + n_cols)
    H[row_start:row_end, ] <- y[, col_indices, drop = FALSE]
  }

  # Add row names for interpretability
  if (!is.null(rownames(y))) {
    base_names <- rownames(y)
  } else {
    base_names <- paste0("y", seq_len(n_obs))
  }

  row_names <- base_names
  for (d in seq_len(delays)) {
    row_names <- c(row_names, paste0(base_names, "(t-", d, ")"))
  }
  rownames(H) <- row_names

  H
}


#' Compute Companion Matrix
#'
#' Computes the companion matrix representation for Krylov subspace DMD.
#'
#' @param H Hankel matrix.
#'
#' @return Companion matrix C such that H2 ≈ H1 C.
#'
#' @details
#' The companion matrix has the structure (eq. 83):
#' \deqn{C = \begin{bmatrix}
#'   0 & 0 & \cdots & 0 & c_1 \\
#'   1 & 0 & \cdots & 0 & c_2 \\
#'   0 & 1 & \cdots & 0 & c_3 \\
#'   \vdots & & \ddots & & \vdots \\
#'   0 & 0 & \cdots & 1 & c_n
#' \end{bmatrix}}
#'
#' @keywords internal
compute_companion_matrix <- function(H) {

  n_cols <- ncol(H)
  n_rows <- nrow(H)

  H1 <- H[, -n_cols, drop = FALSE]
  h_next <- H[, n_cols]

  # Solve for companion coefficients: H1 %*% c = h_next
  # This is overdetermined if n_rows > n_cols - 1
  c_coef <- safe_solve(t(H1) %*% H1, t(H1) %*% h_next)

  # Build companion matrix
  n <- length(c_coef)
  C <- matrix(0, nrow = n, ncol = n)

  # Subdiagonal of ones
  if (n > 1) {
    for (i in 2:n) {
      C[i, i - 1] <- 1
    }
  }

  # Last column is c_coef

  C[, n] <- c_coef

  C
}


#' Print Method for Hankel-DMD Objects
#'
#' @param x A `"hankel_dmd"` object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.hankel_dmd <- function(x, ...) {
  cat("Hankel-DMD (Krylov Subspace) Model\n")
  cat("==================================\n")
  cat(sprintf("  Original observables: %d\n", x$n_obs))
  cat(sprintf("  Time points:          %d\n", x$n_time))
  cat(sprintf("  Delays used:          %d\n", x$delays))
  cat(sprintf("  Hankel dimension:     %d x %d\n", nrow(x$hankel), ncol(x$hankel)))
  cat(sprintf("  Rank used:            %d\n", x$rank))
  cat(sprintf("  Residual norm:        %.6g\n", x$residual))

  # Stability check

  mags <- Mod(x$eigenvalues)
  tol <- 1e-6
  if (all(mags < 1 - tol)) {
    cat("  Stability:            stable (all |eigenvalues| < 1)\n")
  } else if (any(mags > 1 + tol)) {
    n_unstable <- sum(mags > 1 + tol)
    cat(sprintf("  Stability:            unstable (%d mode(s) with |eigenvalue| > 1)\n",
                n_unstable))
  } else {
    cat("  Stability:            marginally stable (eigenvalues near unit circle)\n")
  }

  cat("\nUse summary() for detailed eigenvalue analysis.\n")
  cat("Use dmd_spectrum() with dt parameter for frequency analysis.\n")

  invisible(x)
}


#' Predict Method for Hankel-DMD Objects
#'
#' Forecasts future values using the Hankel-DMD model. Predictions are
#' made in the delay-embedded space and then extracted back to the
#' original observable space.
#'
#' @param object A `"hankel_dmd"` object.
#' @param n_ahead Integer; number of time steps to forecast.
#' @param x0 Initial state. If NULL, uses the last observed Hankel state.
#'   Can be provided as either:
#'   \itemize{
#'     \item A vector of length (delays + 1) * n_obs for the full Hankel state

#'     \item A vector of length n_obs for just the current observation
#'       (delay history taken from training data)
#'   }
#' @param method Prediction method: "modes" (default) or "matrix".
#' @param return_full Logical; if TRUE, return predictions in full Hankel
#'   space. If FALSE (default), extract only the first n_obs rows
#'   (current time predictions).
#' @param ... Additional arguments (unused).
#'
#' @return Matrix of predictions. If return_full = FALSE, dimensions are
#'   (n_obs x n_ahead). If return_full = TRUE, dimensions are
#'   ((delays + 1) * n_obs x n_ahead).
#'
#' @export
predict.hankel_dmd <- function(object, n_ahead = 10, x0 = NULL,
                                method = c("modes", "matrix"),
                                return_full = FALSE, ...) {

  method <- match.arg(method)
  n_ahead <- as.integer(n_ahead)

  if (n_ahead < 1) {
    stop("n_ahead must be a positive integer", call. = FALSE)
  }

  n_obs <- object$n_obs
  n_hankel <- nrow(object$hankel)

  # Determine initial condition in Hankel space

  if (is.null(x0)) {
    x0_hankel <- object$X_last
  } else {
    x0 <- as.numeric(x0)
    if (length(x0) == n_hankel) {
      # Full Hankel state provided
      x0_hankel <- x0
    } else if (length(x0) == n_obs) {
      # Only current observation - construct Hankel state
      # Use last delays from training data
      H <- object$hankel
      x0_hankel <- c(x0, H[(n_obs + 1):n_hankel, ncol(H)])
    } else {
      stop(sprintf("x0 must have length %d (current obs) or %d (full Hankel state)",
                   n_obs, n_hankel), call. = FALSE)
    }
  }

  # Predict in Hankel space
  if (method == "matrix") {
    predictions <- predict_matrix_hankel(object, x0_hankel, n_ahead)
  } else {
    predictions <- predict_modes_hankel(object, x0_hankel, n_ahead)
  }

  # Extract original observables (first n_obs rows)
  if (!return_full) {
    predictions <- predictions[1:n_obs, , drop = FALSE]
  }

  colnames(predictions) <- paste0("t+", seq_len(n_ahead))

  predictions
}


#' Matrix-based prediction for Hankel-DMD (internal)
#' @keywords internal
predict_matrix_hankel <- function(object, x0, n_ahead) {
  n_vars <- length(x0)
  predictions <- matrix(0, nrow = n_vars, ncol = n_ahead)

  current <- x0
  A <- object$A

  for (i in seq_len(n_ahead)) {
    current <- Re(A %*% current)
    predictions[, i] <- current
  }

  predictions
}


#' Mode-based prediction for Hankel-DMD (internal)
#' @keywords internal
predict_modes_hankel <- function(object, x0, n_ahead) {
  n_vars <- length(x0)
  predictions <- matrix(0, nrow = n_vars, ncol = n_ahead)

  modes <- object$modes
  lambdas <- object$eigenvalues

  # Project initial condition onto modes
  b <- safe_solve(modes, x0)

  for (k in seq_len(n_ahead)) {
    evolution <- lambdas^k
    predictions[, k] <- Re(modes %*% (b * evolution))
  }

  predictions
}


#' Reconstruct Time Series from Hankel-DMD
#'
#' Reconstructs the original time series using the Hankel-DMD model.
#'
#' @param object A `"hankel_dmd"` object.
#' @param n_steps Number of steps to reconstruct. Default is original length.
#'
#' @return Reconstructed time series matrix (n_obs x n_steps).
#'
#' @export
hankel_reconstruct <- function(object, n_steps = NULL) {

  if (!inherits(object, "hankel_dmd"))
 {
    stop("object must be a hankel_dmd object", call. = FALSE)
  }

  if (is.null(n_steps)) {
    n_steps <- ncol(object$hankel)
  }

  n_obs <- object$n_obs
  n_hankel <- nrow(object$hankel)

  # Reconstruct Hankel states
  modes <- object$modes
  lambdas <- object$eigenvalues
  b <- object$amplitudes

  reconstruction <- matrix(0, nrow = n_obs, ncol = n_steps)

  for (k in seq_len(n_steps)) {
    evolution <- lambdas^(k - 1)
    hankel_state <- Re(modes %*% (b * evolution))
    # Extract first n_obs elements (current time values)
    reconstruction[, k] <- hankel_state[1:n_obs]
  }

  colnames(reconstruction) <- paste0("t", seq_len(n_steps))

  reconstruction
}
