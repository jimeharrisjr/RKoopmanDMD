# Internal utility functions for RKoopmanDMD
# These are not exported


#' Null-coalescing operator
#'
#' Returns the left-hand side if not NULL, otherwise the right-hand side.
#'
#' @param x Left-hand side value
#' @param y Right-hand side value (default)
#'
#' @return x if not NULL, otherwise y
#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Validate input data matrix
#'
#' Checks that input is a valid matrix for DMD analysis.
#'
#' @param X Input data (matrix or data.frame)
#' @param min_rows Minimum number of rows required
#' @param min_cols Minimum number of columns required
#'
#' @return A numeric matrix
#' @keywords internal
validate_matrix <- function(X, min_rows = 2, min_cols = 3) {
  # Convert data.frame to matrix

if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  # Check if matrix
  if (!is.matrix(X)) {
    stop("X must be a matrix or data.frame", call. = FALSE)
  }

  # Check numeric
  if (!is.numeric(X)) {
    stop("X must contain numeric values", call. = FALSE)
  }

  # Check dimensions
  if (nrow(X) < min_rows) {
    stop(sprintf("X must have at least %d rows (state variables)", min_rows),
         call. = FALSE)
  }

  if (ncol(X) < min_cols) {
    stop(sprintf("X must have at least %d columns (time snapshots)", min_cols),
         call. = FALSE)
  }

  # Check for NA/NaN/Inf
  if (any(!is.finite(X))) {
    stop("X contains non-finite values (NA, NaN, or Inf)", call. = FALSE)
  }

  X
}


#' Determine optimal rank for SVD truncation
#'
#' Automatically selects the number of singular values to retain based on
#' the specified criteria.
#'
#' @param singular_values Vector of singular values from SVD
#' @param rank User-specified rank (NULL for auto, or integer)
#' @param threshold Variance threshold for auto selection (default 0.99)
#' @param tol Tolerance for near-zero singular values
#'
#' @return Integer rank value
#' @keywords internal
determine_rank <- function(singular_values, rank = NULL, threshold = 0.99,
                           tol = .Machine$double.eps * max(singular_values) * 100) {

  n <- length(singular_values)

  # Remove near-zero singular values
  valid_sv <- singular_values > tol
  max_rank <- sum(valid_sv)

  if (max_rank == 0) {
    stop("All singular values are near zero; data may be constant or degenerate",
         call. = FALSE)
  }

  if (!is.null(rank)) {
    # User-specified rank
    r <- as.integer(rank)
    if (r < 1) {
      warning("rank must be at least 1; setting to 1", call. = FALSE)
      r <- 1
    }
    if (r > max_rank) {
      warning(sprintf("Requested rank %d exceeds effective rank %d; using %d",
                      r, max_rank, max_rank), call. = FALSE)
      r <- max_rank
    }
  } else {
    # Auto-select based on variance explained
    variance_explained <- cumsum(singular_values^2) / sum(singular_values^2)
    r <- which(variance_explained >= threshold)[1]
    if (is.na(r)) {
      r <- max_rank
    }
    r <- min(r, max_rank)
  }

  r
}


#' Compute pseudo-inverse using SVD
#'
#' Computes the Moore-Penrose pseudo-inverse of a matrix.
#' Used as a fallback when solve() fails.
#'
#' @param A Input matrix
#' @param tol Tolerance for singular values
#'
#' @return Pseudo-inverse of A
#' @keywords internal
pinv <- function(A, tol = .Machine$double.eps * max(dim(A)) * max(abs(A))) {
  s <- svd(A)
  # Invert only non-zero singular values
  d_inv <- ifelse(s$d > tol, 1 / s$d, 0)
  s$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(s$u)
}


#' Safe matrix solve with pseudo-inverse fallback
#'
#' Attempts to solve a linear system, falling back to pseudo-inverse
#' if the system is singular or near-singular.
#'
#' @param A Coefficient matrix
#' @param b Right-hand side vector or matrix
#'
#' @return Solution x such that A %*% x ≈ b
#' @keywords internal
safe_solve <- function(A, b) {
  tryCatch(
    solve(A, b),
    error = function(e) {
      pinv(A) %*% b
    }
  )
}


#' Check if eigenvalues indicate stability
#'
#' For discrete-time systems, eigenvalues inside the unit circle indicate
#' stability.
#'
#' @param eigenvalues Complex vector of eigenvalues
#' @param tol Tolerance for boundary (default 1e-10)
#'
#' @return List with stability classification
#' @keywords internal
check_stability <- function(eigenvalues, tol = 1e-10) {
  magnitudes <- Mod(eigenvalues)

  list(
    stable = all(magnitudes < 1 - tol),
    marginal = any(abs(magnitudes - 1) < tol),
    unstable = any(magnitudes > 1 + tol),
    max_magnitude = max(magnitudes),
    n_unstable = sum(magnitudes > 1 + tol)
  )
}
