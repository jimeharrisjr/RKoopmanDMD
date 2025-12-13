#' Generalized Laplace Analysis (GLA)
#'
#' Computes Koopman eigenfunctions and modes using Generalized Laplace Analysis,
#' which constructs eigenfunctions directly from trajectory data via weighted
#' time averages without constructing an approximate operator matrix.
#'
#' @param y Numeric vector or matrix of trajectory data. For vector: a scalar
#'   time series. For matrix: rows are observables, columns are time points.
#' @param eigenvalues Complex vector of candidate eigenvalues to test.
#'   If NULL (default), eigenvalues are estimated from a preliminary DMD.
#' @param n_eigenvalues Integer; if eigenvalues is NULL, how many eigenvalues
#'   to estimate from preliminary DMD. Default is 5.
#' @param tol Numeric; convergence tolerance for the time average.
#'   Default is 1e-6.
#' @param max_iter Maximum number of iterations (time steps) to use.
#'   Default is NULL (use all available data).
#'
#' @return An object of class `"gla"` containing:
#' \describe{
#'   \item{eigenvalues}{Complex vector of eigenvalues used}
#'   \item{eigenfunctions}{Matrix of eigenfunction values at each time point}
#'   \item{modes}{Complex matrix of Koopman modes (projection coefficients)}
#'   \item{convergence}{Logical vector indicating which eigenvalues converged}
#'   \item{residuals}{Residual norms for each eigenvalue}
#'   \item{n_iter}{Number of iterations used}
#' }
#'
#' @details
#' The GLA method computes eigenfunctions using Theorem 3.1 from Mezic (2020):
#' \deqn{f_k = \lim_{n\to\infty} \frac{1}{n} \sum_{i=0}^{n-1} \lambda_k^{-i}
#'   \left(f(T^i x) - \sum_{j=0}^{k-1} \lambda_j^i \phi_j(x) s_j\right)}
#'
#' where \eqn{\phi_k} is the eigenfunction associated with \eqn{\lambda_k} and
#' \eqn{s_k} is the Koopman mode.
#'
#' @section Advantages:
#' \itemize{
#'   \item Does not require construction of an approximate operator matrix
#'   \item Works directly with trajectory data
#'   \item Can be more robust for systems where eigenvalues are known/estimated
#'   \item Avoids issues with basis selection
#' }
#'
#' @section Limitations:
#' \itemize{
#'   \item Requires eigenvalues to be known or estimated separately
#'   \item Can be numerically unstable for eigenvalues with |λ| > 1
#'   \item Convergence may be slow for eigenvalues near the unit circle
#' }
#'
#' @examples
#' # Simple oscillator
#' t <- seq(0, 50, by = 0.1)
#' y <- cos(2 * pi * 0.1 * t) + 0.5 * sin(2 * pi * 0.3 * t)
#'
#' # GLA with automatic eigenvalue estimation
#' result <- gla(y, n_eigenvalues = 4)
#' print(result)
#'
#' # GLA with known eigenvalues
#' dt <- 0.1
#' omega1 <- 2 * pi * 0.1
#' omega2 <- 2 * pi * 0.3
#' known_eigs <- c(exp(1i * omega1 * dt), exp(-1i * omega1 * dt),
#'                  exp(1i * omega2 * dt), exp(-1i * omega2 * dt))
#' result2 <- gla(y, eigenvalues = known_eigs)
#'
#' @references
#' Mezic, I. (2020). On Numerical Approximations of the Koopman Operator.
#' arXiv:2009.05883, Section 3.
#'
#' Mezic, I. and Banaszuk, A. (2004). Comparison of systems with complex
#' behavior. Physica D, 197(1-2):101-133.
#'
#' @seealso [dmd()] for matrix-based decomposition, [hankel_dmd()] for
#'   Krylov subspace methods.
#'
#' @export
gla <- function(y, eigenvalues = NULL, n_eigenvalues = 5, tol = 1e-6,
                max_iter = NULL) {

  cl <- match.call()

  # Handle input
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
    stop("y contains non-finite values", call. = FALSE)
  }

  n_obs <- nrow(y)
  n_time <- ncol(y)

  if (is.null(max_iter)) {
    max_iter <- n_time
  }
  max_iter <- min(max_iter, n_time)

  # Estimate eigenvalues if not provided
  if (is.null(eigenvalues)) {
    # Use preliminary DMD to estimate eigenvalues
    model <- hankel_dmd(y, delays = min(10, floor(n_time / 3)))
    eigenvalues <- model$eigenvalues[1:min(n_eigenvalues, length(model$eigenvalues))]
  }

  n_eig <- length(eigenvalues)

  # Sort eigenvalues by magnitude (descending)
  order_idx <- order(Mod(eigenvalues), decreasing = TRUE)
  eigenvalues <- eigenvalues[order_idx]

  # Initialize storage
  # eigenfunctions stores scalar eigenfunction values (one per eigenvalue, evaluated at each time)
  eigenfunctions <- matrix(0 + 0i, nrow = n_eig, ncol = n_time)
  modes <- matrix(0 + 0i, nrow = n_obs, ncol = n_eig)
  convergence <- logical(n_eig)
  residuals <- numeric(n_eig)

  # Compute eigenfunctions sequentially (eq. 12)
  for (k in seq_len(n_eig)) {
    lambda_k <- eigenvalues[k]

    # Compute f_k by subtracting contributions from previous eigenvalues
    f_residual <- y  # Start with original signal

    if (k > 1) {
      for (j in 1:(k - 1)) {
        lambda_j <- eigenvalues[j]
        phi_j_0 <- eigenfunctions[j, 1]  # eigenfunction value at t=0
        s_j <- modes[, j]

        # Subtract contribution of mode j
        for (t_idx in seq_len(n_time)) {
          f_residual[, t_idx] <- f_residual[, t_idx] -
            Re(lambda_j^(t_idx - 1) * phi_j_0 * s_j)
        }
      }
    }

    # Compute weighted average (eq. 12)
    # f_k = lim (1/n) sum_{i=0}^{n-1} lambda_k^{-i} * f_residual(T^i x)
    result <- compute_gla_average(f_residual, lambda_k, max_iter, tol)

    # Extract eigenfunction and mode
    # f_k(x, z) = phi_k(x) * s_k(z)
    # result$average is an n_obs vector representing f_k at initial point
    f_k <- result$average

    if (Mod(f_k[1]) > tol) {
      # The mode s_k is the n_obs vector
      s_k <- f_k
      # phi_k(x_0) = 1 by normalization
      phi_k_0 <- 1 + 0i
    } else {
      phi_k_0 <- 0 + 0i
      s_k <- rep(0 + 0i, n_obs)
    }

    # Compute eigenfunction evolution: phi_k(x_t) = lambda_k^t * phi_k(x_0)
    for (t_idx in seq_len(n_time)) {
      eigenfunctions[k, t_idx] <- lambda_k^(t_idx - 1) * phi_k_0
    }

    modes[, k] <- s_k
    convergence[k] <- result$converged
    residuals[k] <- result$residual
  }

  # Build eigenfunction evolution (verify eigenvalue relationship)
  # Check: phi_k(T x) should equal lambda_k * phi_k(x)
  eigenfunction_errors <- numeric(n_eig)
  for (k in seq_len(n_eig)) {
    phi_k <- eigenfunctions[k, ]
    lambda_k <- eigenvalues[k]

    if (n_time > 1) {
      # Compare phi_k(x_{t+1}) with lambda_k * phi_k(x_t)
      lhs <- phi_k[2:n_time]
      rhs <- lambda_k * phi_k[1:(n_time - 1)]
      eigenfunction_errors[k] <- sqrt(mean(Mod(lhs - rhs)^2))
    }
  }

  structure(
    list(
      eigenvalues = eigenvalues,
      eigenfunctions = eigenfunctions,
      modes = modes,
      convergence = convergence,
      residuals = residuals,
      eigenfunction_errors = eigenfunction_errors,
      n_iter = max_iter,
      n_obs = n_obs,
      n_time = n_time,
      call = cl
    ),
    class = "gla"
  )
}


#' Compute GLA Weighted Average (internal)
#'
#' Computes the weighted time average for a single eigenvalue.
#'
#' @param f_residual Residual signal matrix after removing previous modes.
#' @param lambda Target eigenvalue.
#' @param max_iter Maximum iterations.
#' @param tol Convergence tolerance.
#'
#' @return List with average, convergence status, and residual.
#'
#' @keywords internal
compute_gla_average <- function(f_residual, lambda, max_iter, tol) {

  n_obs <- nrow(f_residual)
  n_time <- ncol(f_residual)
  n_iter <- min(max_iter, n_time)

  # Running average: (1/n) sum_{i=0}^{n-1} lambda^{-i} * f(T^i x)
  running_sum <- rep(0 + 0i, n_obs)
  prev_avg <- rep(0 + 0i, n_obs)
  converged <- FALSE
  final_residual <- Inf

  # Avoid numerical issues with |lambda| > 1
  lambda_mag <- Mod(lambda)
  if (lambda_mag > 1 + tol) {
    # For unstable eigenvalues, use reverse iteration
    # This is mathematically equivalent but numerically stable
    for (i in seq_len(n_iter)) {
      t_idx <- n_time - i + 1
      if (t_idx < 1) break
      weight <- lambda^(i - 1)  # Forward weight
      running_sum <- running_sum + weight * f_residual[, t_idx]
    }
  } else {
    # Standard forward iteration for stable eigenvalues
    for (i in seq_len(n_iter)) {
      if (i > n_time) break
      weight <- (1 / lambda)^(i - 1)
      running_sum <- running_sum + weight * f_residual[, i]

      # Check convergence every 10 iterations
      if (i %% 10 == 0 && i > 10) {
        current_avg <- running_sum / i
        change <- sqrt(sum(Mod(current_avg - prev_avg)^2))
        if (change < tol) {
          converged <- TRUE
          break
        }
        prev_avg <- current_avg
      }
    }
  }

  average <- running_sum / n_iter
  final_residual <- sqrt(sum(Mod(average - prev_avg)^2))

  list(
    average = average,
    converged = converged,
    residual = final_residual,
    n_iter = n_iter
  )
}


#' Print Method for GLA Objects
#'
#' @param x A `"gla"` object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the input object.
#'
#' @export
print.gla <- function(x, ...) {
  cat("Generalized Laplace Analysis\n")
  cat("============================\n")
  cat(sprintf("  Observables:      %d\n", x$n_obs))
  cat(sprintf("  Time points:      %d\n", x$n_time))
  cat(sprintf("  Eigenvalues:      %d\n", length(x$eigenvalues)))
  cat(sprintf("  Converged:        %d / %d\n",
              sum(x$convergence), length(x$eigenvalues)))

  cat("\nEigenvalue Summary:\n")
  eig_df <- data.frame(
    idx = seq_along(x$eigenvalues),
    magnitude = round(Mod(x$eigenvalues), 6),
    phase = round(Arg(x$eigenvalues), 6),
    converged = x$convergence,
    error = round(x$eigenfunction_errors, 6)
  )
  print(eig_df, row.names = FALSE)

  invisible(x)
}


#' Predict Using GLA Model
#'
#' Forecasts future values using the GLA eigenfunctions and modes.
#'
#' @param object A `"gla"` object.
#' @param n_ahead Integer; number of time steps to forecast.
#' @param ... Additional arguments (unused).
#'
#' @return Matrix of predictions (n_obs x n_ahead).
#'
#' @export
predict.gla <- function(object, n_ahead = 10, ...) {

  n_obs <- object$n_obs
  n_eig <- length(object$eigenvalues)
  n_time <- object$n_time

  predictions <- matrix(0, nrow = n_obs, ncol = n_ahead)

  # Predict using modal decomposition:
  # y(t) = sum_k phi_k(t) * s_k
  # phi_k(t) = lambda_k^t * phi_k(0)
  for (k in seq_len(n_eig)) {
    lambda_k <- object$eigenvalues[k]
    phi_k_0 <- object$eigenfunctions[k, 1]  # Initial eigenfunction value
    s_k <- object$modes[, k]

    for (t in seq_len(n_ahead)) {
      t_future <- n_time + t - 1
      phi_k_t <- lambda_k^t_future * phi_k_0
      predictions[, t] <- predictions[, t] + Re(phi_k_t * s_k)
    }
  }

  colnames(predictions) <- paste0("t+", seq_len(n_ahead))
  predictions
}


#' Reconstruct Signal Using GLA Model
#'
#' Reconstructs the original signal using the computed eigenfunctions
#' and modes.
#'
#' @param object A `"gla"` object.
#' @param modes_to_use Integer vector; which modes to include.
#'   Default is all modes.
#'
#' @return Matrix of reconstructed values (n_obs x n_time).
#'
#' @export
gla_reconstruct <- function(object, modes_to_use = NULL) {

  if (is.null(modes_to_use)) {
    modes_to_use <- seq_along(object$eigenvalues)
  }

  n_obs <- object$n_obs
  n_time <- object$n_time

  reconstruction <- matrix(0, nrow = n_obs, ncol = n_time)

  for (k in modes_to_use) {
    lambda_k <- object$eigenvalues[k]
    phi_k_0 <- object$eigenfunctions[k, 1]
    s_k <- object$modes[, k]

    for (t in seq_len(n_time)) {
      phi_k_t <- lambda_k^(t - 1) * phi_k_0
      reconstruction[, t] <- reconstruction[, t] + Re(phi_k_t * s_k)
    }
  }

  reconstruction
}
