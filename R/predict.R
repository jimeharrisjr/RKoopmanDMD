#' Predict Future States Using DMD
#'
#' Forecasts future states of a dynamical system using a fitted DMD model.
#' Supports two prediction methods: direct matrix multiplication or mode-based
#' eigenvalue evolution.
#'
#' @param object A `"dmd"` object created by [dmd()].
#' @param n_ahead Integer; number of time steps to forecast. Default is 10.
#' @param x0 Numeric vector; initial state for prediction. If `NULL` (default),
#'   uses the last observed state from the training data. Can be any state
#'   vector of the same dimension as the original state variables.
#' @param method Character; prediction method to use:
#'   \describe{
#'     \item{"modes"}{(Default) Uses eigenvalue evolution of modes. More
#'       accurate for long-term predictions and provides insight into
#'       modal contributions.
#'     \item{"matrix"}{Direct iterated matrix multiplication. Faster but
#'       may accumulate numerical errors over many steps.}
#'   }
#' @param ... Additional arguments (currently unused, for S3 compatibility).
#'
#' @return A numeric matrix with dimensions (n_vars x n_ahead), where each
#'   column represents the predicted state at successive future time steps.
#'
#' @details
#' ## Mode-based Prediction (method = "modes")
#' Each DMD mode evolves independently according to its eigenvalue:
#' \deqn{x(t) = \sum_i b_i \lambda_i^t \phi_i}
#' where \eqn{b_i} are the mode amplitudes, \eqn{\lambda_i} are eigenvalues,
#' and \eqn{\phi_i} are the modes.
#'
#' ## Matrix-based Prediction (method = "matrix")
#' Iteratively applies the DMD matrix:
#' \deqn{x(t+1) = A x(t)}
#'
#' The mode-based method is generally preferred because it:
#' - Provides insight into which modes dominate at different times
#' - Is more numerically stable for long predictions
#' - Allows analytical computation without iteration
#'
#' @examples
#' # Create oscillatory data
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#'
#' # Fit DMD model
#' model <- dmd(X)
#'
#' # Predict next 50 time steps from last observed state
#' future <- predict(model, n_ahead = 50)
#'
#' # Predict from a custom initial condition
#' custom_start <- c(0.5, 0.5)
#' future_custom <- predict(model, n_ahead = 50, x0 = custom_start)
#'
#' # Compare prediction methods
#' pred_modes <- predict(model, n_ahead = 20, method = "modes")
#' pred_matrix <- predict(model, n_ahead = 20, method = "matrix")
#'
#' @seealso [dmd()] for fitting the model, [dmd_reconstruct()] for
#'   reconstructing training data.
#'
#' @export
predict.dmd <- function(object, n_ahead = 10, x0 = NULL,
                        method = c("modes", "matrix"), ...) {

  # Validate inputs
  method <- match.arg(method)
  n_ahead <- as.integer(n_ahead)

  if (n_ahead < 1) {
    stop("n_ahead must be a positive integer", call. = FALSE)
  }

  n_vars <- object$data_dim[1]

  # Determine initial condition
  if (is.null(x0)) {
    x0 <- object$X_last
    # If data was centered, we need the uncentered version
    if (object$center) {
      x0 <- x0 - object$X_mean
    }
  } else {
    x0 <- as.numeric(x0)
    if (length(x0) != n_vars) {
      stop(sprintf("x0 must have length %d (number of state variables)",
                   n_vars), call. = FALSE)
    }
    # If centered, need to adjust x0
    if (object$center) {
      x0 <- x0 - object$X_mean
    }
  }

  # Perform prediction
  if (method == "matrix") {
    predictions <- predict_matrix(object, x0, n_ahead)
  } else {
    predictions <- predict_modes(object, x0, n_ahead)
  }

  # Add back the mean if data was centered
  if (object$center) {
    predictions <- predictions + object$X_mean
  }

  # Add column names for time steps
  colnames(predictions) <- paste0("t+", seq_len(n_ahead))

  predictions
}


#' Matrix-based prediction (internal)
#' @keywords internal
predict_matrix <- function(object, x0, n_ahead) {
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


#' Mode-based prediction (internal)
#' @keywords internal
predict_modes <- function(object, x0, n_ahead) {
  n_vars <- length(x0)
  predictions <- matrix(0, nrow = n_vars, ncol = n_ahead)

  modes <- object$modes
  lambdas <- object$eigenvalues

  # Project initial condition onto modes
  b <- safe_solve(modes, x0)

  for (k in seq_len(n_ahead)) {
    # Each mode evolves as lambda^k
    evolution <- lambdas^k
    predictions[, k] <- Re(modes %*% (b * evolution))
  }

  predictions
}


#' Forecast with Confidence Bounds (Experimental)
#'
#' Generates predictions with uncertainty estimates based on reconstruction
#' error from the training data.
#'
#' @param object A `"dmd"` object.
#' @param n_ahead Number of time steps to forecast.
#' @param x0 Initial state (optional).
#' @param level Confidence level (default 0.95 for 95% bounds).
#'
#' @return A list containing:
#' \describe{
#'   \item{fit}{Point predictions (matrix)}
#'   \item{lower}{Lower confidence bound (matrix)}
#'   \item{upper}{Upper confidence bound (matrix)}
#'   \item{level}{Confidence level used}
#' }
#'
#' @details
#' This is an experimental feature. The confidence bounds are based on the
#' reconstruction error from the training data and assume that prediction
#' error grows proportionally with forecast horizon. More sophisticated
#' uncertainty quantification may be added in future versions.
#'
#' @examples
#' \dontrun{
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t)) + rnorm(2 * length(t), sd = 0.1)
#' model <- dmd(X)
#' forecast <- dmd_forecast(model, n_ahead = 20)
#' }
#'
#' @export
dmd_forecast <- function(object, n_ahead = 10, x0 = NULL, level = 0.95) {

  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  # Get point predictions
  fit <- predict(object, n_ahead = n_ahead, x0 = x0, method = "modes")

  # Estimate prediction uncertainty from reconstruction error
  recon <- dmd_reconstruct(object)
  error <- object$X_first + if (object$center) object$X_mean else 0

  # This is a simplified uncertainty estimate
  # In practice, more sophisticated methods could be used
  n_time <- object$data_dim[2]
  X_orig <- matrix(0, nrow = object$data_dim[1], ncol = n_time)

  # Reconstruct original data for error calculation
  for (k in seq_len(n_time)) {
    evolution <- object$eigenvalues^(k - 1)
    X_orig[, k] <- Re(object$modes %*% (object$amplitudes * evolution))
  }
  if (object$center) {
    X_orig <- X_orig + object$X_mean
  }

  # Calculate per-variable RMSE (simplified)
  rmse <- sqrt(mean((recon - X_orig[, seq_len(ncol(recon))])^2))

  # Scale uncertainty with forecast horizon (simple linear growth model)
  z <- stats::qnorm((1 + level) / 2)
  horizon_scale <- sqrt(seq_len(n_ahead))

  lower <- fit - z * rmse * matrix(horizon_scale, nrow = nrow(fit),
                                    ncol = n_ahead, byrow = TRUE)
  upper <- fit + z * rmse * matrix(horizon_scale, nrow = nrow(fit),
                                    ncol = n_ahead, byrow = TRUE)

  list(
    fit = fit,
    lower = lower,
    upper = upper,
    level = level
  )
}
