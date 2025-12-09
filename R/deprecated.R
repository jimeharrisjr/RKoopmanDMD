#' Legacy DMD Functions (Deprecated)
#'
#' These functions are deprecated and will be removed in a future version.
#' Please use the new API instead.
#'
#' @name deprecated
#' @keywords internal
NULL


#' Create a kdmd Object (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated. Use [dmd()] instead.
#'
#' @param x A matrix.
#'
#' @return An object of class "kdmd".
#'
#' @examples
#' \dontrun{
#' # Old way (deprecated)
#' m <- matrix(1:20, nrow = 2)
#' A <- getAMatrix(m)
#'
#' # New way
#' model <- dmd(m)
#' A <- model$A
#' }
#'
#' @keywords internal
#' @export
kdmd <- function(x) {
  .Deprecated("dmd", package = "RKoopmanDMD",
              msg = "kdmd() is deprecated. Use dmd() instead.")

  if (!is.matrix(x)) stop("X must be matrix")
  structure(as.matrix(x), class = c("kdmd", "matrix", "array"))
}


#' Compute Koopman Matrix (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated. Use [dmd()] instead, which returns
#' a richer object containing the A matrix along with modes, eigenvalues,
#' and other useful components.
#'
#' @param data The matrix data for which you wish to predict future columns
#'   expressed as a 2D matrix or data frame.
#' @param p The percentage of explanation of the SVD required (default = 1).
#' @param comp The number of components of the SVD to keep (> 1).
#'   Default is NULL, overrides p if set.
#'
#' @return An object of class "kdmd" containing the Koopman matrix.
#'
#' @examples
#' \dontrun{
#' # Old way (deprecated)
#' m <- matrix(1:30, nrow = 3)
#' A <- getAMatrix(m)
#'
#' # New way
#' model <- dmd(m)
#' A <- model$A
#' }
#'
#' @seealso [dmd()] for the new API.
#'
#' @keywords internal
#' @export
getAMatrix <- function(data, p = 1, comp = NA) {
  .Deprecated("dmd", package = "RKoopmanDMD",
              msg = paste("getAMatrix() is deprecated.",
                          "Use dmd() instead, which returns a richer object.",
                          "Access the A matrix via model$A."))

  # Convert p to rank for new API
  if (!is.na(comp)) {
    rank <- as.integer(comp)
  } else if (p == 1) {
    rank <- NULL
  } else {
    # The new API uses variance threshold, not percentage
    # This is an approximation
    rank <- NULL
  }

  # Use new dmd function
  model <- dmd(data, rank = rank)

  # Return A matrix wrapped in legacy class for backward compatibility
  A <- model$A
  structure(A, class = c("kdmd", "matrix", "array"))
}


#' Predict from kdmd Object (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This method is deprecated. Use [predict.dmd()] with a dmd object instead.
#'
#' @param object A kdmd object generated with getAMatrix.
#' @param data The matrix data for which you wish to predict future columns.
#' @param l The length of the predictions (number of columns to predict, default = 1).
#' @param ... Additional arguments (unused).
#'
#' @return A matrix with l additional columns appended.
#'
#' @examples
#' \dontrun{
#' # Old way (deprecated)
#' m <- matrix(1:30, nrow = 3)
#' A <- getAMatrix(m)
#' result <- predict(A, m, l = 5)
#'
#' # New way
#' model <- dmd(m)
#' future <- predict(model, n_ahead = 5)
#' result <- cbind(m, future)
#' }
#'
#' @seealso [predict.dmd()] for the new API.
#'
#' @keywords internal
#' @export
predict.kdmd <- function(object, data, l = 1, ...) {
  .Deprecated("predict.dmd", package = "RKoopmanDMD",
              msg = paste("predict.kdmd() is deprecated.",
                          "Use dmd() and predict.dmd() instead."))

  if (!"kdmd" %in% class(object)) {
    stop("object must be a kdmd object generated with getAMatrix",
         call. = FALSE)
  }

  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("data must be a numeric matrix or a data frame", call. = FALSE)
  }

  data <- as.matrix(data)
  A <- object

  len_predict <- as.integer(l)
  t <- ncol(data)
  ynew <- data

  for (st in 0:(len_predict - 1)) {
    b <- Re(A %*% matrix(ynew[, t + st], ncol = 1))
    ynew <- cbind(ynew, b)
  }

  ynew
}


#' Legacy koopman_dmd Function (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated. Use [dmd()] instead.
#'
#' @param X Data matrix where columns are snapshots over time.
#' @param r Rank for truncation (if NULL, uses full rank).
#'
#' @return A list containing modes, eigenvalues, amplitudes, and rank.
#'
#' @keywords internal
#' @export
koopman_dmd <- function(X, r = NULL) {
  .Deprecated("dmd", package = "RKoopmanDMD",
              msg = "koopman_dmd() is deprecated. Use dmd() instead.")

  model <- dmd(X, rank = r)

  # Return in old format for backward compatibility
  list(
    modes = model$modes,
    eigenvalues = model$eigenvalues,
    amplitudes = model$amplitudes,
    rank = model$rank
  )
}


#' Legacy predict_dmd Function (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is deprecated. Use [predict.dmd()] instead.
#'
#' @param dmd_result Output from koopman_dmd().
#' @param n_steps Number of time steps to predict.
#' @param x0 Initial condition.
#'
#' @return Matrix of predictions.
#'
#' @keywords internal
#' @export
predict_dmd <- function(dmd_result, n_steps, x0 = NULL) {
  .Deprecated("predict.dmd", package = "RKoopmanDMD",
              msg = "predict_dmd() is deprecated. Use predict.dmd() instead.")

  modes <- dmd_result$modes
  lambdas <- dmd_result$eigenvalues

  if (is.null(x0)) {
    b <- dmd_result$amplitudes
  } else {
    b <- solve(modes, x0)
  }

  predictions <- matrix(0, nrow = nrow(modes), ncol = n_steps)

  for (k in 1:n_steps) {
    time_evolution <- lambdas^(k - 1)
    predictions[, k] <- Re(modes %*% (b * time_evolution))
  }

  predictions
}
