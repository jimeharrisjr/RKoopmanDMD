#' Lifting Functions for Extended DMD
#'
#' Functions for transforming state space to higher-dimensional feature space
#' where nonlinear dynamics become approximately linear.
#'
#' @name lifting
#' @keywords internal
NULL


#' Apply Lifting Transformation to Data
#'
#' Transforms data matrix using a lifting function.
#'
#' @param X Numeric matrix (n_vars x n_time).
#' @param lifting_fn A lifting function that takes a matrix and returns a lifted matrix.
#' @param allow_col_reduction Logical; if TRUE, allow lifting functions that reduce
#'   the number of columns (e.g., time-delay embedding). Default FALSE.
#'
#' @return Lifted data matrix (n_lifted x n_time or n_lifted x reduced_time).
#'
#' @keywords internal
lift_data <- function(X, lifting_fn, allow_col_reduction = FALSE) {
  if (is.null(lifting_fn)) {
    return(X)
  }

  lifted <- lifting_fn(X)

  # Validate output
  if (!is.matrix(lifted)) {
    lifted <- as.matrix(lifted)
  }

  if (!allow_col_reduction && ncol(lifted) != ncol(X)) {
    stop("Lifting function must preserve the number of columns (time points). ",
         "For time-delay embedding, use allow_col_reduction = TRUE.",
         call. = FALSE)
  }

  if (ncol(lifted) > ncol(X)) {
    stop("Lifting function cannot increase the number of columns (time points)",
         call. = FALSE)
  }

  lifted
}


#' Create Lifting Function from Specification
#'
#' Converts a lifting specification (string, function, or list) into a
#' callable lifting function.
#'
#' @param lifting Lifting specification: NULL, character, function, or list.
#' @param X Reference data matrix (used to determine dimensions and validate).
#'
#' @return A lifting function or NULL.
#'
#' @keywords internal
make_lifting_fn <- function(lifting, X = NULL) {

  if (is.null(lifting)) {
    return(NULL)
  }

  # Already a function

if (is.function(lifting)) {
    validate_lifting_fn(lifting, X)
    return(lifting)
  }

  # Character shorthand
  if (is.character(lifting)) {
    lifting <- match.arg(lifting, c("poly2", "poly3", "poly4",
                                     "poly_cross2", "poly_cross3",
                                     "trig", "trig2",
                                     "delay2", "delay3", "delay5"))

    fn <- switch(lifting,
      "poly2" = function(X) lift_polynomial(X, degree = 2),
      "poly3" = function(X) lift_polynomial(X, degree = 3),
      "poly4" = function(X) lift_polynomial(X, degree = 4),
      "poly_cross2" = function(X) lift_poly_cross(X, degree = 2),
      "poly_cross3" = function(X) lift_poly_cross(X, degree = 3),
      "trig" = function(X) lift_trigonometric(X, harmonics = 1),
      "trig2" = function(X) lift_trigonometric(X, harmonics = 2),
      "delay2" = function(X) lift_delay(X, delays = 2),
      "delay3" = function(X) lift_delay(X, delays = 3),
      "delay5" = function(X) lift_delay(X, delays = 5)
    )

    return(fn)
  }

  # List specification with parameters
  if (is.list(lifting)) {
    type <- lifting$type
    if (is.null(type)) {
      stop("List lifting specification must include 'type' element", call. = FALSE)
    }

    fn <- switch(type,
      "poly" = {
        degree <- lifting$degree %||% 2
        function(X) lift_polynomial(X, degree = degree)
      },
      "poly_cross" = {
        degree <- lifting$degree %||% 2
        function(X) lift_poly_cross(X, degree = degree)
      },
      "trig" = {
        harmonics <- lifting$harmonics %||% 1
        function(X) lift_trigonometric(X, harmonics = harmonics)
      },
      "delay" = {
        delays <- lifting$delays %||% 2
        function(X) lift_delay(X, delays = delays)
      },
      "rbf" = {
        centers <- lifting$centers
        sigma <- lifting$sigma %||% 1
        if (is.null(centers)) {
          stop("RBF lifting requires 'centers' matrix", call. = FALSE)
        }
        function(X) lift_rbf(X, centers = centers, sigma = sigma)
      },
      stop(sprintf("Unknown lifting type: '%s'", type), call. = FALSE)
    )

    return(fn)
  }

  stop("lifting must be NULL, a character string, a function, or a list",
       call. = FALSE)
}


#' Validate a Lifting Function
#'
#' Checks that a user-provided lifting function works correctly.
#'
#' @param fn A lifting function.
#' @param X Reference data matrix for testing.
#'
#' @return TRUE invisibly if valid, otherwise throws an error.
#'
#' @keywords internal
validate_lifting_fn <- function(fn, X = NULL) {
  if (!is.function(fn)) {
    stop("Lifting function must be a function", call. = FALSE)
  }

  if (!is.null(X)) {
    # Test on a small subset
    test_X <- X[, 1:min(3, ncol(X)), drop = FALSE]

    result <- tryCatch(
      fn(test_X),
      error = function(e) {
        stop(sprintf("Lifting function failed: %s", e$message), call. = FALSE)
      }
    )

    if (!is.matrix(result) && !is.numeric(result)) {
      stop("Lifting function must return a numeric matrix", call. = FALSE)
    }

    if (ncol(as.matrix(result)) != ncol(test_X)) {
      stop("Lifting function must preserve the number of columns", call. = FALSE)
    }

    if (nrow(as.matrix(result)) < nrow(test_X)) {
      stop("Lifting function should not reduce dimensionality", call. = FALSE)
    }
  }

  invisible(TRUE)
}


# ============================================================================
# Built-in Lifting Functions
# ============================================================================

#' Polynomial Lifting
#'
#' Lifts state space by adding polynomial terms up to specified degree.
#'
#' @param X Numeric matrix (n_vars x n_time).
#' @param degree Maximum polynomial degree (default 2).
#'
#' @return Lifted matrix with dimensions (n_vars * degree) x n_time.
#'   Contains original X and powers X^2, X^3, ..., X^degree.
#'
#' @examples
#' X <- matrix(1:10, nrow = 2)
#' lifted <- RKoopmanDMD:::lift_polynomial(X, degree = 3)
#' # Result has 6 rows: x1, x2, x1^2, x2^2, x1^3, x2^3
#'
#' @keywords internal
lift_polynomial <- function(X, degree = 2) {
  if (degree < 1) {
    stop("degree must be at least 1", call. = FALSE)
  }

  if (degree == 1) {
    return(X)
  }

  lifted <- X
  for (d in 2:degree) {
    lifted <- rbind(lifted, X^d)
  }

  # Add row names for clarity
  n_vars <- nrow(X)
  base_names <- if (!is.null(rownames(X))) rownames(X) else paste0("x", seq_len(n_vars))
  new_names <- base_names
  for (d in 2:degree) {
    new_names <- c(new_names, paste0(base_names, "^", d))
  }
  rownames(lifted) <- new_names

  lifted
}


#' Polynomial Lifting with Cross-terms
#'
#' Lifts state space by adding all polynomial terms including cross-products
#' up to specified total degree.
#'
#' @param X Numeric matrix (n_vars x n_time).
#' @param degree Maximum total polynomial degree (default 2).
#'
#' @return Lifted matrix including all monomials up to total degree.
#'   For n=2, degree=2: x1, x2, x1^2, x1*x2, x2^2
#'
#' @examples
#' X <- matrix(1:10, nrow = 2)
#' lifted <- RKoopmanDMD:::lift_poly_cross(X, degree = 2)
#' # Result: x1, x2, x1^2, x1*x2, x2^2
#'
#' @keywords internal
lift_poly_cross <- function(X, degree = 2) {
  if (degree < 1) {
    stop("degree must be at least 1", call. = FALSE)
  }

  n_vars <- nrow(X)
  n_time <- ncol(X)

  if (degree == 1) {
    return(X)
  }

  # Generate all monomial exponent combinations
  # For efficiency, we'll use a recursive approach
  exponents <- generate_monomials(n_vars, degree)

  # Exclude the constant term (all zeros) - keep only degree >= 1
  exponents <- exponents[rowSums(exponents) >= 1, , drop = FALSE]

  # Compute each monomial
  n_terms <- nrow(exponents)
  lifted <- matrix(0, nrow = n_terms, ncol = n_time)

  for (i in seq_len(n_terms)) {
    exp_vec <- exponents[i, ]
    term <- rep(1, n_time)
    for (j in seq_len(n_vars)) {
      if (exp_vec[j] > 0) {
        term <- term * (X[j, ]^exp_vec[j])
      }
    }
    lifted[i, ] <- term
  }

  # Add row names
  base_names <- if (!is.null(rownames(X))) rownames(X) else paste0("x", seq_len(n_vars))
  row_names <- apply(exponents, 1, function(exp_vec) {
    parts <- character(0)
    for (j in seq_len(n_vars)) {
      if (exp_vec[j] == 1) {
        parts <- c(parts, base_names[j])
      } else if (exp_vec[j] > 1) {
        parts <- c(parts, paste0(base_names[j], "^", exp_vec[j]))
      }
    }
    paste(parts, collapse = "*")
  })
  rownames(lifted) <- row_names

  lifted
}


#' Generate Monomial Exponents
#'
#' Generates all combinations of exponents for monomials up to a given degree.
#'
#' @param n_vars Number of variables.
#' @param max_degree Maximum total degree.
#'
#' @return Matrix where each row is an exponent vector.
#'
#' @keywords internal
generate_monomials <- function(n_vars, max_degree) {
  if (n_vars == 1) {
    return(matrix(0:max_degree, ncol = 1))
  }

  # Recursive generation
  result <- matrix(0, nrow = 0, ncol = n_vars)

  for (d in 0:max_degree) {
    # For each total degree d, generate combinations
    sub_exponents <- generate_monomials(n_vars - 1, d)
    for (i in seq_len(nrow(sub_exponents))) {
      remaining <- d - sum(sub_exponents[i, ])
      if (remaining >= 0) {
        new_row <- c(remaining, sub_exponents[i, ])
        result <- rbind(result, new_row)
      }
    }
  }

  # Remove duplicates and sort
  result <- unique(result)
  result <- result[order(rowSums(result)), , drop = FALSE]

  result
}


#' Trigonometric Lifting
#'
#' Lifts state space by adding sine and cosine terms.
#'
#' @param X Numeric matrix (n_vars x n_time).
#' @param harmonics Number of harmonic frequencies (default 1).
#'
#' @return Lifted matrix with dimensions (n_vars * (1 + 2*harmonics)) x n_time.
#'   Contains: X, sin(X), cos(X), sin(2X), cos(2X), ...
#'
#' @examples
#' X <- matrix(seq(0, 2*pi, length.out = 10), nrow = 1)
#' lifted <- RKoopmanDMD:::lift_trigonometric(X, harmonics = 2)
#'
#' @keywords internal
lift_trigonometric <- function(X, harmonics = 1) {
  if (harmonics < 1) {
    stop("harmonics must be at least 1", call. = FALSE)
  }

  n_vars <- nrow(X)
  lifted <- X

  base_names <- if (!is.null(rownames(X))) rownames(X) else paste0("x", seq_len(n_vars))
  new_names <- base_names

  for (h in seq_len(harmonics)) {
    if (h == 1) {
      lifted <- rbind(lifted, sin(X), cos(X))
      new_names <- c(new_names,
                     paste0("sin(", base_names, ")"),
                     paste0("cos(", base_names, ")"))
    } else {
      lifted <- rbind(lifted, sin(h * X), cos(h * X))
      new_names <- c(new_names,
                     paste0("sin(", h, "*", base_names, ")"),
                     paste0("cos(", h, "*", base_names, ")"))
    }
  }

  rownames(lifted) <- new_names
  lifted
}


#' Time-delay Embedding Lifting
#'
#' Lifts state space by including time-delayed copies of the state.
#' This is useful for scalar time series or when dynamics depend on history.
#'
#' @param X Numeric matrix (n_vars x n_time).
#' @param delays Number of delay steps to include (default 2).
#'
#' @return Lifted matrix with dimensions (n_vars * (delays + 1)) x (n_time - delays).
#'   Note: This reduces the number of usable time points.
#'
#' @details
#' For a scalar time series y(t), delay embedding creates:
#' [y(t), y(t-1), y(t-2), ..., y(t-delays)]
#'
#' This is based on Takens' embedding theorem, which guarantees that
#' for a sufficient embedding dimension, the reconstructed state space
#' is diffeomorphic to the original attractor.
#'
#' @examples
#' y <- sin(seq(0, 4*pi, length.out = 50))
#' Y <- matrix(y, nrow = 1)
#' lifted <- RKoopmanDMD:::lift_delay(Y, delays = 3)
#' # Result has 4 rows and 47 columns
#'
#' @keywords internal
lift_delay <- function(X, delays = 2) {
  if (delays < 1) {
    stop("delays must be at least 1", call. = FALSE)
  }

  n_vars <- nrow(X)
  n_time <- ncol(X)

  if (n_time <= delays) {
    stop(sprintf("Not enough time points (%d) for %d delays", n_time, delays),
         call. = FALSE)
  }

  # Number of valid time points after embedding
  n_valid <- n_time - delays

  # Build the delay-embedded matrix
  lifted <- matrix(0, nrow = n_vars * (delays + 1), ncol = n_valid)

  for (d in 0:delays) {
    row_start <- d * n_vars + 1
    row_end <- (d + 1) * n_vars
    col_start <- delays - d + 1
    col_end <- n_time - d
    lifted[row_start:row_end, ] <- X[, col_start:col_end, drop = FALSE]
  }

  # Add row names
  base_names <- if (!is.null(rownames(X))) rownames(X) else paste0("x", seq_len(n_vars))
  new_names <- base_names
  for (d in seq_len(delays)) {
    new_names <- c(new_names, paste0(base_names, "(t-", d, ")"))
  }
  rownames(lifted) <- new_names

  lifted
}


#' Radial Basis Function Lifting
#'
#' Lifts state space by adding Gaussian radial basis function features.
#'
#' @param X Numeric matrix (n_vars x n_time).
#' @param centers Matrix of RBF centers (n_vars x n_centers).
#' @param sigma Width parameter for Gaussian RBFs (default 1).
#'
#' @return Lifted matrix with dimensions (n_vars + n_centers) x n_time.
#'   Contains original X plus RBF activations for each center.
#'
#' @details
#' Each RBF feature is computed as:
#' phi_i(x) = exp(-||x - c_i||^2 / (2 * sigma^2))
#'
#' where c_i is the i-th center.
#'
#' @examples
#' X <- matrix(rnorm(20), nrow = 2, ncol = 10)
#' centers <- matrix(rnorm(6), nrow = 2, ncol = 3)
#' lifted <- RKoopmanDMD:::lift_rbf(X, centers, sigma = 0.5)
#'
#' @keywords internal
lift_rbf <- function(X, centers, sigma = 1) {
  if (!is.matrix(centers)) {
    centers <- as.matrix(centers)
  }

  if (nrow(centers) != nrow(X)) {
    stop(sprintf("centers must have %d rows (same as X)", nrow(X)), call. = FALSE)
  }

  if (sigma <= 0) {
    stop("sigma must be positive", call. = FALSE)
  }

  n_vars <- nrow(X)
  n_time <- ncol(X)
  n_centers <- ncol(centers)

  # Compute RBF activations
  rbf_features <- matrix(0, nrow = n_centers, ncol = n_time)

  for (i in seq_len(n_centers)) {
    center_i <- centers[, i]
    # Compute squared distance from each column of X to center_i
    diff <- X - center_i  # broadcasts center_i across columns
    sq_dist <- colSums(diff^2)
    rbf_features[i, ] <- exp(-sq_dist / (2 * sigma^2))
  }

  # Combine original data with RBF features
  lifted <- rbind(X, rbf_features)

  # Add row names
  base_names <- if (!is.null(rownames(X))) rownames(X) else paste0("x", seq_len(n_vars))
  rbf_names <- paste0("rbf_", seq_len(n_centers))
  rownames(lifted) <- c(base_names, rbf_names)

  lifted
}


# ============================================================================
# User-facing Functions
# ============================================================================

#' Inspect Lifting Transformation
#'
#' Examines the lifting transformation used in a DMD model, or applies
#' the lifting function to a new state vector.
#'
#' @param object A `"dmd"` object created by [dmd()].
#' @param x0 Optional state vector to transform. If provided, returns the
#'   lifted version of this state.
#'
#' @return If `x0` is NULL, returns a list with lifting information:
#' \describe{
#'   \item{type}{The lifting specification used}
#'   \item{original_dim}{Number of original state variables}
#'   \item{lifted_dim}{Number of lifted state variables}
#'   \item{observables}{Indices of observable (original) states in lifted space}
#' }
#' If `x0` is provided, returns the lifted state vector.
#'
#' @examples
#' # Create data
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#'
#' # Fit DMD with polynomial lifting
#' model <- dmd(X, lifting = "poly2")
#'
#' # Inspect lifting info
#' dmd_lift(model)
#'
#' # Lift a specific state
#' x0 <- c(1, 0)
#' dmd_lift(model, x0)
#'
#' @seealso [dmd()] for fitting models with lifting functions.
#'
#' @export
dmd_lift <- function(object, x0 = NULL) {
  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  if (is.null(object$lifting)) {
    message("No lifting function was used in this model")
    return(invisible(NULL))
  }

  if (!is.null(x0)) {
    # Apply lifting to provided state
    x0 <- as.numeric(x0)
    if (length(x0) != object$n_vars_original) {
      stop(sprintf("x0 must have length %d (original state dimension)",
                   object$n_vars_original), call. = FALSE)
    }
    x0_matrix <- matrix(x0, ncol = 1)
    lifted <- lift_data(x0_matrix, object$lifting_fn)
    return(as.vector(lifted))
  }

  # Return lifting info
  list(
    type = object$lifting,
    original_dim = object$n_vars_original,
    lifted_dim = object$n_vars_lifted,
    observables = object$observables,
    expansion_ratio = object$n_vars_lifted / object$n_vars_original
  )
}


#' List Available Lifting Functions
#'
#' Returns information about built-in lifting function options.
#'
#' @return A data frame describing available lifting functions.
#'
#' @examples
#' list_lifting_functions()
#'
#' @export
list_lifting_functions <- function() {
  data.frame(
    name = c("poly2", "poly3", "poly4",
             "poly_cross2", "poly_cross3",
             "trig", "trig2",
             "delay2", "delay3", "delay5"),
    description = c(
      "Polynomial degree 2 (x, x^2)",
      "Polynomial degree 3 (x, x^2, x^3)",
      "Polynomial degree 4 (x, x^2, x^3, x^4)",
      "Polynomial with cross-terms degree 2",
      "Polynomial with cross-terms degree 3",
      "Trigonometric (x, sin(x), cos(x))",
      "Trigonometric 2 harmonics",
      "Time-delay embedding (2 delays)",
      "Time-delay embedding (3 delays)",
      "Time-delay embedding (5 delays)"
    ),
    parametric = c(
      "list(type='poly', degree=N)",
      "list(type='poly', degree=N)",
      "list(type='poly', degree=N)",
      "list(type='poly_cross', degree=N)",
      "list(type='poly_cross', degree=N)",
      "list(type='trig', harmonics=N)",
      "list(type='trig', harmonics=N)",
      "list(type='delay', delays=N)",
      "list(type='delay', delays=N)",
      "list(type='delay', delays=N)"
    ),
    stringsAsFactors = FALSE
  )
}


#' Check if Lifting Specification is Delay-based
#'
#' Determines if a lifting specification involves time-delay embedding,
#' which reduces the number of usable time points.
#'
#' @param lifting Lifting specification
#'
#' @return Logical; TRUE if delay-based
#'
#' @keywords internal
is_delay_lifting <- function(lifting) {
  if (is.null(lifting)) {
    return(FALSE)
  }

  if (is.character(lifting)) {
    return(grepl("^delay", lifting))
  }

  if (is.list(lifting) && !is.null(lifting$type)) {
    return(lifting$type == "delay")
  }

  # For custom functions, we can't know - assume not delay
  FALSE
}


# Note: %||% operator is defined in utils.R
