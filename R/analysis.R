#' DMD Eigenvalue Spectrum Analysis
#'
#' Analyzes the eigenvalue spectrum of a DMD model to characterize
#' system dynamics including stability, oscillation frequencies, and
#' growth/decay rates.
#'
#' @param object A `"dmd"` object.
#' @param dt Numeric; time step between snapshots. Default is 1.
#'   Used to convert eigenvalue phases to physical frequencies.
#'
#' @return A data frame with one row per mode containing:
#' \describe{
#'   \item{mode}{Mode index}
#'   \item{eigenvalue}{Complex eigenvalue}
#'   \item{magnitude}{Eigenvalue magnitude (|lambda|)}
#'   \item{phase}{Eigenvalue phase in radians}
#'   \item{frequency}{Oscillation frequency (cycles per time unit)}
#'   \item{period}{Oscillation period (time units per cycle)}
#'   \item{growth_rate}{Continuous-time growth rate}
#'   \item{half_life}{Time for mode to decay to half (if decaying)}
#'   \item{doubling_time}{Time for mode to double (if growing)}
#'   \item{stability}{Character: "decaying", "neutral", or "growing"}
#' }
#'
#' @details
#' For discrete-time systems with time step dt:
#' - Frequency: \eqn{f = \arg(\lambda) / (2\pi \cdot dt)}
#' - Growth rate: \eqn{\sigma = \log(|\lambda|) / dt}
#' - Half-life (for decaying modes): \eqn{t_{1/2} = -\log(2) / \sigma}
#'
#' @examples
#' # Create damped oscillator data
#' t <- seq(0, 10, by = 0.1)
#' x <- exp(-0.1 * t) * cos(2 * pi * 0.5 * t)
#' y <- exp(-0.1 * t) * sin(2 * pi * 0.5 * t)
#' X <- rbind(x, y)
#'
#' model <- dmd(X)
#' spectrum <- dmd_spectrum(model, dt = 0.1)
#' print(spectrum)
#'
#' @seealso [dmd()], [summary.dmd()]
#'
#' @export
dmd_spectrum <- function(object, dt = 1) {

  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  lambdas <- object$eigenvalues
  n_modes <- length(lambdas)

  magnitudes <- Mod(lambdas)
  phases <- Arg(lambdas)

  # Convert to physical quantities
  frequencies <- phases / (2 * pi * dt)
  periods <- ifelse(abs(frequencies) > 1e-10, 1 / abs(frequencies), Inf)

  growth_rates <- log(magnitudes) / dt

  # Half-life for decaying modes, doubling time for growing
  half_life <- ifelse(growth_rates < 0, -log(2) / growth_rates, NA)
  doubling_time <- ifelse(growth_rates > 0, log(2) / growth_rates, NA)

  # Stability classification
  stability <- ifelse(magnitudes < 0.9999, "decaying",
                      ifelse(magnitudes > 1.0001, "growing", "neutral"))

  data.frame(
    mode = seq_len(n_modes),
    eigenvalue = lambdas,
    magnitude = magnitudes,
    phase = phases,
    frequency = frequencies,
    period = periods,
    growth_rate = growth_rates,
    half_life = half_life,
    doubling_time = doubling_time,
    stability = stability,
    stringsAsFactors = FALSE
  )
}


#' Reconstruct Training Data from DMD Model
#'
#' Reconstructs the original time series using the DMD modes and
#' eigenvalues. Useful for assessing reconstruction accuracy and
#' understanding mode contributions.
#'
#' @param object A `"dmd"` object.
#' @param n_steps Integer; number of time steps to reconstruct.
#'   Default is the original number of time steps.
#' @param modes Integer vector; which modes to include in reconstruction.
#'   Default is all modes. Use this to see contributions of specific modes.
#'
#' @return A matrix with the same dimensions as the training data
#'   (or specified n_steps), representing the reconstructed time series.
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#'
#' model <- dmd(X)
#'
#' # Full reconstruction
#' X_recon <- dmd_reconstruct(model)
#'
#' # Reconstruction error
#' error <- sqrt(mean((X - X_recon)^2))
#' cat("RMSE:", error, "\n")
#'
#' # Reconstruct using only first mode
#' X_mode1 <- dmd_reconstruct(model, modes = 1)
#'
#' @seealso [dmd()], [dmd_error()]
#'
#' @export
dmd_reconstruct <- function(object, n_steps = NULL, modes = NULL) {

  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  if (is.null(n_steps)) {
    n_steps <- object$data_dim[2]
  }

  n_vars <- object$data_dim[1]
  all_modes <- object$modes
  all_lambdas <- object$eigenvalues
  all_b <- object$amplitudes

  # Select modes
  if (is.null(modes)) {
    modes_idx <- seq_len(object$rank)
  } else {
    modes_idx <- as.integer(modes)
    if (any(modes_idx < 1 | modes_idx > object$rank)) {
      stop(sprintf("modes must be integers between 1 and %d", object$rank),
           call. = FALSE)
    }
  }

  selected_modes <- all_modes[, modes_idx, drop = FALSE]
  selected_lambdas <- all_lambdas[modes_idx]
  selected_b <- all_b[modes_idx]

  # Reconstruct
  reconstruction <- matrix(0, nrow = n_vars, ncol = n_steps)

  for (k in seq_len(n_steps)) {
    evolution <- selected_lambdas^(k - 1)
    reconstruction[, k] <- Re(selected_modes %*% (selected_b * evolution))
  }

  # Add back mean if centered
  if (object$center) {
    reconstruction <- reconstruction + object$X_mean
  }

  colnames(reconstruction) <- paste0("t", seq_len(n_steps))

  reconstruction
}


#' Compute DMD Reconstruction Error
#'
#' Calculates various error metrics comparing the DMD reconstruction
#' to the original data.
#'
#' @param object A `"dmd"` object.
#' @param X_original Optional; the original data matrix. If not provided,
#'   uses the first column and predicts forward.
#'
#' @return A list containing:
#' \describe{
#'   \item{rmse}{Root mean squared error}
#'   \item{mae}{Mean absolute error}
#'   \item{mape}{Mean absolute percentage error (if data contains no zeros)}
#'   \item{relative_error}{Frobenius norm of error / Frobenius norm of data}
#'   \item{max_error}{Maximum absolute error}
#'   \item{per_variable}{Named vector of RMSE per state variable}
#' }
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#'
#' model <- dmd(X)
#' errors <- dmd_error(model, X)
#' print(errors)
#'
#' @seealso [dmd_reconstruct()]
#'
#' @export
dmd_error <- function(object, X_original = NULL) {

  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  if (is.null(X_original)) {
    stop("X_original must be provided to compute error", call. = FALSE)
  }

  X_original <- validate_matrix(X_original, min_rows = 1, min_cols = 2)

  # Reconstruct
  n_steps <- ncol(X_original)
  X_recon <- dmd_reconstruct(object, n_steps = n_steps)

  # Compute errors
  diff <- X_original - X_recon

  rmse <- sqrt(mean(diff^2))
  mae <- mean(abs(diff))

  # Relative error (Frobenius norm)
  rel_error <- sqrt(sum(diff^2)) / sqrt(sum(X_original^2))

  # MAPE (avoid division by zero)
  if (any(X_original == 0)) {
    mape <- NA
  } else {
    mape <- mean(abs(diff / X_original)) * 100
  }

  # Per-variable RMSE
  per_var_rmse <- sqrt(rowMeans(diff^2))
  names(per_var_rmse) <- paste0("var", seq_len(nrow(X_original)))

  list(
    rmse = rmse,
    mae = mae,
    mape = mape,
    relative_error = rel_error,
    max_error = max(abs(diff)),
    per_variable = per_var_rmse
  )
}


#' Extract Dominant DMD Modes
#'
#' Identifies and extracts the most important DMD modes based on
#' amplitude or other criteria.
#'
#' @param object A `"dmd"` object.
#' @param n Integer; number of dominant modes to return. Default is 3.
#' @param criterion Character; how to rank modes:
#'   \describe{
#'     \item{"amplitude"}{(Default) By magnitude of initial amplitudes}
#'     \item{"energy"}{By mode energy (amplitude * magnitude)}
#'     \item{"stability"}{By eigenvalue magnitude (least stable first)}
#'   }
#'
#' @return A list containing:
#' \describe{
#'   \item{indices}{Integer vector of mode indices (ranked)}
#'   \item{eigenvalues}{Complex eigenvalues of dominant modes}
#'   \item{modes}{Matrix of dominant modes (columns)}
#'   \item{amplitudes}{Complex amplitudes of dominant modes}
#'   \item{criterion_values}{Values used for ranking}
#' }
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t), 0.1 * cos(3 * t))
#'
#' model <- dmd(X)
#' dominant <- dmd_dominant_modes(model, n = 2)
#' print(dominant$indices)
#'
#' @export
dmd_dominant_modes <- function(object, n = 3,
                               criterion = c("amplitude", "energy", "stability")) {

  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  criterion <- match.arg(criterion)
  n <- min(n, object$rank)

  # Compute ranking values
  amps <- Mod(object$amplitudes)
  mags <- Mod(object$eigenvalues)

  ranking_values <- switch(criterion,
                           amplitude = amps,
                           energy = amps * mags,
                           stability = mags)

  # Get top n indices
  if (criterion == "stability") {
    # For stability, we want largest magnitudes (most unstable/persistent)
    indices <- order(ranking_values, decreasing = TRUE)[seq_len(n)]
  } else {
    indices <- order(ranking_values, decreasing = TRUE)[seq_len(n)]
  }

  list(
    indices = indices,
    eigenvalues = object$eigenvalues[indices],
    modes = object$modes[, indices, drop = FALSE],
    amplitudes = object$amplitudes[indices],
    criterion_values = ranking_values[indices]
  )
}


#' Check DMD Model Stability
#'
#' Performs a detailed stability analysis of the DMD model.
#'
#' @param object A `"dmd"` object.
#' @param tol Numeric; tolerance for classifying eigenvalues on the
#'   unit circle as "neutral". Default is 1e-6.
#'
#' @return A list with stability information:
#' \describe{
#'   \item{is_stable}{Logical; TRUE if all eigenvalues inside unit circle}
#'   \item{is_unstable}{Logical; TRUE if any eigenvalue outside unit circle}
#'   \item{is_marginal}{Logical; TRUE if any eigenvalue on unit circle}
#'   \item{spectral_radius}{Maximum eigenvalue magnitude}
#'   \item{n_stable}{Number of decaying modes}
#'   \item{n_unstable}{Number of growing modes}
#'   \item{n_neutral}{Number of neutral modes}
#'   \item{classification}{Character vector classifying each mode}
#' }
#'
#' @examples
#' # Stable system
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(exp(-0.1 * t) * cos(t), exp(-0.1 * t) * sin(t))
#' model <- dmd(X)
#' dmd_stability(model)
#'
#' @export
dmd_stability <- function(object, tol = 1e-6) {

  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  mags <- Mod(object$eigenvalues)
  spectral_radius <- max(mags)

  classification <- ifelse(mags < 1 - tol, "stable",
                           ifelse(mags > 1 + tol, "unstable", "neutral"))

  list(
    is_stable = all(mags < 1 - tol),
    is_unstable = any(mags > 1 + tol),
    is_marginal = any(abs(mags - 1) <= tol),
    spectral_radius = spectral_radius,
    n_stable = sum(classification == "stable"),
    n_unstable = sum(classification == "unstable"),
    n_neutral = sum(classification == "neutral"),
    classification = classification
  )
}
