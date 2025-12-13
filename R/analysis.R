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


#' Compute DMD Residual for Error Assessment
#'
#' Computes the residual of the DMD approximation, which measures how well
#' the finite-dimensional approximation captures the true Koopman operator.
#' Based on the error analysis in Mezic (2020), Section 4.3 and 5.2.
#'
#' @param object A `"dmd"` or `"hankel_dmd"` object.
#' @param X_original Optional; the original data matrix. If not provided,
#'   uses stored data from the model (if available).
#'
#' @return A list containing:
#' \describe{
#'   \item{residual_norm}{Overall Frobenius norm of residual}
#'   \item{residual_relative}{Residual norm relative to data norm}
#'   \item{per_step_residual}{Vector of residual norms at each time step}
#'   \item{per_mode_residual}{Vector of residual contributions per mode}
#'   \item{pseudospectral_bound}{Upper bound on pseudospectral radius (epsilon)}
#'   \item{residual_matrix}{Full residual matrix (X_pred - X_true)}
#' }
#'
#' @details
#' The residual measures the error in the finite section approximation:
#' \deqn{r = U\tilde{\phi} - \tilde{\lambda}\tilde{\phi}}
#'
#' From Proposition 4.1 in Mezic (2020):
#' \deqn{U\tilde{\phi} - \tilde{\lambda}\tilde{\phi} = \tilde{e} \cdot (U\tilde{f} - P_{\tilde{F}}U\tilde{f})}
#'
#' The pseudospectral bound \eqn{\epsilon} satisfies:
#' \deqn{|U\tilde{e} \cdot \tilde{f} - \tilde{\lambda}\tilde{e} \cdot \tilde{f}| = |e_N r|}
#'
#' This means the computed eigenvalues are in the \eqn{\epsilon}-pseudospectrum
#' of the true Koopman operator.
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#' model <- dmd(X)
#' res <- dmd_residual(model, X)
#' cat("Relative residual:", res$residual_relative, "\n")
#' cat("Pseudospectral bound:", res$pseudospectral_bound, "\n")
#'
#' @references
#' Mezic, I. (2020). On Numerical Approximations of the Koopman Operator.
#' arXiv:2009.05883, Sections 4.3 and 5.2.
#'
#' @seealso [dmd_error()] for reconstruction error metrics,
#'   [dmd_pseudospectrum()] for pseudospectral analysis.
#'
#' @export
dmd_residual <- function(object, X_original = NULL) {

  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  # Get data dimensions
  n_vars <- object$data_dim[1]
  n_time <- object$data_dim[2]

  # Reconstruct predictions using the DMD matrix
  # For each time step k, predict: x_{k+1} = A x_k
  if (inherits(object, "hankel_dmd")) {
    H <- object$hankel
    X_data <- H
  } else if (!is.null(X_original)) {
    X_data <- validate_matrix(X_original, min_rows = 1, min_cols = 2)
    # Apply lifting if needed
    if (!is.null(object$lifting_fn)) {
      X_data <- lift_data(X_data, object$lifting_fn,
                          allow_col_reduction = is_delay_lifting(object$lifting))
    }
    # Center if needed
    if (object$center) {
      X_data <- X_data - object$X_mean
    }
  } else {
    stop("X_original must be provided for standard DMD objects", call. = FALSE)
  }

  n_cols <- ncol(X_data)

  # Compute one-step predictions
  X1 <- X_data[, -n_cols, drop = FALSE]
  X2_true <- X_data[, -1, drop = FALSE]
  X2_pred <- object$A %*% X1

  # Residual matrix
  residual_matrix <- X2_pred - X2_true

  # Overall residual norm
  residual_norm <- sqrt(sum(residual_matrix^2))
  data_norm <- sqrt(sum(X2_true^2))
  residual_relative <- residual_norm / data_norm

  # Per-step residual
  per_step_residual <- sqrt(colSums(residual_matrix^2))

  # Per-mode residual contribution
  # Project residual onto each mode
  modes <- object$modes
  per_mode_residual <- numeric(object$rank)

  for (i in seq_len(object$rank)) {
    mode_i <- modes[, i]
    # Project residual onto this mode direction
    mode_contrib <- abs(t(Conj(mode_i)) %*% residual_matrix)
    per_mode_residual[i] <- sqrt(sum(Mod(mode_contrib)^2)) / Mod(sqrt(sum(mode_i * Conj(mode_i))))
  }

  # Pseudospectral bound (eq. 102)
  # epsilon = |e_N * r| where e_N is last component of eigenvector
  # For each eigenvalue, compute the bound
  eig_vecs <- safe_solve(modes, diag(n_vars))  # columns are eigenvector coefficients
  pseudospectral_bounds <- numeric(object$rank)

  for (i in seq_len(object$rank)) {
    # Get the eigenvector in coefficient space
    e_coef <- safe_solve(modes, modes[, i])
    e_N <- abs(e_coef[length(e_coef)])
    # Average residual contribution
    r_avg <- mean(per_step_residual)
    pseudospectral_bounds[i] <- e_N * r_avg
  }

  pseudospectral_bound <- max(pseudospectral_bounds)

  list(
    residual_norm = residual_norm,
    residual_relative = residual_relative,
    per_step_residual = per_step_residual,
    per_mode_residual = per_mode_residual,
    pseudospectral_bound = pseudospectral_bound,
    pseudospectral_bounds = pseudospectral_bounds,
    residual_matrix = residual_matrix
  )
}


#' Pseudospectrum Analysis for DMD
#'
#' Computes and optionally visualizes the pseudospectrum of the DMD operator.
#' The pseudospectrum reveals regions where eigenvalue estimates are reliable.
#'
#' @param object A `"dmd"` object.
#' @param epsilon Numeric vector of epsilon values for pseudospectrum contours.
#'   Default is `c(0.01, 0.05, 0.1, 0.2)`.
#' @param grid_n Integer; number of grid points in each direction for
#'   computing the pseudospectrum. Default is 100.
#' @param xlim,ylim Limits for the complex plane grid. If NULL, automatically
#'   determined from eigenvalues.
#' @param plot Logical; if TRUE, produce a plot. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{z_grid}{Complex grid points}
#'   \item{sigma_min}{Matrix of minimum singular values of (zI - A)}
#'   \item{eigenvalues}{DMD eigenvalues for reference}
#'   \item{epsilon}{Epsilon values used}
#' }
#'
#' @details
#' The \eqn{\epsilon}-pseudospectrum of an operator A is defined as:
#' \deqn{\sigma_\epsilon(A) = \{z \in \mathbb{C} : \|(zI - A)^{-1}\| \geq 1/\epsilon\}}
#'
#' Equivalently, it is the set of z where the minimum singular value of
#' (zI - A) is at most \eqn{\epsilon}.
#'
#' From Theorem 5.1 in Mezic (2020), the Krylov subspace approximation
#' converges in the pseudospectral sense: for any \eqn{\epsilon > 0}, for
#' large enough n, the approximate eigenfunctions are in the
#' \eqn{\epsilon}-pseudospectrum.
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#' model <- dmd(X)
#'
#' # Compute pseudospectrum (with plot)
#' ps <- dmd_pseudospectrum(model)
#'
#' # Compute without plotting
#' ps <- dmd_pseudospectrum(model, plot = FALSE)
#'
#' @references
#' Mezic, I. (2020). On Numerical Approximations of the Koopman Operator.
#' arXiv:2009.05883, Section 5.2.
#'
#' Trefethen, L.N. and Embree, M. (2005). Spectra and Pseudospectra.
#' Princeton University Press.
#'
#' @export
dmd_pseudospectrum <- function(object, epsilon = c(0.01, 0.05, 0.1, 0.2),
                                grid_n = 100, xlim = NULL, ylim = NULL,
                                plot = TRUE) {

  if (!inherits(object, "dmd")) {
    stop("object must be a dmd object", call. = FALSE)
  }

  A <- object$A_tilde  # Use reduced matrix for efficiency
  lambdas <- object$eigenvalues
  n <- nrow(A)

  # Determine grid limits
  if (is.null(xlim)) {
    re_range <- range(Re(lambdas))
    re_margin <- max(0.5, diff(re_range) * 0.3)
    xlim <- c(re_range[1] - re_margin, re_range[2] + re_margin)
  }

  if (is.null(ylim)) {
    im_range <- range(Im(lambdas))
    im_margin <- max(0.5, diff(im_range) * 0.3)
    ylim <- c(im_range[1] - im_margin, im_range[2] + im_margin)
  }

  # Create grid
  x_seq <- seq(xlim[1], xlim[2], length.out = grid_n)
  y_seq <- seq(ylim[1], ylim[2], length.out = grid_n)

  # Compute minimum singular value at each grid point
  sigma_min <- matrix(0, nrow = grid_n, ncol = grid_n)

  I_n <- diag(n)

  for (i in seq_len(grid_n)) {
    for (j in seq_len(grid_n)) {
      z <- complex(real = x_seq[i], imaginary = y_seq[j])
      M <- z * I_n - A
      sv <- svd(M, nu = 0, nv = 0)$d
      sigma_min[i, j] <- min(sv)
    }
  }

  result <- list(
    x = x_seq,
    y = y_seq,
    sigma_min = sigma_min,
    eigenvalues = lambdas,
    epsilon = epsilon
  )

  if (plot) {
    plot_pseudospectrum(result)
  }

  invisible(result)
}


#' Plot Pseudospectrum (internal)
#' @keywords internal
plot_pseudospectrum <- function(ps_result) {

  x <- ps_result$x
  y <- ps_result$y
  sigma_min <- ps_result$sigma_min
  lambdas <- ps_result$eigenvalues
  epsilon <- ps_result$epsilon

  # Draw unit circle
  theta <- seq(0, 2 * pi, length.out = 100)
  circle_x <- cos(theta)
  circle_y <- sin(theta)

  # Set up plot
  graphics::plot(range(x), range(y), type = "n", asp = 1,
                 xlab = "Real", ylab = "Imaginary",
                 main = "DMD Pseudospectrum")

  # Add contours for each epsilon
  colors <- grDevices::heat.colors(length(epsilon) + 2)
  colors <- colors[seq_len(length(epsilon))]

  for (k in seq_along(epsilon)) {
    graphics::contour(x, y, sigma_min, levels = epsilon[k],
                      add = TRUE, col = colors[k], lwd = 1.5,
                      drawlabels = TRUE)
  }

  # Add unit circle
  graphics::lines(circle_x, circle_y, lty = 2, col = "gray50")

  # Add axes
  graphics::abline(h = 0, v = 0, col = "gray80", lty = 3)

  # Add eigenvalues
  graphics::points(Re(lambdas), Im(lambdas), pch = 19, col = "blue", cex = 1.2)

  # Legend
  graphics::legend("topright",
                   legend = c(paste0("eps=", epsilon), "Eigenvalues", "Unit circle"),
                   col = c(colors, "blue", "gray50"),
                   lty = c(rep(1, length(epsilon)), NA, 2),
                   pch = c(rep(NA, length(epsilon)), 19, NA),
                   bty = "n", cex = 0.8)
}


#' Estimate Convergence Rate
#'
#' Estimates the convergence rate of the DMD approximation by computing
#' eigenvalues at increasing sample sizes and measuring convergence.
#'
#' @param X Data matrix (n_vars x n_time).
#' @param sample_fractions Numeric vector of fractions of data to use.
#'   Default is `c(0.25, 0.5, 0.75, 1.0)`.
#' @param rank Rank for DMD (NULL for auto).
#' @param ... Additional arguments passed to [dmd()].
#'
#' @return A list containing:
#' \describe{
#'   \item{sample_sizes}{Vector of sample sizes used}
#'   \item{eigenvalues}{List of eigenvalue vectors at each sample size}
#'   \item{eigenvalue_changes}{Vector of max eigenvalue changes between successive fits}
#'   \item{convergence_estimate}{Estimated convergence rate}
#' }
#'
#' @details
#' From Theorem 4.3 in Mezic (2020), for systems with pure point spectrum,
#' the finite section approximation converges at rate O(1/m) where m is
#' the number of samples. This function empirically estimates this rate.
#'
#' @examples
#' t <- seq(0, 20, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#' conv <- dmd_convergence(X)
#' cat("Estimated convergence rate:", conv$convergence_estimate, "\n")
#'
#' @export
dmd_convergence <- function(X, sample_fractions = c(0.25, 0.5, 0.75, 1.0),
                             rank = NULL, ...) {

  X <- validate_matrix(X, min_rows = 1, min_cols = 10)
  n_time <- ncol(X)

  sample_sizes <- floor(sample_fractions * n_time)
  sample_sizes <- pmax(sample_sizes, 5)  # minimum 5 samples
  sample_sizes <- unique(sample_sizes)

  eigenvalues_list <- list()
  models <- list()

  for (i in seq_along(sample_sizes)) {
    m <- sample_sizes[i]
    X_sub <- X[, 1:m, drop = FALSE]

    model <- dmd(X_sub, rank = rank, ...)
    models[[i]] <- model
    eigenvalues_list[[i]] <- model$eigenvalues
  }

  # Compute eigenvalue changes
  eigenvalue_changes <- numeric(length(sample_sizes) - 1)

  for (i in 2:length(sample_sizes)) {
    # Match eigenvalues by nearest neighbor
    prev_eig <- eigenvalues_list[[i - 1]]
    curr_eig <- eigenvalues_list[[i]]

    # Use minimum of two lengths
    n_compare <- min(length(prev_eig), length(curr_eig))

    if (n_compare > 0) {
      # Sort by magnitude for comparison
      prev_sorted <- sort(Mod(prev_eig), decreasing = TRUE)[1:n_compare]
      curr_sorted <- sort(Mod(curr_eig), decreasing = TRUE)[1:n_compare]
      eigenvalue_changes[i - 1] <- max(abs(prev_sorted - curr_sorted))
    }
  }

  # Estimate convergence rate using log-log regression
  # If rate is O(1/m^alpha), then log(change) ~ -alpha * log(m)
  if (length(sample_sizes) >= 3 && any(eigenvalue_changes > 0)) {
    valid_idx <- eigenvalue_changes > 0
    if (sum(valid_idx) >= 2) {
      log_m <- log(sample_sizes[-1][valid_idx])
      log_change <- log(eigenvalue_changes[valid_idx])
      fit <- stats::lm(log_change ~ log_m)
      convergence_rate <- -stats::coef(fit)[2]
    } else {
      convergence_rate <- NA
    }
  } else {
    convergence_rate <- NA
  }

  list(
    sample_sizes = sample_sizes,
    eigenvalues = eigenvalues_list,
    eigenvalue_changes = eigenvalue_changes,
    convergence_estimate = convergence_rate,
    models = models
  )
}
