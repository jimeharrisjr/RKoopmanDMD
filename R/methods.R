#' Print Method for DMD Objects
#'
#' Displays a concise summary of a DMD model.
#'
#' @param x A `"dmd"` object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#' model <- dmd(X)
#' print(model)
#'
#' @export
print.dmd <- function(x, ...) {
  cat("Dynamic Mode Decomposition Model\n")
  cat("================================\n")

  # Show lifting info if present
  if (!is.null(x$lifting)) {
    cat(sprintf("  Original vars:    %d\n", x$n_vars_original))
    cat(sprintf("  Lifted vars:      %d\n", x$n_vars_lifted))
    lift_desc <- if (is.character(x$lifting)) {
      x$lifting
    } else if (is.function(x$lifting)) {
      "custom function"
    } else if (is.list(x$lifting)) {
      paste0(x$lifting$type, " (parametric)")
    } else {
      "unknown"
    }
    cat(sprintf("  Lifting:          %s\n", lift_desc))
  } else {
    cat(sprintf("  State variables:  %d\n", x$data_dim[1]))
  }

  cat(sprintf("  Time snapshots:   %d\n", x$data_dim[2]))
  cat(sprintf("  Rank used:        %d\n", x$rank))
  cat(sprintf("  Data centered:    %s\n", if (x$center) "yes" else "no"))

  # Quick stability check (use same tolerance as dmd_stability)
  mags <- Mod(x$eigenvalues)
  tol <- 1e-6
  if (all(mags < 1 - tol)) {
    cat("  Stability:        stable (all |eigenvalues| < 1)\n")
  } else if (any(mags > 1 + tol)) {
    n_unstable <- sum(mags > 1 + tol)
    cat(sprintf("  Stability:        unstable (%d mode(s) with |eigenvalue| > 1)\n",
                n_unstable))
  } else {
    cat("  Stability:        marginally stable (eigenvalues near unit circle)\n")
  }

  cat("\nUse summary() for detailed eigenvalue analysis.\n")
  cat("Use plot() to visualize the eigenvalue spectrum.\n")
  if (!is.null(x$lifting)) {
    cat("Use dmd_lift() to inspect lifting transformation.\n")
  }

  invisible(x)
}


#' Summary Method for DMD Objects
#'
#' Provides detailed summary statistics for a DMD model, including
#' eigenvalue analysis and mode characteristics.
#'
#' @param object A `"dmd"` object.
#' @param ... Additional arguments (unused).
#'
#' @return A list of class `"summary.dmd"` containing:
#' \describe{
#'   \item{n_vars}{Number of state variables}
#'   \item{n_time}{Number of time snapshots}
#'   \item{rank}{Rank used for decomposition}
#'   \item{eigenvalues}{Complex eigenvalues}
#'   \item{eigenvalue_table}{Data frame with eigenvalue properties}
#'   \item{stability}{Stability classification}
#'   \item{variance_explained}{Proportion of variance explained by retained modes}
#' }
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#' model <- dmd(X)
#' summary(model)
#'
#' @export
summary.dmd <- function(object, ...) {

  lambdas <- object$eigenvalues
  n_modes <- length(lambdas)

  # Compute eigenvalue properties
  magnitudes <- Mod(lambdas)
  phases <- Arg(lambdas)  # in radians
  frequencies <- phases / (2 * pi * object$dt)  # cycles per unit time
  growth_rates <- log(magnitudes) / object$dt

  # Stability classification
  stability <- ifelse(magnitudes < 0.9999, "decaying",
                      ifelse(magnitudes > 1.0001, "growing", "neutral"))

  # Mode amplitudes (magnitude)
  amp_mags <- Mod(object$amplitudes)

  # Create summary table
  eig_table <- data.frame(
    mode = seq_len(n_modes),
    magnitude = round(magnitudes, 6),
    phase = round(phases, 6),
    frequency = round(frequencies, 6),
    growth_rate = round(growth_rates, 6),
    stability = stability,
    amplitude = round(amp_mags, 6)
  )

  # Overall stability
  overall_stability <- check_stability(lambdas)

  # Variance explained by retained singular values
  total_var <- sum(object$svd$S^2)
  # We'd need the full SVD to compute fraction, so just report singular values
  sv_info <- data.frame(
    component = seq_along(object$svd$S),
    singular_value = object$svd$S,
    variance_prop = object$svd$S^2 / total_var
  )

  result <- structure(
    list(
      n_vars = object$data_dim[1],
      n_time = object$data_dim[2],
      rank = object$rank,
      centered = object$center,
      eigenvalues = lambdas,
      eigenvalue_table = eig_table,
      stability = overall_stability,
      singular_values = sv_info,
      call = object$call,
      # Lifting info
      lifting = object$lifting,
      n_vars_original = object$n_vars_original,
      n_vars_lifted = object$n_vars_lifted
    ),
    class = "summary.dmd"
  )

  result
}


#' Print Method for DMD Summary
#'
#' @param x A `"summary.dmd"` object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the input object.
#' @export
print.summary.dmd <- function(x, ...) {
  cat("Dynamic Mode Decomposition - Summary\n")
  cat("====================================\n\n")

  cat("Data:\n")
  if (!is.null(x$lifting)) {
    cat(sprintf("  Original vars:    %d\n", x$n_vars_original))
    cat(sprintf("  Lifted vars:      %d\n", x$n_vars_lifted))
    lift_desc <- if (is.character(x$lifting)) {
      x$lifting
    } else if (is.function(x$lifting)) {
      "custom function"
    } else if (is.list(x$lifting)) {
      paste0(x$lifting$type, " (parametric)")
    } else {
      "unknown"
    }
    cat(sprintf("  Lifting:          %s\n", lift_desc))
  } else {
    cat(sprintf("  State variables:  %d\n", x$n_vars))
  }
  cat(sprintf("  Time snapshots:   %d\n", x$n_time))
  cat(sprintf("  Rank used:        %d\n", x$rank))
  cat(sprintf("  Data centered:    %s\n\n", if (x$centered) "yes" else "no"))

  cat("Stability Analysis:\n")
  if (x$stability$stable) {
    cat("  System is STABLE (all eigenvalues inside unit circle)\n")
  } else if (x$stability$unstable) {
    cat(sprintf("  System is UNSTABLE (%d mode(s) growing)\n",
                x$stability$n_unstable))
  } else {
    cat("  System is MARGINALLY STABLE (eigenvalues on unit circle)\n")
  }
  cat(sprintf("  Maximum eigenvalue magnitude: %.6f\n\n", x$stability$max_magnitude))

  cat("Eigenvalue Analysis:\n")
  print(x$eigenvalue_table, row.names = FALSE)

  cat("\nSingular Values (retained):\n")
  print(x$singular_values, row.names = FALSE)

  invisible(x)
}


#' Plot Method for DMD Objects
#'
#' Creates diagnostic plots for a DMD model, including the eigenvalue
#' spectrum and mode amplitudes.
#'
#' @param x A `"dmd"` object.
#' @param which Integer vector specifying which plots to produce:
#'   \describe{
#'     \item{1}{Eigenvalue spectrum in complex plane with unit circle}
#'     \item{2}{Eigenvalue magnitudes (bar plot)}
#'     \item{3}{Mode amplitudes (bar plot)}
#'     \item{4}{Singular value decay}
#'   }
#'   Default is `c(1, 2)`.
#' @param ... Additional graphical parameters passed to plotting functions.
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#' model <- dmd(X)
#'
#' # Default plots
#' plot(model)
#'
#' # All diagnostic plots
#' par(mfrow = c(2, 2))
#' plot(model, which = 1:4)
#'
#' @export
plot.dmd <- function(x, which = c(1, 2), ...) {

  lambdas <- x$eigenvalues
  mags <- Mod(lambdas)
  amps <- Mod(x$amplitudes)

  for (w in which) {
    if (w == 1) {
      # Eigenvalue spectrum in complex plane
      plot_eigenvalue_spectrum(lambdas, ...)
    } else if (w == 2) {
      # Eigenvalue magnitudes
      plot_eigenvalue_magnitudes(mags, ...)
    } else if (w == 3) {
      # Mode amplitudes
      plot_mode_amplitudes(amps, ...)
    } else if (w == 4) {
      # Singular value decay
      plot_singular_values(x$svd$S, ...)
    }
  }

  invisible(x)
}


#' Plot eigenvalue spectrum (internal)
#' @keywords internal
plot_eigenvalue_spectrum <- function(eigenvalues, ...) {
  # Draw unit circle
  theta <- seq(0, 2 * pi, length.out = 100)
  circle_x <- cos(theta)
  circle_y <- sin(theta)

  # Set up plot
  re <- Re(eigenvalues)
  im <- Im(eigenvalues)

  xlim <- range(c(re, circle_x)) * 1.1
  ylim <- range(c(im, circle_y)) * 1.1

  graphics::plot(circle_x, circle_y, type = "l", lty = 2, col = "gray50",
                 xlim = xlim, ylim = ylim, asp = 1,
                 xlab = "Real", ylab = "Imaginary",
                 main = "DMD Eigenvalue Spectrum")
  graphics::abline(h = 0, v = 0, col = "gray80", lty = 3)
  graphics::points(re, im, pch = 19, col = "steelblue", cex = 1.5)

  # Color unstable eigenvalues
  unstable <- Mod(eigenvalues) > 1
  if (any(unstable)) {
    graphics::points(re[unstable], im[unstable], pch = 19, col = "red", cex = 1.5)
  }

  graphics::legend("topright",
                   legend = c("Eigenvalues", "Unstable", "Unit circle"),
                   pch = c(19, 19, NA),
                   lty = c(NA, NA, 2),
                   col = c("steelblue", "red", "gray50"),
                   bty = "n")
}


#' Plot eigenvalue magnitudes (internal)
#' @keywords internal
plot_eigenvalue_magnitudes <- function(magnitudes, ...) {
  n <- length(magnitudes)
  cols <- ifelse(magnitudes > 1, "red",
                 ifelse(magnitudes > 0.99, "orange", "steelblue"))

  graphics::barplot(magnitudes,
                    names.arg = seq_len(n),
                    col = cols,
                    xlab = "Mode",
                    ylab = "Magnitude",
                    main = "Eigenvalue Magnitudes")
  graphics::abline(h = 1, lty = 2, col = "gray50")
}


#' Plot mode amplitudes (internal)
#' @keywords internal
plot_mode_amplitudes <- function(amplitudes, ...) {
  n <- length(amplitudes)
  graphics::barplot(amplitudes,
                    names.arg = seq_len(n),
                    col = "steelblue",
                    xlab = "Mode",
                    ylab = "Amplitude",
                    main = "Mode Amplitudes")
}


#' Plot singular value decay (internal)
#' @keywords internal
plot_singular_values <- function(singular_values, ...) {
  n <- length(singular_values)
  sv_norm <- singular_values / max(singular_values)

  graphics::plot(seq_len(n), sv_norm, type = "b",
                 pch = 19, col = "steelblue",
                 xlab = "Component",
                 ylab = "Normalized Singular Value",
                 main = "Singular Value Decay",
                 log = "y")
  graphics::grid()
}


#' Extract Coefficients from DMD Model
#'
#' Extracts eigenvalues or mode amplitudes from a DMD model.
#'
#' @param object A `"dmd"` object.
#' @param type Character; what to extract:
#'   \describe{
#'     \item{"eigenvalues"}{(Default) Complex eigenvalues}
#'     \item{"amplitudes"}{Complex mode amplitudes}
#'     \item{"modes"}{DMD modes matrix}
#'   }
#' @param ... Additional arguments (unused).
#'
#' @return Depending on `type`: complex vector or complex matrix.
#'
#' @examples
#' t <- seq(0, 10, by = 0.1)
#' X <- rbind(cos(t), sin(t))
#' model <- dmd(X)
#'
#' # Get eigenvalues
#' coef(model)
#'
#' # Get amplitudes
#' coef(model, type = "amplitudes")
#'
#' @export
coef.dmd <- function(object, type = c("eigenvalues", "amplitudes", "modes"), ...) {
  type <- match.arg(type)

  switch(type,
         eigenvalues = object$eigenvalues,
         amplitudes = object$amplitudes,
         modes = object$modes)
}
