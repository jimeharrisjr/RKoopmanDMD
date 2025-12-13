#' RKoopmanDMD: Koopman Operator and Dynamic Mode Decomposition
#'
#' @description
#' The RKoopmanDMD package provides tools for analyzing dynamical systems using
#' Dynamic Mode Decomposition (DMD) and the Koopman operator framework. DMD is
#' a data-driven method that extracts spatiotemporal coherent structures from
#' time-series data without requiring knowledge of the underlying equations.
#'
#' @section Main Functions:
#' \describe{
#'   \item{[dmd()]}{Perform Dynamic Mode Decomposition on time-series data}
#'   \item{[hankel_dmd()]}{Krylov subspace DMD using time-delay embedding (avoids curse of dimensionality)}
#'   \item{[gla()]}{Generalized Laplace Analysis - compute eigenfunctions via weighted time averages}
#'   \item{[predict.dmd()]}{Forecast future states using a fitted DMD model}
#'   \item{[dmd_spectrum()]}{Analyze the eigenvalue spectrum for stability}
#'   \item{[dmd_reconstruct()]}{Reconstruct training data from DMD modes}
#'   \item{[dmd_residual()]}{Compute residual for error assessment}
#'   \item{[dmd_pseudospectrum()]}{Pseudospectral analysis for eigenvalue reliability}
#'   \item{[dmd_convergence()]}{Estimate convergence rate of DMD approximation}
#' }
#'
#' @section The Koopman Operator:
#' The Koopman operator is an infinite-dimensional linear operator that
#' describes the evolution of observables (functions of state) for a nonlinear
#' dynamical system. DMD provides a finite-dimensional approximation to this
#' operator by computing a best-fit linear map between time-shifted snapshots
#' of the data.
#'
#' @section Numerical Methods:
#' This package implements several numerical approaches based on Mezic (2020):
#' \describe{
#'   \item{Standard DMD}{Finite section method using SVD-based projection}
#'   \item{Hankel-DMD}{Krylov subspace method using time-delayed observables.
#'     Avoids the curse of dimensionality by letting dynamics select the basis.}
#'   \item{GLA}{Generalized Laplace Analysis for direct eigenfunction computation
#'     without constructing an operator matrix.}
#' }
#'
#' @section Typical Workflow:
#' ```
#' # 1. Organize data as matrix (rows = variables, columns = time)
#' X <- matrix(...)
#'
#' # 2. Fit DMD model (choose method based on data)
#' model <- dmd(X)           # Standard DMD
#' model <- hankel_dmd(y)    # For scalar time series
#' model <- dmd(X, lifting = "poly2")  # Extended DMD with lifting
#'
#' # 3. Examine eigenvalues for stability
#' summary(model)
#' plot(model)
#'
#' # 4. Assess approximation quality
#' res <- dmd_residual(model, X)
#' ps <- dmd_pseudospectrum(model)
#'
#' # 5. Predict future states
#' future <- predict(model, n_ahead = 50)
#' ```
#'
#' @references
#' Mezic, I. (2020). On Numerical Approximations of the Koopman Operator.
#' arXiv:2009.05883.
#'
#' Schmid, P. J. (2010). Dynamic mode decomposition of numerical and
#' experimental data. Journal of Fluid Mechanics, 656, 5-28.
#'
#' Kutz, J. N., Brunton, S. L., Brunton, B. W., & Proctor, J. L. (2016).
#' Dynamic Mode Decomposition: Data-Driven Modeling of Complex Systems. SIAM.
#'
#' Williams, M. O., Kevrekidis, I. G., & Rowley, C. W. (2015). A data-driven
#' approximation of the Koopman operator: Extending dynamic mode decomposition.
#' Journal of Nonlinear Science, 25(6), 1307-1346.
#'
#' @keywords internal
"_PACKAGE"
