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
#'   \item{[predict.dmd()]}{Forecast future states using a fitted DMD model}
#'   \item{[dmd_spectrum()]}{Analyze the eigenvalue spectrum for stability}
#'   \item{[dmd_reconstruct()]}{Reconstruct training data from DMD modes}
#' }
#'
#' @section The Koopman Operator:
#' The Koopman operator is an infinite-dimensional linear operator that
#' describes the evolution of observables (functions of state) for a nonlinear
#' dynamical system. DMD provides a finite-dimensional approximation to this
#' operator by computing a best-fit linear map between time-shifted snapshots
#' of the data.
#'
#' @section Typical Workflow:
#' ```
#' # 1. Organize data as matrix (rows = variables, columns = time)
#' X <- matrix(...)
#'
#' # 2. Fit DMD model
#' model <- dmd(X)
#'
#' # 3. Examine eigenvalues for stability
#' summary(model)
#' plot(model)
#'
#' # 4. Predict future states
#' future <- predict(model, n_ahead = 50)
#' ```
#'
#' @references
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
