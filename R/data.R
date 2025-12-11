#' International Data Analytics Olympiad 2020 (IDAO 2020) Competition Satellite Data
#'
#' A subset of data from the IDAO 2020 Competition containing satellite positions
#'
#'
#' @format ## `sat0`
#' A data frame with 741 rows and 15 columns:
#' \describe{
#'   \item{id}{Measurement ID}
#'   \item{epoch}{timestamp of measurement}
#'   \item{sat_id}{ID of satellite}
#'   \item{x,y,z}{real x, y, and z coordinate of the satellite}
#'   \item{Vx,Vy,Vz}{x-y-z velocity components}
#'   \item{x_sim,y_sim,z_sim}{simulated x, y, and z coordinate of the satellite}
#'   \item{Vx_sim,Vy_sim,Vz_sim}{simulated x-y-z velocity components}
#'   ...
#' }
#' @source <https://www.kaggle.com/datasets/idawoodjee/predict-the-positions-and-speeds-of-600-satellites>
"sat0"
