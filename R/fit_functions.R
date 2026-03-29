# Fit functions and objective function for B LEAST (ISO 6143:2001)
# Each fit function returns a list with three elements:
#   x      - vector of fitted x values
#   dx_dy  - vector of partial derivatives dx/dy
#   dx_db  - list of vectors, one per coefficient: partial derivatives dx/db_j

#' Linear calibration function
#'
#' Evaluates \eqn{x = b_1 + b_2 y} and its partial derivatives.
#'
#' @param y numeric vector of instrument responses.
#' @param b numeric vector of two coefficients \eqn{(b_1, b_2)}.
#' @return A list with elements \code{x}, \code{dx_dy}, and \code{dx_db}.
#' @export
b_linear_func <- function(y, b) {
  x     <- b[1] + b[2] * y
  dx_dy <- b[2] + 0 * y
  dx_db <- list(1 + 0 * y, y)
  list(x = x, dx_dy = dx_dy, dx_db = dx_db)
}

#' Second-order polynomial calibration function
#'
#' Evaluates \eqn{x = b_1 + b_2 y + b_3 y^2} and its partial derivatives.
#'
#' @param y numeric vector of instrument responses.
#' @param b numeric vector of three coefficients \eqn{(b_1, b_2, b_3)}.
#' @return A list with elements \code{x}, \code{dx_dy}, and \code{dx_db}.
#' @export
b_second_order_poly <- function(y, b) {
  x     <- b[1] + b[2] * y + b[3] * y^2
  dx_dy <- b[2] + 2 * b[3] * y
  dx_db <- list(1 + 0 * y, y, y^2)
  list(x = x, dx_dy = dx_dy, dx_db = dx_db)
}

#' Third-order polynomial calibration function
#'
#' Evaluates \eqn{x = b_1 + b_2 y + b_3 y^2 + b_4 y^3} and its partial
#' derivatives.
#'
#' @param y numeric vector of instrument responses.
#' @param b numeric vector of four coefficients \eqn{(b_1, b_2, b_3, b_4)}.
#' @return A list with elements \code{x}, \code{dx_dy}, and \code{dx_db}.
#' @export
b_third_order_poly <- function(y, b) {
  x     <- b[1] + b[2] * y + b[3] * y^2 + b[4] * y^3
  dx_dy <- b[2] + 2 * b[3] * y + 3 * b[4] * y^2
  dx_db <- list(1 + 0 * y, y, y^2, y^3)
  list(x = x, dx_dy = dx_dy, dx_db = dx_db)
}

#' Power calibration function
#'
#' Evaluates \eqn{x = b_1 + b_2 y^{(1 + b_3)}} and its partial derivatives.
#'
#' @param y numeric vector of instrument responses (must be positive).
#' @param b numeric vector of three coefficients \eqn{(b_1, b_2, b_3)}.
#' @return A list with elements \code{x}, \code{dx_dy}, and \code{dx_db}.
#' @export
b_power_func <- function(y, b) {
  x     <- b[1] + b[2] * y^(1 + b[3])
  dx_dy <- b[2] * (b[3] + 1) * y^b[3]
  dx_db <- list(1 + 0 * y,
                y^(1 + b[3]),
                b[2] * y^(1 + b[3]) * log(y))
  list(x = x, dx_dy = dx_dy, dx_db = dx_db)
}

#' Exponential calibration function
#'
#' Evaluates \eqn{x = b_1 + b_2 e^{b_3 y}} and its partial derivatives.
#'
#' @param y numeric vector of instrument responses.
#' @param b numeric vector of three coefficients \eqn{(b_1, b_2, b_3)}.
#' @return A list with elements \code{x}, \code{dx_dy}, and \code{dx_db}.
#' @export
b_exp_func <- function(y, b) {
  e     <- exp(b[3] * y)
  x     <- b[1] + b[2] * e
  dx_dy <- b[2] * b[3] * e
  dx_db <- list(1 + 0 * y, e, b[2] * y * e)
  list(x = x, dx_dy = dx_dy, dx_db = dx_db)
}

# Weighted residuals for calibration data (errors-in-both-variables).
# Returns a vector of length 2n: weighted x deviations followed by
# weighted y deviations.
.b_objective_func <- function(x, ux, y, uy, y2, b, func) {
  f   <- func(y2, b)
  x2  <- f$x
  wdx <- (x2 - x) / ux
  wdy <- (y2 - y) / uy
  c(wdx, wdy)
}
