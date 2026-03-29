# Evaluate fitted calibration function at measurement data.

#' Evaluate fitted function at measurement points
#'
#' Given fitted coefficients \code{b} and their covariance \code{b_cov},
#' evaluates the calibration function at each measurement and propagates
#' uncertainty.
#'
#' @param meas_data numeric matrix with \eqn{m} rows and two columns:
#'   instrument response \eqn{y}, standard uncertainty \eqn{u(y)}.
#' @param b numeric vector of fitted coefficients from \code{\link{b_least}}.
#' @param b_cov covariance matrix of \code{b} from \code{\link{b_least}}.
#' @param func calibration function used during fitting.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{x}{numeric vector of assigned values \eqn{x}.}
#'     \item{x_cov}{covariance matrix of \code{x}.}
#'   }
#' @export
b_eval <- function(meas_data, b, b_cov, func) {
  result <- b_eval_xy(meas_data, b, b_cov, func)
  list(x = result$x, x_cov = result$x_cov)
}

#' Evaluate fitted function and return full x-y covariance structure
#'
#' Like \code{\link{b_eval}} but also returns the covariance of the
#' (adjusted) instrument responses and the cross-covariance between
#' assigned values and responses.
#'
#' @inheritParams b_eval
#'
#' @return A list with elements:
#'   \describe{
#'     \item{x}{numeric vector of assigned values.}
#'     \item{y}{numeric vector of instrument responses (same as
#'       \code{meas_data[, 1]}).}
#'     \item{x_cov}{covariance matrix of \code{x} (\eqn{m \times m}).}
#'     \item{y_cov}{covariance matrix of \code{y} (\eqn{m \times m}).}
#'     \item{xy_cov}{cross-covariance matrix of \code{x} and \code{y}
#'       (\eqn{m \times m}).}
#'   }
#' @export
b_eval_xy <- function(meas_data, b, b_cov, func) {
  y   <- meas_data[, 1]
  uy  <- meas_data[, 2]
  ny  <- length(y)
  nb  <- length(b)
  f       <- func(y, b)
  x       <- f$x
  dx_dy   <- f$dx_dy
  dx_db   <- do.call(cbind, f$dx_db)                        # ny x nb
  # diag(scalar) in R creates a scalar-sized identity, not a 1x1 matrix;
  # diag(x, nrow, ncol) avoids this when ny == 1.
  jx      <- cbind(diag(x = dx_dy, nrow = ny, ncol = ny), dx_db)
  jy      <- cbind(diag(ny), matrix(0, ny, nb))
  j       <- rbind(jx, jy)             # 2*ny x (ny + nb)
  cv_in   <- matrix(0, ny + nb, ny + nb)
  cv_in[seq_len(ny), seq_len(ny)]                       <- diag(x = uy^2, nrow = ny, ncol = ny)
  cv_in[(ny + 1):(ny + nb), (ny + 1):(ny + nb)]        <- b_cov
  cov_out <- j %*% cv_in %*% t(j)     # 2*ny x 2*ny
  idx_x <- seq_len(ny)
  idx_y <- (ny + 1):(2 * ny)
  list(
    x      = x,
    y      = y,
    x_cov  = cov_out[idx_x, idx_x, drop = FALSE],
    y_cov  = cov_out[idx_y, idx_y, drop = FALSE],
    xy_cov = cov_out[idx_x, idx_y, drop = FALSE]
  )
}
