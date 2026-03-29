# S3 methods for class "B_Least".

#' Print a B_Least calibration result
#'
#' Displays the fitted coefficients, their standard uncertainties, the full
#' covariance matrix, the sum of squared weighted residuals, and the maximum
#' absolute weighted deviation.
#'
#' @param x a \code{B_Least} object returned by \code{\link{b_least}}.
#' @param ... ignored.
#' @return Invisibly \code{x}.
#' @export
print.B_Least <- function(x, ...) {
  b_disp_cal_results(x$b, x$b_cov, x$b_res)
  invisible(x)
}

#' Plot a B_Least calibration result
#'
#' Draws the fitted calibration curve with uncertainty band, calibration
#' reference points, and measurement points with their propagated
#' uncertainties.  The calibration data are taken from the stored
#' \code{B_Least} object; only the measurement data need to be supplied.
#'
#' @param x a \code{B_Least} object returned by \code{\link{b_least}}.
#' @param meas_data numeric matrix as returned by
#'   \code{\link{b_read_meas_data}}.
#' @param ... ignored.
#' @return Invisibly \code{NULL}.
#' @export
plot.B_Least <- function(x, meas_data, ...) {
  b_plot(x$cal_data, meas_data, x$b, x$b_cov, x$func)
}
