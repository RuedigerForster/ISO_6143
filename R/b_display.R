# Display / print helper functions.

#' Display calibration data
#'
#' Prints the calibration data matrix to the console.
#'
#' @param cal_data numeric matrix as returned by \code{\link{b_read_cal_data}}.
#' @return Invisibly \code{NULL}.
#' @export
b_disp_cal_data <- function(cal_data) {
  cat("Calibration data:\n")
  print(cal_data)
  invisible(NULL)
}

#' Display calibration results
#'
#' Prints the fitted coefficients, their standard uncertainties, the full
#' covariance matrix, the sum of squared weighted residuals, and the maximum
#' absolute weighted deviation.
#'
#' @param b numeric vector of fitted coefficients.
#' @param b_cov covariance matrix of \code{b}.
#' @param b_res weighted residual vector (or \code{NA} for Monte Carlo output).
#' @return Invisibly \code{NULL}.
#' @export
b_disp_cal_results <- function(b, b_cov, b_res) {
  cat("Coefficients b\n")
  print(b)
  cat("Uncertainties u(b)\n")
  print(sqrt(diag(b_cov)))
  cat("Covariance cov(b)\n")
  print(b_cov)
  cat("Residual\n")
  print(sum(b_res^2))
  cat("Maximum absolute value of weighted deviations\n")
  print(max(abs(b_res)))
  invisible(NULL)
}

#' Display measurement results
#'
#' Prints a table of assigned values, their standard uncertainties, and the
#' original measurement data.  For more than one measurement, also prints the
#' full covariance matrix.
#'
#' @param x numeric vector of assigned values from \code{\link{b_eval}}.
#' @param x_cov covariance matrix of \code{x}.
#' @param meas_data numeric matrix as returned by
#'   \code{\link{b_read_meas_data}}.
#' @return Invisibly \code{NULL}.
#' @export
b_disp_meas_results <- function(x, x_cov, meas_data) {
  cat("Measurement results:\n")
  ux <- sqrt(diag(x_cov))
  print(cbind(x = x, ux = ux, meas_data))
  if (length(ux) > 1) {
    cat("Covariance cov(x)\n")
    print(x_cov)
  }
  invisible(NULL)
}

#' Fit and evaluate in one step with printed output
#'
#' Convenience wrapper: fits \code{func} to \code{cal_data}, evaluates at
#' \code{meas_data}, and prints all intermediate results.
#'
#' @param cal_data numeric matrix as returned by \code{\link{b_read_cal_data}}.
#' @param meas_data numeric matrix as returned by
#'   \code{\link{b_read_meas_data}}.
#' @param func calibration function.
#'
#' @return Invisibly a list with elements \code{b}, \code{b_cov},
#'   \code{b_res}, \code{x}, \code{x_cov}.
#' @export
b_test <- function(cal_data, meas_data, func) {
  b_disp_cal_data(cal_data)
  cal_res  <- b_least(cal_data, func)
  b_disp_cal_results(cal_res$b, cal_res$b_cov, cal_res$b_res)
  meas_res <- b_eval(meas_data, cal_res$b, cal_res$b_cov, func)
  b_disp_meas_results(meas_res$x, meas_res$x_cov, meas_data)
  invisible(c(cal_res, meas_res))
}
