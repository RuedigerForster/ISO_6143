# Monte Carlo supplement for B LEAST.

#' Sample calibration data for Monte Carlo
#'
#' Draws \code{nsamples} independent samples from the calibration data by
#' treating each \eqn{(x_i, y_i)} pair as normally distributed with standard
#' deviations \eqn{u(x_i)} and \eqn{u(y_i)}.
#'
#' @param cal_data numeric matrix as returned by \code{\link{b_read_cal_data}}.
#' @param nsamples integer number of Monte Carlo samples (default 10000).
#' @param seed optional integer seed for reproducibility.
#'
#' @return A three-dimensional array of shape
#'   \code{[nsamples, n, 2]}, where \code{[i, j, 1]} is the \eqn{j}-th
#'   sampled \eqn{x} value in the \eqn{i}-th trial, and
#'   \code{[i, j, 2]} is the corresponding sampled \eqn{y} value.
#' @export
b_sample_cal_data_mc <- function(cal_data, nsamples = 10000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n           <- nrow(cal_data)
  cal_samples <- array(0, dim = c(nsamples, n, 2))
  for (j in seq_len(n)) {
    cal_samples[, j, 1] <- rnorm(nsamples, cal_data[j, 1], cal_data[j, 2])
    cal_samples[, j, 2] <- rnorm(nsamples, cal_data[j, 3], cal_data[j, 4])
  }
  cal_samples
}

#' Sample measurement data for Monte Carlo
#'
#' Draws \code{nsamples} independent samples from the measurement data by
#' treating each \eqn{y_i} as normally distributed with standard deviation
#' \eqn{u(y_i)}.
#'
#' @param meas_data numeric matrix as returned by
#'   \code{\link{b_read_meas_data}}.
#' @param nsamples integer number of Monte Carlo samples (default 10000).
#' @param seed optional integer seed for reproducibility.
#'
#' @return A matrix of shape \code{[nsamples, m]}.
#' @export
b_sample_meas_data_mc <- function(meas_data, nsamples = 10000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n            <- nrow(meas_data)
  meas_samples <- matrix(0, nrow = nsamples, ncol = n)
  for (j in seq_len(n)) {
    meas_samples[, j] <- rnorm(nsamples, meas_data[j, 1], meas_data[j, 2])
  }
  meas_samples
}

#' Fit calibration function to each Monte Carlo sample
#'
#' For each of the \code{nsamples} calibration trials, runs
#' Levenberg-Marquardt (without analytic Jacobian for speed) and collects the
#' fitted coefficients.
#'
#' @param cal_samples array as returned by
#'   \code{\link{b_sample_cal_data_mc}}.
#' @param func calibration function.
#'
#' @return A matrix of shape \code{[nsamples, nb]} where \code{nb} is the
#'   number of coefficients.
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @export
b_least_mc <- function(cal_samples, func) {
  nsamples   <- dim(cal_samples)[1]
  cal_data   <- .b_cal_samples_to_cal_data(cal_samples)
  n          <- nrow(cal_data)
  b_start    <- .b_least_start(cal_data, func)
  y2_b_start <- c(cal_data[, 3], b_start)
  y2_b_scale <- y2_b_start
  y2_b_scale[y2_b_scale == 0] <- 1
  y2_b_start2 <- y2_b_start / y2_b_scale
  cal_data_i  <- cal_data
  b_samples   <- matrix(0, nrow = nsamples, ncol = length(b_start))
  for (i in seq_len(nsamples)) {
    cal_data_i[, 1] <- cal_samples[i, , 1]
    cal_data_i[, 3] <- cal_samples[i, , 2]
    lm_i <- nls.lm(
      par     = y2_b_start2,
      fn      = .b_residuals,
      control = nls.lm.control(maxiter = 1000),
      cal_data   = cal_data_i,
      y2_b_scale = y2_b_scale,
      func       = func
    )
    b_samples[i, ] <- (lm_i$par * y2_b_scale)[-seq_len(n)]
  }
  b_samples
}

#' Evaluate calibration function for each Monte Carlo sample
#'
#' For each sample, applies \code{func} with the corresponding row of
#' \code{b_samples} to the corresponding row of \code{meas_samples}.
#'
#' @param meas_samples matrix as returned by
#'   \code{\link{b_sample_meas_data_mc}}.
#' @param b_samples matrix as returned by \code{\link{b_least_mc}}.
#' @param func calibration function.
#'
#' @return A matrix of shape \code{[nsamples, m]} of assigned values.
#' @export
b_eval_mc <- function(meas_samples, b_samples, func) {
  nsamples  <- nrow(meas_samples)
  x_samples <- matrix(0, nrow = nsamples, ncol = ncol(meas_samples))
  for (i in seq_len(nsamples)) {
    x_samples[i, ] <- func(meas_samples[i, ], b_samples[i, ])$x
  }
  x_samples
}

#' Mean and covariance of Monte Carlo samples
#'
#' @param samples numeric matrix where each row is one sample.
#' @return A list with elements \code{mean} and \code{cov}.
#' @export
b_mean_cov_mc <- function(samples) {
  list(mean = colMeans(samples), cov = cov(samples))
}

#' Display calibration data summary from Monte Carlo samples
#'
#' @param cal_samples array as returned by
#'   \code{\link{b_sample_cal_data_mc}}.
#' @return Invisibly \code{NULL}.
#' @export
b_disp_cal_data_mc <- function(cal_samples) {
  b_disp_cal_data(.b_cal_samples_to_cal_data(cal_samples))
}

#' Display calibration results from Monte Carlo samples
#'
#' @param b_samples matrix as returned by \code{\link{b_least_mc}}.
#' @return Invisibly \code{NULL}.
#' @export
b_disp_cal_results_mc <- function(b_samples) {
  mc <- b_mean_cov_mc(b_samples)
  b_disp_cal_results(mc$mean, mc$cov, NA_real_)
}

#' Display measurement results from Monte Carlo samples
#'
#' @param x_samples matrix as returned by \code{\link{b_eval_mc}}.
#' @param meas_samples matrix as returned by
#'   \code{\link{b_sample_meas_data_mc}}.
#' @return Invisibly \code{NULL}.
#' @export
b_disp_meas_results_mc <- function(x_samples, meas_samples) {
  mc       <- b_mean_cov_mc(x_samples)
  meas_sum <- cbind(y  = colMeans(meas_samples),
                    uy = apply(meas_samples, 2, sd))
  b_disp_meas_results(mc$mean, mc$cov, meas_sum)
}

#' Full Monte Carlo calibration and evaluation with printed output
#'
#' @param cal_data numeric matrix as returned by \code{\link{b_read_cal_data}}.
#' @param meas_data numeric matrix as returned by
#'   \code{\link{b_read_meas_data}}.
#' @param func calibration function.
#' @param nsamples integer number of Monte Carlo samples (default 10000).
#'
#' @return Invisibly a list with elements \code{b_samples} and
#'   \code{x_samples}.
#' @export
b_test_mc <- function(cal_data, meas_data, func, nsamples = 10000L) {
  cal_samples  <- b_sample_cal_data_mc(cal_data,  nsamples)
  meas_samples <- b_sample_meas_data_mc(meas_data, nsamples)
  b_disp_cal_data_mc(cal_samples)
  b_samples  <- b_least_mc(cal_samples, func)
  b_disp_cal_results_mc(b_samples)
  x_samples  <- b_eval_mc(meas_samples, b_samples, func)
  b_disp_meas_results_mc(x_samples, meas_samples)
  invisible(list(b_samples = b_samples, x_samples = x_samples))
}

# Convert 3D calibration sample array to a cal_data matrix (mean + sd).
.b_cal_samples_to_cal_data <- function(cal_samples) {
  cal_data       <- matrix(0, nrow = dim(cal_samples)[2], ncol = 4)
  cal_data[, 1]  <- colMeans(cal_samples[, , 1])
  cal_data[, 2]  <- apply(cal_samples[, , 1], 2, sd)
  cal_data[, 3]  <- colMeans(cal_samples[, , 2])
  cal_data[, 4]  <- apply(cal_samples[, , 2], 2, sd)
  cal_data
}
