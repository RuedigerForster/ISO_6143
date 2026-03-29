# Plot calibration curve and data points.

#' Plot calibration curve with uncertainty bands
#'
#' Draws the fitted calibration curve with a \eqn{k=2} coverage band,
#' calibration reference points with error bars and uncertainty ellipses, and
#' measurement points with their propagated uncertainties.
#'
#' Colour scheme: axes and grid in steel blue, fitted curve in poison green
#' (\code{"#39FF14"}), data points and uncertainty indicators in dark orange.
#' Calibration reference points are drawn as filled circles (\code{pch = 16});
#' measurement points as filled triangles (\code{pch = 17}).
#'
#' @param cal_data numeric matrix as returned by \code{\link{b_read_cal_data}}.
#' @param meas_data numeric matrix as returned by
#'   \code{\link{b_read_meas_data}}.
#' @param b numeric vector of fitted coefficients from \code{\link{b_least}}.
#' @param b_cov covariance matrix of \code{b} from \code{\link{b_least}}.
#' @param func calibration function used during fitting.
#'
#' @return Invisibly \code{NULL}.
#' @export
b_plot <- function(cal_data, meas_data, b, b_cov, func) {
  k         <- 2
  col_fit   <- "#39FF14"     # poison green — fitted curve
  col_ax    <- "steelblue"   # axes, grid, frame
  col_data  <- "darkorange"  # data points and uncertainties
  # Fit curve over the response range
  ymin   <- min(cal_data[, 3], meas_data[, 1])
  ymax   <- max(cal_data[, 3], meas_data[, 1])
  fy     <- seq(ymin, ymax, length.out = 100)
  f_data <- cbind(fy, 0)
  fit    <- b_eval(f_data, b, b_cov, func)
  fx     <- fit$x
  ufx    <- sqrt(diag(fit$x_cov))
  # Measurement results with full xy covariance
  xy      <- b_eval_xy(meas_data, b, b_cov, func)
  x_meas  <- xy$x
  xc_meas <- xy$x_cov
  yc_meas <- xy$y_cov
  xyc     <- xy$xy_cov
  # Axis limits
  xlim <- range(c(fx - k * ufx, fx + k * ufx,
                  cal_data[, 1] - k * cal_data[, 2],
                  cal_data[, 1] + k * cal_data[, 2],
                  x_meas - k * sqrt(diag(xc_meas)),
                  x_meas + k * sqrt(diag(xc_meas))))
  ylim <- range(c(fy,
                  cal_data[, 3] - k * cal_data[, 4],
                  cal_data[, 3] + k * cal_data[, 4],
                  meas_data[, 1] - k * sqrt(diag(yc_meas)),
                  meas_data[, 1] + k * sqrt(diag(yc_meas))))
  # Apply colour theme to frame and axis text; restore on exit
  old_par <- par(fg = col_ax, col.axis = col_ax, col.lab = col_ax)
  on.exit(par(old_par), add = TRUE)
  # Blank canvas
  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = "Assigned value x", ylab = "Instrument response y")
  # Grid
  grid(col = adjustcolor(col_ax, alpha.f = 0.4), lty = "dotted", lwd = 0.8)
  # Uncertainty band around fit curve
  polygon(c(fx - k * ufx, rev(fx + k * ufx)),
          c(fy, rev(fy)),
          col = adjustcolor(col_fit, alpha.f = 0.2), border = NA)
  lines(fx, fy, col = col_fit, lwd = 2)
  # Calibration reference points (circles)
  for (i in seq_len(nrow(cal_data))) {
    cv_i <- matrix(c(cal_data[i, 2]^2, 0, 0, cal_data[i, 4]^2), 2, 2)
    .b_plot_ellipse(cal_data[i, 1], cal_data[i, 3], cv_i, col_data)
  }
  arrows(cal_data[, 1] - k * cal_data[, 2], cal_data[, 3],
         cal_data[, 1] + k * cal_data[, 2], cal_data[, 3],
         angle = 90, code = 3, length = 0.05, col = col_data)
  arrows(cal_data[, 1], cal_data[, 3] - k * cal_data[, 4],
         cal_data[, 1], cal_data[, 3] + k * cal_data[, 4],
         angle = 90, code = 3, length = 0.05, col = col_data)
  points(cal_data[, 1], cal_data[, 3], pch = 16, col = col_data)
  # Measurement points (triangles)
  for (i in seq_len(nrow(meas_data))) {
    cv_i <- matrix(c(xc_meas[i, i], xyc[i, i], xyc[i, i], yc_meas[i, i]), 2, 2)
    .b_plot_ellipse(x_meas[i], meas_data[i, 1], cv_i, col_data)
  }
  arrows(x_meas - k * sqrt(diag(xc_meas)), meas_data[, 1],
         x_meas + k * sqrt(diag(xc_meas)), meas_data[, 1],
         angle = 90, code = 3, length = 0.05, col = col_data)
  arrows(x_meas, meas_data[, 1] - k * sqrt(diag(yc_meas)),
         x_meas, meas_data[, 1] + k * sqrt(diag(yc_meas)),
         angle = 90, code = 3, length = 0.05, col = col_data)
  points(x_meas, meas_data[, 1], pch = 17, col = col_data)
  legend("topleft",
         legend  = c("Fit x = f(y)", "Reference points", "Measurement points"),
         col     = c(col_fit, col_data, col_data),
         lty     = c(1, NA, NA),
         pch     = c(NA, 16, 17),
         text.col = col_ax,
         box.col  = col_ax)
  invisible(NULL)
}

# Draw a filled uncertainty ellipse (k = 2.45 coverage factor).
.b_plot_ellipse <- function(px, py, cv, color) {
  k   <- 2.45
  eig <- eigen(k^2 * cv)
  d   <- pmax(eig$values, 0)   # guard against tiny negative eigenvalues
  t   <- seq(0, 2 * pi, length.out = 100)
  e   <- eig$vectors %*% diag(sqrt(d))
  g   <- e %*% rbind(cos(t), sin(t))
  polygon(g[1, ] + px, g[2, ] + py,
          col = adjustcolor(color, alpha.f = 0.4), border = NA)
}
