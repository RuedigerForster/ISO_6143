#' ISO6143.2001: Calibration of Gas Analysers per ISO 6143:2001
#'
#' Implements the B LEAST calibration algorithm of ISO 6143:2001.
#' See \code{\link{b_least}} for the main fitting function and
#' \code{\link{b_eval}} for evaluating fitted functions at measurement points.
#'
#' @keywords internal
#' @importFrom stats coef cov lm poly rnorm sd
#' @importFrom utils read.csv
#' @importFrom grDevices adjustcolor
#' @importFrom graphics arrows grid legend lines par points polygon
"_PACKAGE"
