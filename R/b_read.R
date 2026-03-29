# Read calibration / measurement data from tab-separated text files.

#' Read calibration data from a text file
#'
#' Reads a tab-separated file with four columns: assigned value \eqn{x},
#' standard uncertainty \eqn{u(x)}, instrument response \eqn{y}, standard
#' uncertainty \eqn{u(y)}.  No header line is expected.
#'
#' @param filepath character string giving the path to the file.
#' @param sep field separator; defaults to \code{"\\t"}.
#'
#' @return A numeric matrix with \eqn{n} rows and four columns named
#'   \code{x}, \code{ux}, \code{y}, \code{uy}.
#' @export
b_read_cal_data <- function(filepath, sep = "\t") {
  as.matrix(read.table(filepath, header = FALSE, sep = sep,
                       col.names = c("x", "ux", "y", "uy")))
}

#' Read measurement data from a text file
#'
#' Reads a tab-separated file with two columns: instrument response \eqn{y}
#' and standard uncertainty \eqn{u(y)}.  No header line is expected.
#'
#' @param filepath character string giving the path to the file.
#' @param sep field separator; defaults to \code{"\\t"}.
#'
#' @return A numeric matrix with \eqn{m} rows and two columns named
#'   \code{y} and \code{uy}.
#' @export
b_read_meas_data <- function(filepath, sep = "\t") {
  as.matrix(read.table(filepath, header = FALSE, sep = sep,
                       col.names = c("y", "uy")))
}

#' Path to a package example data file
#'
#' Returns the full path to one of the six example data files from Annex B
#' of ISO 6143:2001 that are shipped with the package.
#'
#' @param example integer, 1, 2, or 3.
#' @param type character, \code{"cal"} for calibration data or
#'   \code{"meas"} for measurement data.
#'
#' @return A character string giving the full path to the file.
#' @export
b_example_data_path <- function(example, type) {
  stopifnot(example %in% 1:3, type %in% c("cal", "meas"))
  system.file("extdata",
              sprintf("b_least_%d_data_%s.txt", example, type),
              package = "ISO6143.2001",
              mustWork = TRUE)
}
