# B LEAST core: starting values, residuals, Jacobian, covariance, and fit.

# Moore-Penrose pseudoinverse via SVD (base R, no extra dependency).
.pinv <- function(A) {
  s     <- svd(A)
  tol   <- max(dim(A)) * max(s$d) * .Machine$double.eps
  d_inv <- ifelse(s$d > tol, 1 / s$d, 0)
  s$v %*% diag(d_inv) %*% t(s$u)
}

# Compute OLS starting values for the b coefficients.
.b_least_start <- function(cal_data, func) {
  x <- cal_data[, 1]
  y <- cal_data[, 3]
  if (identical(func, b_linear_func)) {
    as.numeric(coef(lm(x ~ y)))
  } else if (identical(func, b_second_order_poly)) {
    as.numeric(coef(lm(x ~ poly(y, 2, raw = TRUE))))
  } else if (identical(func, b_third_order_poly)) {
    as.numeric(coef(lm(x ~ poly(y, 3, raw = TRUE))))
  } else if (identical(func, b_power_func)) {
    c(as.numeric(coef(lm(x ~ y))), 0)
  } else if (identical(func, b_exp_func)) {
    # x = b0 + b1*exp(b2*y) approximated as 2nd-order Taylor:
    # x ≈ (b0+b1) + b1*b2*y + b1*b2^2/2*y^2
    cc <- as.numeric(coef(lm(x ~ poly(y, 2, raw = TRUE))))
    b2 <- 2 * cc[3] / cc[2]
    b1 <- cc[2] / b2
    b0 <- cc[1] - b1
    c(b0, b1, b2)
  } else {
    stop("Unknown fit function. Use one of b_linear_func, b_second_order_poly, ",
         "b_third_order_poly, b_power_func, b_exp_func.")
  }
}

# Scaled residuals passed to nls.lm.
.b_residuals <- function(params, cal_data, y2_b_scale, func) {
  n    <- nrow(cal_data)
  y2_b <- params * y2_b_scale
  y2   <- y2_b[seq_len(n)]
  b    <- y2_b[-seq_len(n)]
  .b_objective_func(cal_data[, 1], cal_data[, 2],
                    cal_data[, 3], cal_data[, 4],
                    y2, b, func)
}

# Jacobian of the scaled residuals (2n x n+nb matrix).
.b_jacobian <- function(params, cal_data, y2_b_scale, func) {
  n        <- nrow(cal_data)
  nb       <- length(y2_b_scale) - n
  y2_b     <- params * y2_b_scale
  y2       <- y2_b[seq_len(n)]
  b        <- y2_b[-seq_len(n)]
  y2_scale <- y2_b_scale[seq_len(n)]
  b_scale  <- y2_b_scale[-seq_len(n)]
  ux       <- cal_data[, 2]
  uy       <- cal_data[, 4]
  f        <- func(y2, b)
  dx2_dy2  <- f$dx_dy
  dx2_db   <- f$dx_db
  jacobi   <- matrix(0, nrow = 2 * n, ncol = n + nb)
  for (i in seq_len(n)) {
    jacobi[i, i] <- dx2_dy2[i] * y2_scale[i] / ux[i]
    for (j in seq_len(nb)) {
      jacobi[i, n + j] <- dx2_db[[j]][i] * b_scale[j] / ux[i]
    }
    jacobi[n + i, i] <- y2_scale[i] / uy[i]
  }
  jacobi
}

# GUM uncertainty propagation for the b coefficients.
.b_covariance <- function(params, cal_data, y2_b_scale, func) {
  n               <- nrow(cal_data)
  b_scale         <- y2_b_scale[-seq_len(n)]
  dg              <- .b_jacobian(params, cal_data, y2_b_scale, func)
  dg_inv          <- .pinv(dg)
  db_dg           <- dg_inv[-seq_len(n), , drop = FALSE]
  j               <- diag(b_scale) %*% (-db_dg)
  j %*% t(j)
}

#' Fit calibration function coefficients by errors-in-both-variables least squares
#'
#' Uses Levenberg-Marquardt (MINPACK) with an analytic Jacobian to fit the
#' coefficients \eqn{b} of \code{func} to the calibration data.  Uncertainty
#' is propagated from the input uncertainties via a GUM-compatible analytic
#' covariance calculation.
#'
#' @param cal_data numeric matrix with \eqn{n} rows and four columns:
#'   assigned value \eqn{x}, standard uncertainty \eqn{u(x)}, instrument
#'   response \eqn{y}, standard uncertainty \eqn{u(y)}.
#' @param func calibration function; one of \code{\link{b_linear_func}},
#'   \code{\link{b_second_order_poly}}, \code{\link{b_third_order_poly}},
#'   \code{\link{b_power_func}}, \code{\link{b_exp_func}}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{b}{numeric vector of fitted coefficients.}
#'     \item{b_cov}{covariance matrix of \code{b}.}
#'     \item{b_res}{weighted residual vector of length \eqn{2n}.}
#'   }
#' @importFrom minpack.lm nls.lm nls.lm.control
#' @export
b_least <- function(cal_data, func) {
  n            <- nrow(cal_data)
  y2_start     <- cal_data[, 3]
  b_start      <- .b_least_start(cal_data, func)
  y2_b_start   <- c(y2_start, b_start)
  y2_b_scale   <- y2_b_start
  y2_b_scale[y2_b_scale == 0] <- 1
  y2_b_start2  <- y2_b_start / y2_b_scale
  lm_result <- nls.lm(
    par     = y2_b_start2,
    fn      = .b_residuals,
    jac     = .b_jacobian,
    control = nls.lm.control(maxiter = 1000),
    cal_data    = cal_data,
    y2_b_scale  = y2_b_scale,
    func        = func
  )
  y2_b_opt <- lm_result$par * y2_b_scale
  b_opt    <- y2_b_opt[-seq_len(n)]
  b_cov    <- .b_covariance(lm_result$par, cal_data, y2_b_scale, func)
  b_res    <- .b_residuals(lm_result$par, cal_data, y2_b_scale, func)
  list(b = b_opt, b_cov = b_cov, b_res = b_res)
}
