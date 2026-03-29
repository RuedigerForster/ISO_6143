# Tests for B LEAST (ISO 6143:2001), mirroring the METAS Python test suite.
# Reference values are from Annex B of the standard and from
# metas-b-least/tests/test_example_*.py.
#
# Tolerances are expressed as maximum absolute deviation, converted to
# relative tolerance for expect_equal() where the sign of the expected value
# matters.  Helper:
near <- function(actual, expected, abs_tol) {
  testthat::expect_true(
    abs(actual - expected) <= abs_tol,
    label = sprintf("abs(%g - %g) = %g, tol = %g",
                    actual, expected, abs(actual - expected), abs_tol)
  )
}

# ---------------------------------------------------------------------------
# Example 1  (Annex B.2.1, pages 23-25): linear function
# ---------------------------------------------------------------------------
test_that("Example 1 linear: coefficients and uncertainties", {
  cal_data  <- b_read_cal_data(b_example_data_path(1, "cal"))
  meas_data <- b_read_meas_data(b_example_data_path(1, "meas"))
  res       <- b_least(cal_data, b_linear_func)
  b         <- res$b
  b_cov     <- res$b_cov

  near(b[1],              -3.5747e-1, 0.0001e-1)
  near(b[2],               2.4612e1,  0.0001e1)
  near(sqrt(b_cov[1, 1]),  1.5716e-1, 0.0003e-1)
  near(sqrt(b_cov[2, 2]),  4.8048e-1, 0.0013e-1)
  near(b_cov[1, 2],       -5.6921e-2, 0.0031e-2)
})

test_that("Example 1 linear: residuals", {
  cal_data <- b_read_cal_data(b_example_data_path(1, "cal"))
  res      <- b_least(cal_data, b_linear_func)
  near(sum(res$b_res^2),    0.6743, 0.0001)
  near(max(abs(res$b_res)), 0.568,  0.001)
})

test_that("Example 1 linear: assigned values and covariances", {
  cal_data  <- b_read_cal_data(b_example_data_path(1, "cal"))
  meas_data <- b_read_meas_data(b_example_data_path(1, "meas"))
  res       <- b_least(cal_data, b_linear_func)
  ev        <- b_eval(meas_data, res$b, res$b_cov, b_linear_func)
  x         <- ev$x
  xc        <- ev$x_cov

  near(x[1],          5.9923,    0.0001)
  near(sqrt(xc[1,1]), 1.6377e-1, 0.0001e-1)
  near(x[2],          1.4409e1,  0.0001e1)
  near(sqrt(xc[2,2]), 3.5599e-1, 0.0003e-1)
  near(x[3],          4.3943e1,  0.0001e1)
  near(sqrt(xc[3,3]), 1.1631,    0.0002)
  near(xc[1, 2],      1.16e-2,   0.01e-2)
  near(xc[1, 3],      1.48e-2,   0.01e-2)
  near(xc[2, 3],      1.37e-1,   0.01e-1)
})

# ---------------------------------------------------------------------------
# Example 2  (Annex B.2.2, pages 26-28): linear and 2nd-order polynomial
# ---------------------------------------------------------------------------
test_that("Example 2 linear: coefficients and uncertainties", {
  cal_data <- b_read_cal_data(b_example_data_path(2, "cal"))
  res      <- b_least(cal_data, b_linear_func)
  b        <- res$b
  b_cov    <- res$b_cov

  near(b[1],              3.9189e-4,  0.0622e-4)
  near(b[2],              2.4286e-5,  0.0001e-5)
  near(sqrt(b_cov[1, 1]), 1.1458e-3,  0.0002e-3)
  near(sqrt(b_cov[2, 2]), 2.4161e-8,  0.0002e-8)
  near(b_cov[1, 2],      -7.2747e-12, 0.0098e-12)
})

test_that("Example 2 linear: residuals and assigned values", {
  cal_data  <- b_read_cal_data(b_example_data_path(2, "cal"))
  meas_data <- b_read_meas_data(b_example_data_path(2, "meas"))
  res       <- b_least(cal_data, b_linear_func)
  ev        <- b_eval(meas_data, res$b, res$b_cov, b_linear_func)

  near(sum(res$b_res^2),    6.1697, 0.1253)
  near(max(abs(res$b_res)), 1.6322, 0.0057)
  near(ev$x[1],              1.7004,    0.0001)
  near(sqrt(ev$x_cov[1,1]), 2.0244e-3, 0.0002e-3)
  near(ev$x[2],              8.9863,    0.0005)
  near(sqrt(ev$x_cov[2,2]), 9.9718e-3, 0.0001e-3)
})

test_that("Example 2 second-order poly: coefficients", {
  cal_data <- b_read_cal_data(b_example_data_path(2, "cal"))
  res      <- b_least(cal_data, b_second_order_poly)
  b        <- res$b
  b_cov    <- res$b_cov

  near(b[1],              -1.4037e-4,  0.0927e-4)
  near(sqrt(b_cov[1, 1]),  1.175e-3,   0.001e-3)
  near(b[2],               2.4403e-5,  0.0002e-5)
  near(sqrt(b_cov[2, 2]),  5.901e-8,   0.001e-8)
  near(b[3],              -4.1096e-13, 0.0231e-13)
  near(sqrt(b_cov[3, 3]),  1.895e-13,  0.001e-13)
  near(b_cov[1, 2],       -2.057e-11,  0.002e-11)
  near(b_cov[1, 3],        4.667e-17,  0.003e-17)
  near(b_cov[2, 3],       -1.020e-20,  0.001e-20)
})

test_that("Example 2 second-order poly: residuals and assigned values", {
  cal_data  <- b_read_cal_data(b_example_data_path(2, "cal"))
  meas_data <- b_read_meas_data(b_example_data_path(2, "meas"))
  res       <- b_least(cal_data, b_second_order_poly)
  ev        <- b_eval(meas_data, res$b, res$b_cov, b_second_order_poly)

  near(sum(res$b_res^2),    1.4687, 0.0724)
  near(max(abs(res$b_res)), 0.8678, 0.0014)
  near(ev$x[1],              1.7061,    0.0002)
  near(sqrt(ev$x_cov[1,1]), 3.2910e-3, 0.0006e-3)
  near(ev$x[2],              8.9727,    0.0004)
  near(sqrt(ev$x_cov[2,2]), 1.1762e-2, 0.0001e-2)
})

# ---------------------------------------------------------------------------
# Example 3  (Annex B.2.3, pages 28-30): linear, power, exponential
# ---------------------------------------------------------------------------
test_that("Example 3 linear: residuals only (poor fit)", {
  cal_data <- b_read_cal_data(b_example_data_path(3, "cal"))
  res      <- b_least(cal_data, b_linear_func)
  near(sum(res$b_res^2),    272.6392, 0.0001)
  near(max(abs(res$b_res)),   6.8352, 0.0010)
})

test_that("Example 3 power function: coefficients", {
  cal_data <- b_read_cal_data(b_example_data_path(3, "cal"))
  res      <- b_least(cal_data, b_power_func)
  b        <- res$b
  b_cov    <- res$b_cov

  near(b[1],              1.2128e-1,  0.0002e-1)
  near(sqrt(b_cov[1, 1]), 1.7821e-2,  0.0432e-2)
  near(b[2],              5.1213e-4,  0.0003e-4)
  near(sqrt(b_cov[2, 2]), 2.3693e-5,  0.0657e-5)
  near(b[3],              8.4986e-2,  0.0005e-2)
  near(sqrt(b_cov[3, 3]), 4.9745e-3,  0.1428e-3)
  near(b_cov[1, 2],      -3.8265e-7,  0.2213e-7)
  near(b_cov[1, 3],       7.9326e-5,  0.4641e-5)
  near(b_cov[2, 3],      -1.1780e-7,  0.0662e-7)
})

test_that("Example 3 power function: residuals and assigned values", {
  cal_data  <- b_read_cal_data(b_example_data_path(3, "cal"))
  meas_data <- b_read_meas_data(b_example_data_path(3, "meas"))
  res       <- b_least(cal_data, b_power_func)
  ev        <- b_eval(meas_data, res$b, res$b_cov, b_power_func)

  near(sum(res$b_res^2),    8.3804, 0.0001)
  near(max(abs(res$b_res)), 1.1594, 0.0001)
  near(ev$x[1],              5.3456,    0.0001)
  near(sqrt(ev$x_cov[1,1]), 1.4141e-2, 0.0041e-2)
})

test_that("Example 3 exponential function: coefficients and residuals", {
  cal_data  <- b_read_cal_data(b_example_data_path(3, "cal"))
  meas_data <- b_read_meas_data(b_example_data_path(3, "meas"))
  res       <- b_least(cal_data, b_exp_func)
  ev        <- b_eval(meas_data, res$b, res$b_cov, b_exp_func)

  near(res$b[1],            -4.8019e1, 0.0577e1)
  near(res$b[2],             4.8024e1, 0.0556e1)
  near(res$b[3],             2.1261e-5, 0.0001e-2)
  near(sum(res$b_res^2),    0.6581, 0.0009)
  near(max(abs(res$b_res)), 0.3552, 0.0023)
  near(ev$x[1],              5.3357, 0.0001)
})
