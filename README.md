# ISO6143.2001

An R implementation of the B LEAST calibration algorithm defined in
**ISO 6143:2001** — *Gas analysis — Comparison methods for determining and
checking the composition of calibration gas mixtures*.

## Overview

The package fits polynomial, power, and exponential calibration functions to
reference gas standards using errors-in-both-variables least squares
(Levenberg–Marquardt), propagates measurement uncertainty through the fitted
function via analytical Jacobians, and provides optional Monte Carlo
verification. Numerical examples from Annex B of the standard are included.

The algorithm is a port of the
[METAS B LEAST Python implementation](https://github.com/METAS-MSL/metas-b-least)
by Michael Wollensack (METAS).

## Disclaimer

> **This software is provided without warranty of any kind.** The authors
> accept no responsibility for results obtained with it. Users are solely
> responsible for verifying that outputs are fit for their intended purpose.

The implementation deviates in parts from the published standard. All known
deviations are flagged with `# DEVIATION` comments in the source code and
documented in the accompanying note in each affected function.

## Installation

The package is not yet on CRAN. Install directly from GitHub:

```r
# install.packages("remotes")   # if not already installed
remotes::install_github("RuedigerForster/ISO_6143")
```

Or clone the repository and install from a local source:

```r
# install.packages("devtools")  # if not already installed
devtools::install("path/to/ISO_6143")
```

## Quick start

```r
library(ISO6143.2001)

# Load one of the three Annex B example datasets (1, 2, or 3)
cal <- b_read_cal_data(b_example_data_path(1, "cal"))

# Fit a linear calibration function
fit <- b_least(cal, b_linear_func)
print(fit)
plot(fit)

# Evaluate the fitted function at a measured value
meas <- b_read_meas_data(b_example_data_path(1, "meas"))
result <- b_eval(fit, meas)
print(result)
```

## Dependencies

| Package | Role |
|---------|------|
| `minpack.lm` | Levenberg–Marquardt solver (MINPACK back-end, same as SciPy) |

## License

GPL (≥ 3) — see the `LICENSE` file.
