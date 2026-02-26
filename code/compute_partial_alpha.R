#' ============================================================================
#' compute_partial_alpha.R
#' ----------------------------------------------------------------------------
#' Compute the sample *partial* autocorrelation at lag s via OLS.
#'
#' The partial autocorrelation at lag s isolates the *direct* linear
#' relationship between y_t and y_{t-s} after removing the influence of
#' intervening lags y_{t-1}, ..., y_{t-s+1}.
#'
#' Method (OLS regression approach):
#'   1.  Regress y_t on y_{t-1}, y_{t-2}, ..., y_{t-s}.
#'   2.  The partial autocorrelation at lag s is the OLS coefficient on y_{t-s}.
#'
#' Formally, estimate:
#'
#'   y_t = phi_1 * y_{t-1} + phi_2 * y_{t-2} + ... + phi_s * y_{t-s} + u_t
#'
#' Then  PACF(s) = hat{phi}_s.
#'
#' @param ts_vec A numeric vector (the time series).
#' @param s      Integer lag (s >= 0).
#'
#' @return Scalar partial autocorrelation at lag s, or NA if s >= length(ts_vec).
#' ============================================================================

compute_partial_alpha <- function(ts_vec, s) {

    # By convention, partial autocorrelation at lag 0 is 1
    if (s == 0L) return(1)

    n <- length(ts_vec)

    # Guard against lags longer than the series
    if (s >= n) return(NA_real_)

    # embed() creates a matrix:
    #   column 1 = y_t,  column 2 = y_{t-1}, ..., column s+1 = y_{t-s}
    mat <- embed(ts_vec, s + 1L)

    y <- mat[, 1]                 # dependent variable: y_t
    X <- mat[, -1, drop = FALSE]  # regressors: y_{t-1} ... y_{t-s}

    fit <- lm(y ~ X)

    # Return the coefficient on y_{t-s} (the last regressor)
    # Coefficients vector: (Intercept), X1, X2, ..., Xs  -> element s+1
    unname(coef(fit)[s + 1L])
}
