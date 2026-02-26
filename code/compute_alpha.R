#' ============================================================================
#' compute_alpha.R
#' ----------------------------------------------------------------------------
#' Compute the sample autocorrelation (correlogram) at lag s.
#'
#' Definition (sample autocorrelation at lag s):
#'
#'   alpha(s) = Cov(y_t, y_{t-s}) / [ SD(y_t) * SD(y_{t-s}) ]
#'
#' where the covariance and standard deviations are computed over the
#' overlapping portion of the series.
#'
#' NOTE: For a stationary series the population autocorrelation simplifies to
#'   rho(s) = gamma(s) / gamma(0)
#' because SD(y_t) = SD(y_{t-s}).  The function below uses the more general
#' formula (ratio of covariance to product of SDs) which also works for
#' non-stationary series where the two sub-vectors may differ in spread.
#'
#' @param ts_vec A numeric vector (the time series).
#' @param s      Integer lag (s >= 0).
#'
#' @return Scalar autocorrelation at lag s, or NA if s >= length(ts_vec).
#' ============================================================================

compute_alpha <- function(ts_vec, s) {

    # By convention, autocorrelation at lag 0 is 1
    if (s == 0L) return(1)

    n <- length(ts_vec)

    # Guard against lags longer than the series
    if (s >= n) return(NA_real_)

    # x = y_{s+1}, ..., y_T   (the "current" portion)
    # y = y_1, ..., y_{T-s}   (the "lagged" portion)
    x <- ts_vec[(s + 1L):n]
    y <- ts_vec[1L:(n - s)]

    cov_xy <- cov(x, y, use = "complete.obs")
    sd_x   <- sd(x, na.rm = TRUE)
    sd_y   <- sd(y, na.rm = TRUE)

    # Avoid division by zero
    if (sd_x == 0 || sd_y == 0) return(NA_real_)

    cov_xy / (sd_x * sd_y)
}
