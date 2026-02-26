#' ============================================================================
#' fit_actual.R
#' ----------------------------------------------------------------------------
#' Plot two time-series objects on the same axes (actual vs. fitted).
#'
#' Purpose:
#'   Overlay the *actual* series (blue) and the *fitted / predicted* series
#'   (red) so we can visually assess goodness-of-fit.
#'
#' @param ts1  A `ts` object – typically the actual (observed) series.
#' @param ts2  A `ts` object – typically the fitted (predicted) series.
#' @param main Character string for the plot title (default "").
#' @param ts1n Legend label for ts1 (default "").
#' @param ts2n Legend label for ts2 (default "").
#'
#' @return     Invisible NULL – called for its side-effect (a base-R plot).
#' ============================================================================

fit_actual <- function(ts1, ts2, main = "", ts1n = "", ts2n = "") {

    # --- Input validation ----------------------------------------------------
    # Use inherits() instead of class() == "ts" — safer with S4/multi-class
    if (!inherits(ts1, "ts")) stop("ts1 must be of class 'ts'.")
    if (!inherits(ts2, "ts")) stop("ts2 must be of class 'ts'.")

    # --- Compute common y-axis limits ----------------------------------------
    y_range <- range(ts1, ts2, na.rm = TRUE)

    # --- Plot actual series (blue) -------------------------------------------
    par(new = FALSE, lwd = 0.1, bty = "o", xaxt = "s")
    plot.ts(
        ts1,
        col  = "blue",
        ylim = y_range,
        ylab = "",
        xlab = "",
        xaxt = "n",
        main = main
    )

    # --- Overlay fitted series (red) -----------------------------------------
    par(new = TRUE)
    plot.ts(
        ts2,
        col  = "red",
        ylim = y_range,
        ylab = "",
        xlab = "",
        xaxt = "n",
        main = ""
    )

    # --- Add legend and zero-line --------------------------------------------
    legend(
        x      = par("usr")[1],
        y      = par("usr")[4],
        legend = c(ts1n, ts2n),
        lty    = c(1, 1),
        lwd    = c(1, 1),
        cex    = 0.75,
        col    = c("blue", "red")
    )

    abline(h = 0)

    invisible(NULL)
}
