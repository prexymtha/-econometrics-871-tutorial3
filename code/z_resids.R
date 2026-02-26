#' ============================================================================
#' z_resids.R
#' ----------------------------------------------------------------------------
#' Plot standardised (z-score) residuals from a fitted model.
#'
#' The standardised residual for observation i is:
#'
#'   z_i = (e_i - mean(e)) / sd(e)
#'
#' where e_i is the i-th raw residual.  Horizontal reference lines at
#' +/- 2 highlight observations that may be outliers (roughly outside the
#' 95 % region under normality).
#'
#' @param model A fitted model object that has a `$residuals` component
#'              (e.g. output from `lm()` or `dynlm()`).
#'
#' @return Invisible NULL – called for its side-effect (a base-R plot).
#' ============================================================================

z_resids <- function(model) {

    # --- Extract and standardise residuals -----------------------------------
    r <- as.numeric(model$residuals)

    # z_i = (e_i - mean(e)) / sd(e)
    r_std <- (r - mean(r)) / sd(r)

    # --- Set up the plot canvas ----------------------------------------------
    n <- length(r_std)

    par(las = 1, bty = "o", cex = 0.75)
    plot(
        r_std,
        pch  = "",
        col  = "red",
        main = "Standardised residuals",
        xlab = "Index",
        ylab = "z"
    )

    # --- Draw vertical segments from zero to each residual -------------------
    segments(
        x0  = seq_len(n),
        y0  = 0,
        x1  = seq_len(n),
        y1  = r_std,
        col = "red"
    )

    # --- Reference lines -----------------------------------------------------
    abline(h = 0)
    abline(h = c(2, -2), lty = 2, col = "blue")

    invisible(NULL)
}
