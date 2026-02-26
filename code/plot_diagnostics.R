#' ============================================================================
#' plot_diagnostics.R
#' ----------------------------------------------------------------------------
#' Produce a 2 x 2 panel of diagnostic plots for a time-series regression:
#'
#'   (1) Normal Q-Q plot        – checks normality of residuals
#'   (2) Actual vs. Fitted      – visual goodness-of-fit
#'   (3) Standardised residuals – flags outliers (|z| > 2)
#'   (4) ACF of residuals       – checks serial correlation
#'
#' Requires: fit_actual() and z_resids() to be sourced beforehand.
#'
#' @param model  A fitted model object (e.g. from `lm()` or `dynlm()`).
#' @param actual A `ts` object – the original dependent variable.
#' @param title  Character string; overall panel title (default "Model Diagnostics").
#'
#' @return Invisible NULL – called for its plotting side-effect.
#' ============================================================================

plot_diagnostics <- function(model, actual, title = "Model Diagnostics") {

    # Save current graphics parameters so we can restore them
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)

    par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))

    # Panel 1: Normal Q-Q plot
    plot(model, which = 2)

    # Panel 2: Actual vs. Fitted
    fit_actual(
        actual,
        model$fitted.values,
        main = "Actual and Fitted",
        ts1n = "Actual",
        ts2n = "Fitted"
    )

    # Panel 3: Standardised residuals
    z_resids(model)

    # Panel 4: ACF of residuals
    acf(model$residuals, main = "ACF of Residuals")

    # Overall title
    mtext(title, outer = TRUE, cex = 1.2, font = 2)

    invisible(NULL)
}
