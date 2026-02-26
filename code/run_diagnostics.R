#' ============================================================================
#' run_diagnostics.R
#' ----------------------------------------------------------------------------
#' Run a standard battery of misspecification tests on a fitted linear model
#' and return a tidy summary data-frame.
#'
#' Tests included:
#'
#'   1. Breusch–Godfrey (BG) test for serial correlation up to order `bg_order`
#'      H0: no serial correlation in residuals up to lag p.
#'
#'   2. Breusch–Pagan (BP) test for heteroskedasticity
#'      H0: residual variance is constant (homoskedastic).
#'
#'   3. Jarque–Bera (JB) test for normality
#'      H0: residuals are normally distributed.
#'      JB = (T/6) * [ S^2 + (1/4)(K - 3)^2 ]
#'      where S = skewness, K = kurtosis.
#'
#'   4. RESET test for functional-form misspecification
#'      H0: the model is correctly specified (no omitted non-linearities).
#'
#' @param model    A fitted model object (e.g. from `lm()` or `dynlm()`).
#' @param bg_order Integer; maximum lag order for the BG test (default 5).
#'
#' @return A data.frame with columns: Test, Statistic, p_value.
#' ============================================================================

run_diagnostics <- function(model, bg_order = 5L) {

    # --- Load required packages (should already be loaded) -------------------
    requireNamespace("lmtest",  quietly = TRUE)
    requireNamespace("moments", quietly = TRUE)

    # --- Breusch–Godfrey (serial correlation) --------------------------------
    bg <- lmtest::bgtest(formula = model, order = bg_order, type = "F")

    # --- Breusch–Pagan (heteroskedasticity) ----------------------------------
    bp <- lmtest::bptest(formula = model)

    # --- Jarque–Bera (normality) ---------------------------------------------
    jb <- moments::jarque.test(as.vector(model$residuals))

    # --- RESET (functional form) ---------------------------------------------
    rt <- lmtest::resettest(model)

    # --- Assemble results ----------------------------------------------------
    results <- data.frame(
        Test      = c(
            paste0("Breusch-Godfrey (order ", bg_order, ")"),
            "Breusch-Pagan",
            "Jarque-Bera",
            "RESET"
        ),
        Statistic = round(c(bg$statistic, bp$statistic, jb$statistic, rt$statistic), 4),
        p_value   = round(c(bg$p.value, bp$p.value, jb$p.value, rt$p.value), 4),
        stringsAsFactors = FALSE
    )

    rownames(results) <- NULL
    results
}
