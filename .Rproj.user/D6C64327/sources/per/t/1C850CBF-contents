#' ============================================================================
#' main.R
#' ============================================================================
#' Econometrics 871 – Tutorial 3
#' Consequences of Non-Stationarity and Unit Root Testing
#'
#' Master setup script:
#'   • Clears the workspace
#'   • Installs / loads all required packages
#'   • Sources helper functions from code/
#'   • Reads and prepares the two datasets from data/ (ukm1q, fulton)
#'
#' Run this file FIRST, or let the R Markdown document source it.
#' ============================================================================


# =========================================================================== #
# 0.  House-keeping                                                           #
# =========================================================================== #
# Clear the console and the environment
rm(list = ls())  # Remove all objects from the environment
graphics.off()  # Turn off all graphics devices
cat('\f')  # Clear the console ( in the original code not in markdown because it doesn't clear when knitting )



# =========================================================================== #
# 1.  Package management (using pacman for convenience)                       #
# =========================================================================== #
#
#   dynlm      – dynamic linear models with ts lag/diff/trend notation
#   lmtest     – Breusch-Godfrey, Breusch-Pagan, RESET, LR tests
#   moments    – Jarque-Bera test, skewness, kurtosis
#   fUnitRoots – automated Augmented Dickey-Fuller (ADF) tests
#   xts / zoo  – extended time-series classes and utilities
#   ggplot2    – grammar-of-graphics plotting
#   ggfortify  – autoplot() for ts objects
#   viridis    – colour-blind-friendly palettes
#   reshape2   – melt() for wide -> long conversion (used in correlogram plots)
#   gridExtra  – grid.arrange() to combine ggplots
#
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(
    dynlm,
    lmtest,
    moments,
    fUnitRoots,
    xts,
    zoo,
    ggplot2,
    ggfortify,
    viridis,
    reshape2,
    gridExtra
)


# =========================================================================== #
# 2.  Source helper functions  from the code folder                                               #
# =========================================================================== #
#
#   fit_actual.R           – overlay actual vs. fitted time series
#   z_resids.R             – standardised-residual spike plot
#   compute_alpha.R        – sample autocorrelation at lag s
#   compute_partial_alpha.R – sample partial autocorrelation at lag s (OLS)
#   run_diagnostics.R      – consolidated misspecification-test table
#   plot_diagnostics.R     – 2x2 diagnostic panel
#
source("code/fit_actual.R")
source("code/z_resids.R")
source("code/compute_alpha.R")
source("code/compute_partial_alpha.R")
source("code/run_diagnostics.R")
source("code/plot_diagnostics.R")


# =========================================================================== #
# 3.  Load datasets                                                           #
# =========================================================================== #
#
#   ukm1q  – UK macroeconomic data (quarterly, 1963 Q1 – 1989 Q4)
#             Key variables:
#               x85a  – log total expenditure (1985 prices)
#               pa    – log GDP deflator
#
#   fulton – Fulton fish market data (daily)
#             Key variables:
#               LogPrice – log of fish price
#               Stormy   – dummy for stormy weather
#               Mixed    – dummy for mixed weather
#

ukm1q  <- read.delim2("data/ukm1q.txt")
fulton <- read.delim2("data/fulton.txt")


# =========================================================================== #
# 4.  Create xts objects (for xts-based operations later)                     #
# =========================================================================== #

# Fulton: convert Date column to proper Date class; drop index & date columns
fulton_xts <- xts(
    fulton[, c(-1, -6)],
    order.by = as.Date(as.character(fulton$Date), format = "%y%m%d")
)

# UK M1: quarterly index
ukm1q_xts <- xts(
    ukm1q[, -1],
    order.by = as.yearqtr(ukm1q$date)
)


# =========================================================================== #
# 5.  Create core ts objects for dynlm                                        #
# =========================================================================== #
#
#   dynlm requires objects of class `ts` (not xts).
#   We drop the last observation of ukm1q because it is NA.

x85a  <- ts(ukm1q$x85a[-108], start = c(1963, 1), frequency = 4)
dx85a <- diff(x85a)            # first difference: Delta x85a_t
pa    <- ts(ukm1q$pa[-108],   start = c(1963, 1), frequency = 4)


cat("\n--- Tutorial 3 environment ready ---\n")
cat("Objects loaded: ukm1q, fulton, fulton_xts, ukm1q_xts, x85a, dx85a, pa\n")
cat("Helper functions: fit_actual, z_resids, compute_alpha,",
    "compute_partial_alpha, run_diagnostics, plot_diagnostics\n\n")

