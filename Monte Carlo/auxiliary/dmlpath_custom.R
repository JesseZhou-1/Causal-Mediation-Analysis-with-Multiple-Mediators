############################################################
## dmlpath_custom.R
## Thin wrapper around modified cmedR dmlpath()/dmlmed()
##
## Provides a simple interface for experiment scripts while
## delegating all estimation to the cmedR implementation
## (modified for count/discrete treatment D).
##
## Supports three methods:
##   "mr2" - Multiply Robust Type 2 (IPW + Regression, default)
##   "ipw" - Pure IPW with SuperLearner propensities
##   "reg" - Pure Regression Imputation with SuperLearner
############################################################

# Load dependencies for cmedR functions
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(caret)
  library(SuperLearner)
})

# Source the modified cmedR estimators
source("dmlmed.R")
source("dmlpath.R")

#' DML Estimator for Path-Specific Effects (wrapper)
#'
#' @param data   data.frame with all variables
#' @param D      character: treatment column name (count/discrete)
#' @param M      list of character: mediator names in causal order
#' @param Y      character: outcome column name (treated as numeric)
#' @param C      character vector: covariate column names
#' @param d      numeric: treatment value
#' @param dstar  numeric: reference treatment value
#' @param n_folds integer: number of cross-fitting folds (default 5)
#' @param SL_library character vector: SuperLearner library
#' @param trim_low  numeric: lower quantile for weight censoring
#' @param trim_high numeric: upper quantile for weight censoring
#' @param seed   integer or NULL: random seed for fold assignment
#' @param method character: "mr2" (multiply robust), "ipw", or "reg"
#'
#' @return list with ATE, PSE vector, and method info
dmlpath_custom <- function(
    data,
    D = "D",
    M = list("M1", "M2"),
    Y = "Y",
    C = "C",
    d = 2,
    dstar = 0,
    n_folds = 5,
    SL_library = c("SL.mean", "SL.glmnet", "SL.ranger"),
    trim_low = 0.01,
    trim_high = 0.99,
    seed = NULL,
    method = "mr2"
) {
  stopifnot(method %in% c("mr2", "ipw", "reg"))

  # Call the modified cmedR dmlpath
  result <- dmlpath(
    data = data,
    D = D,
    M = M,
    Y = Y,
    C = C,
    d = d,
    dstar = dstar,
    num_folds = n_folds,
    V = 5L,
    seed = seed,
    SL.library = SL_library,
    stratifyCV = TRUE,
    censor = TRUE,
    censor_low = trim_low,
    censor_high = trim_high,
    method = method
  )

  # Convert tibble format to list format for backward compatibility
  get_val <- function(pattern) {
    idx <- grep(pattern, result$Estimand)
    if (length(idx) > 0) {
      as.numeric(unname(result$Mean[idx[1]]))
    } else {
      NA_real_
    }
  }

  ate <- get_val("^ATE")

  # Extract PSE values and clean names:
  #   "PSE:D->Y(2,0)" -> "D->Y"
  pse_idx <- grep("^PSE:", result$Estimand)
  pse_vals <- as.numeric(unname(result$Mean[pse_idx]))
  pse_names <- result$Estimand[pse_idx]
  clean_names <- gsub("^PSE:", "", gsub("\\([^)]+\\)$", "", pse_names))
  clean_names[is.na(clean_names) | clean_names == ""] <- paste0("PSE_", seq_along(clean_names))
  names(pse_vals) <- clean_names

  list(
    ATE = ate,
    PSE = pse_vals,
    n_folds = n_folds,
    method = method
  )
}
