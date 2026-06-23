############################################################
## Run Experiment 1: Compare four estimators for PSE
## This script sources all required functions and runs a single
## demonstration of the DGP with all four estimators.
############################################################

## ----------------------------
## 0) Setup and load packages
## ----------------------------
suppressPackageStartupMessages({
  library(stats)
  library(MASS)
  library(nnet)
  library(dplyr)
  library(tidyr)
  library(paths)    # for pathimp
  library(stringr)  # for pathimp output processing
})

cat("Packages loaded successfully.\n\n")

## ----------------------------
## 1) Define helper functions
## ----------------------------

# Utility functions from cmedR
demean <- function(x, w = rep(1, length(x))) x - weighted.mean(x, w, na.rm = TRUE)

trimQ <- function(x, low = 0.01, high = 0.99) {
  min <- quantile(x, low)
  max <- quantile(x, high)
  x[x < min] <- min
  x[x > max] <- max
  x
}

comb_list_vec <- function(...) {
  mapply(c, ..., SIMPLIFY = FALSE)
}

## ----------------------------
## 2) Define linmed_inner (for linpath)
## ----------------------------
linmed_inner <- function(
    data,
    D,
    M,
    Y,
    C = NULL,
    d = 1,
    dstar = 0,
    m = rep(0, length(M)),
    interaction_DM = FALSE,
    interaction_DC = FALSE,
    interaction_MC = FALSE,
    weights_name = NULL,
    minimal = FALSE
) {
  df <- data

  # check for missing data and create missing summary output
  key_vars <- c(D, M, Y, C)
  if (!minimal) {
    miss_summary <- sapply(
      key_vars,
      FUN = function(v) c(
        nmiss = sum(!is.na(df[[v]])),
        miss = sum(is.na(df[[v]]))
      )
    ) |>
      t() |>
      as.data.frame()
  }

  # assign weights
  if (is.null(weights_name)) {
    weights <- rep(1, nrow(df))
  } else {
    weights <- df[[weights_name]]
  }

  # demean covariates
  for(covariate in C) df[[covariate]] <- demean(df[[covariate]], w = weights)

  # mediator model(s) predictors
  m_preds <- paste(c(C, D), collapse = " + ")
  if (interaction_DC) {
    m_preds <- paste(m_preds, "+", paste(D, C, sep = ":", collapse = " + "))
  }

  # outcome model predictors
  y_preds <- paste(c(C, D, M), collapse = " + ")
  if (interaction_DM) {
    y_preds <- paste(y_preds, "+", paste(D, M, sep = ":", collapse = " + "))
  }
  if (interaction_DC) {
    y_preds <- paste(y_preds, "+", paste(D, C, sep = ":", collapse = " + "))
  }
  if (interaction_MC) {
    y_preds <- paste(y_preds, "+", paste(outer(M, C, FUN = "paste", sep = ":"), collapse = " + "))
  }

  # specify formulas
  m_forms <- lapply(M, function(x) paste(x, "~", m_preds))
  y_form <- paste(Y, "~", y_preds)

  # fit mediator and outcome models
  m_models <- lapply(m_forms, function(x) lm(as.formula(x), data = df, weights = weights))
  names(m_models) <- M
  y_model <- lm(as.formula(y_form), data = df, weights = weights)

  # compute effects
  if (interaction_DM) {
    NDE_part <- mapply(
      function(M_k, M_model_k) y_model$coef[[paste0(D, ":", M_k)]] * (M_model_k$coef[["(Intercept)"]] + M_model_k$coef[[D]] * dstar),
      M, m_models
    )
    NIE_part <- mapply(
      function(M_k, M_model_k) M_model_k$coef[[D]] * (y_model$coef[[M_k]] + y_model$coef[[paste0(D, ":", M_k)]] * d),
      M, m_models
    )
    CDE_part <- mapply(
      function(M_k, m_k) y_model$coef[[paste0(D, ":", M_k)]] * m_k,
      M, m
    )
  } else {
    NDE_part <- 0
    NIE_part <- mapply(
      function(M_k, M_model_k) M_model_k$coef[[D]] * y_model$coef[[M_k]],
      M, m_models
    )
    CDE_part <- 0
  }

  NDE <- (y_model$coef[[D]] + sum(NDE_part)) * (d - dstar)
  NIE <- sum(NIE_part) * (d - dstar)
  ATE <- NDE + NIE
  CDE <- (y_model$coef[[D]] + sum(CDE_part)) * (d - dstar)

  if (minimal) {
    out <- list(ATE = ATE, NDE = NDE, NIE = NIE, CDE = CDE)
  } else {
    out <- list(ATE = ATE, NDE = NDE, NIE = NIE, CDE = CDE,
                model_m = m_models, model_y = y_model, miss_summary = miss_summary)
  }
  return(out)
}

## ----------------------------
## 3) Source the main functions
## ----------------------------
source("ipwmed.R")
source("ipwvent.R")
source("linpath.R")
source("ipwpath.R")
source("medsim.R")
source("pathimp.R")

cat("All functions sourced successfully.\n\n")

## ----------------------------
## 5) DGP for Experiment 1
## ----------------------------
# Original parameters (probabilities near 0.5 - linear approximation works well)
# exp1_params <- list(
#   a0 = -0.2, a1 = 0.7,
#   b10 = -0.1, b11 = 0.8, b12 = 0.6,
#   b20 = -0.2, b21 = 0.6, b22 = 0.9, b23 = 0.4,
#   g0 = 0.0, g1 = 0.5, g2 = 0.6, g3 = 0.7, g4 = 0.4,
#   sigma = 1.0
# )

# Modified parameters: shift intercepts to create more extreme probabilities
# where logistic curves are more nonlinear (far from 0.5)
# b10 = -1.5: baseline P(M1=1|D=0,C=0) ≈ 0.18
# b20 = 1.2: baseline P(M2=1|D=0,M1=0,C=0) ≈ 0.77
# This creates asymmetry and exposes nonlinearity
exp1_params <- list(
  a0 = -0.2, a1 = 0.7,                           # Treatment model (unchanged)
  b10 = -1.5, b11 = 1.2, b12 = 0.8,              # M1: low baseline, strong D effect
  b20 = 1.2, b21 = 0.8, b22 = 1.5, b23 = 0.6,   # M2: high baseline, strong M1 effect
  g0 = 0.0, g1 = 0.5, g2 = 0.6, g3 = 0.7, g4 = 0.4,  # Outcome model (unchanged)
  sigma = 1.0
)

invlogit <- function(x) 1 / (1 + exp(-x))

gen_exp1 <- function(n, par) {
  C <- rnorm(n, 0, 1)
  pD <- invlogit(par$a0 + par$a1 * C)
  D <- rbinom(n, 1, pD)
  pM1 <- invlogit(par$b10 + par$b11 * D + par$b12 * C)
  M1 <- rbinom(n, 1, pM1)
  pM2 <- invlogit(par$b20 + par$b21 * D + par$b22 * M1 + par$b23 * C)
  M2 <- rbinom(n, 1, pM2)
  muY <- par$g0 + par$g1 * D + par$g2 * M1 + par$g3 * M2 + par$g4 * C
  Y <- rnorm(n, muY, par$sigma)
  data.frame(C = C, D = D, M1 = M1, M2 = M2, Y = Y)
}

## ----------------------------
## 6) Compute Monte Carlo truth
## ----------------------------
truth_exp1 <- function(N_truth = 500000L, par = exp1_params) {
  C <- rnorm(N_truth)
  pM1 <- function(d, c) invlogit(par$b10 + par$b11 * d + par$b12 * c)
  pM2 <- function(d, m1, c) invlogit(par$b20 + par$b21 * d + par$b22 * m1 + par$b23 * c)
  EY <- function(d, m1, m2, c) par$g0 + par$g1 * d + par$g2 * m1 + par$g3 * m2 + par$g4 * c

  M1_0 <- rbinom(N_truth, 1, pM1(0, C))
  M1_1 <- rbinom(N_truth, 1, pM1(1, C))
  M2_0_00 <- rbinom(N_truth, 1, pM2(0, M1_0, C))
  M2_1_00 <- rbinom(N_truth, 1, pM2(1, M1_0, C))
  M2_1_11 <- rbinom(N_truth, 1, pM2(1, M1_1, C))

  Y_0_00 <- EY(0, M1_0, M2_0_00, C)
  Y_1_11 <- EY(1, M1_1, M2_1_11, C)
  Y_1_00 <- EY(1, M1_0, M2_0_00, C)
  Y_1_01 <- EY(1, M1_0, M2_1_00, C)

  c(
    ATE = mean(Y_1_11 - Y_0_00),
    `D->Y` = mean(Y_1_00 - Y_0_00),
    `D->M2->Y` = mean(Y_1_01 - Y_1_00),
    `D->M1~>Y` = mean(Y_1_11 - Y_1_01)
  )
}

cat("Computing Monte Carlo truth (500K samples)...\n")
set.seed(12345)
truth <- truth_exp1(N_truth = 500000L)
cat("True PSEs:\n")
print(round(truth, 4))
cat("\n")

## ----------------------------
## 7) Generate one dataset
## ----------------------------
set.seed(2024)
n <- 2000
df <- gen_exp1(n, exp1_params)
cat("Generated dataset with n =", n, "observations.\n\n")

## ----------------------------
## 8) Run Linear estimator
## ----------------------------
cat("Running Linear path-specific estimator...\n")
lin_fit <- tryCatch(
  linpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C", boot = FALSE),
  error = function(e) { cat("  Error:", e$message, "\n"); NULL }
)

if (!is.null(lin_fit)) {
  lin_results <- c(
    ATE = as.numeric(lin_fit$ATE),
    `D->Y` = as.numeric(lin_fit$PSE[["D->Y"]]),
    `D->M2->Y` = as.numeric(lin_fit$PSE[["D->M2->Y"]]),
    `D->M1~>Y` = as.numeric(lin_fit$PSE[["D->M1~>Y"]])
  )
  cat("Linear estimates:\n")
  print(round(lin_results, 4))
} else {
  lin_results <- rep(NA, 4)
  names(lin_results) <- c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y")
}
cat("\n")

## ----------------------------
## 9) Run IPW estimator
## ----------------------------
cat("Running IPW path-specific estimator...\n")
ipw_fit <- tryCatch(
  ipwpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C", boot = FALSE),
  error = function(e) { cat("  Error:", e$message, "\n"); NULL }
)

if (!is.null(ipw_fit)) {
  ipw_results <- c(
    ATE = as.numeric(ipw_fit$ATE),
    `D->Y` = as.numeric(ipw_fit$PSE[["D->Y"]]),
    `D->M2->Y` = as.numeric(ipw_fit$PSE[["D->M2->Y"]]),
    `D->M1~>Y` = as.numeric(ipw_fit$PSE[["D->M1~>Y"]])
  )
  cat("IPW estimates:\n")
  print(round(ipw_results, 4))
} else {
  ipw_results <- rep(NA, 4)
  names(ipw_results) <- c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y")
}
cat("\n")

## ----------------------------
## 10) Run Regression Imputation (pathimp)
## ----------------------------
cat("Running Regression Imputation estimator...\n")
reg_fit <- tryCatch({
  y0 <- glm(Y ~ D + C, data = df)
  y1 <- glm(Y ~ D + C + M1, data = df)
  y2 <- glm(Y ~ D + C + M1 + M2, data = df)
  pathimp(
    D = "D", Y = "Y", M = list("M1", "M2"),
    Y_models = list(y0, y1, y2),
    D_model = NULL, data = df,
    boot_reps = 2, out_ipw = FALSE
  )
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(reg_fit)) {
  pure_df <- reg_fit$org_obj$pure
  pure_t1 <- pure_df[pure_df$decomposition == "Type I", ]
  reg_results <- c(
    ATE = pure_t1$estimate[pure_t1$estimand == "total"],
    `D->Y` = pure_t1$estimate[pure_t1$estimand == "direct"],
    `D->M2->Y` = pure_t1$estimate[pure_t1$estimand == "via M2"],
    `D->M1~>Y` = pure_t1$estimate[pure_t1$estimand == "via M1"]
  )
  cat("RegImp estimates:\n")
  print(round(reg_results, 4))
} else {
  reg_results <- rep(NA, 4)
  names(reg_results) <- c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y")
}
cat("\n")

## ----------------------------
## 11) Run MedSim (simulation estimator)
## ----------------------------
cat("Running MedSim (simulation-based) estimator...\n")
med_fit <- tryCatch({
  spec <- list(
    list(func = "glm", formula = M1 ~ D + C, args = list(family = binomial())),
    list(func = "glm", formula = M2 ~ D + M1 + C, args = list(family = binomial())),
    list(func = "lm", formula = Y ~ D + M1 + M2 + C)
  )
  medsim(
    data = df,
    num_sim = 5000,
    treatment = "D",
    intv_med = NULL,
    model_spec = spec,
    seed = 999,
    boot = FALSE
  )
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(med_fit)) {
  # Extract estimates from medsim output
  nm <- names(med_fit)
  is_scalar <- vapply(med_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))
  scalars <- unlist(med_fit[is_scalar])

  grab <- function(pattern) {
    idx <- grep(pattern, names(scalars), fixed = FALSE)
    if (length(idx) == 0) return(NA_real_)
    scalars[idx[1]]
  }

  med_results <- c(
    ATE = grab("TE\\("),
    `D->Y` = grab("PSE.*D\\s*->\\s*Y\\s*\\)"),
    `D->M2->Y` = grab("PSE.*M2\\s*->\\s*Y"),
    `D->M1~>Y` = grab("PSE.*M1\\s*~>\\s*Y")
  )
  cat("MedSim estimates:\n")
  print(round(med_results, 4))
} else {
  med_results <- rep(NA, 4)
  names(med_results) <- c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y")
}
cat("\n")

## ----------------------------
## 12) Summary comparison
## ----------------------------
cat("=" , rep("=", 60), "\n", sep = "")
cat("SUMMARY: Path-Specific Effects Comparison\n")
cat("=" , rep("=", 60), "\n", sep = "")

results_df <- data.frame(
  Estimand = c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y"),
  Truth = round(truth, 4),
  Linear = round(lin_results, 4),
  IPW = round(ipw_results, 4),
  RegImp = round(reg_results, 4),
  MedSim = round(med_results, 4)
)
print(results_df, row.names = FALSE)

cat("\nBias relative to truth:\n")
bias_df <- data.frame(
  Estimand = c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y"),
  Linear = round(lin_results - truth, 4),
  IPW = round(ipw_results - truth, 4),
  RegImp = round(reg_results - truth, 4),
  MedSim = round(med_results - truth, 4)
)
print(bias_df, row.names = FALSE)

##############################################################
## PART 2: Interventional Effects (IDE, IIE, OE)
## Intervening on M2 (focal mediator), M1 is treatment-induced confounder
##############################################################

cat("\n\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("PART 2: INTERVENTIONAL EFFECTS (intervening on M2)\n")
cat("=" , rep("=", 60), "\n\n")

## ----------------------------
## 13) Monte Carlo truth for Interventional Effects
## ----------------------------
truth_intv_exp1 <- function(N_truth = 500000L, par = exp1_params) {
  C <- rnorm(N_truth)
  pM1 <- function(d, c) invlogit(par$b10 + par$b11 * d + par$b12 * c)
  pM2 <- function(d, m1, c) invlogit(par$b20 + par$b21 * d + par$b22 * m1 + par$b23 * c)
  EY <- function(d, m1, m2, c) par$g0 + par$g1 * d + par$g2 * m1 + par$g3 * m2 + par$g4 * c

  # Draw M1 under D=0 and D=1
  M1_0 <- rbinom(N_truth, 1, pM1(0, C))
  M1_1 <- rbinom(N_truth, 1, pM1(1, C))

  # For OE: E[Y(1,M1(1),M2(1,M1(1)))] - E[Y(0,M1(0),M2(0,M1(0)))]
  # Natural M2 under D=1, M1=M1_1
  M2_1_natural <- rbinom(N_truth, 1, pM2(1, M1_1, C))
  # Natural M2 under D=0, M1=M1_0
  M2_0_natural <- rbinom(N_truth, 1, pM2(0, M1_0, C))

  Y_1_natural <- EY(1, M1_1, M2_1_natural, C)
  Y_0_natural <- EY(0, M1_0, M2_0_natural, C)

  # For IDE: E[Y(1,M1(1),M2*(0))] - E[Y(0,M1(0),M2(0,M1(0)))]
  # M2*(0) is drawn from marginal distribution of M2|D=0 (integrating over M1)
  # We need to draw M1 from f(M1|D=0) and then M2 from f(M2|D=0,M1), independently for each unit
  # For the counterfactual, we need M2 drawn as if D=0, marginalizing over M1|D=0
  M1_for_M2star <- rbinom(N_truth, 1, pM1(0, C))
  M2_star_0 <- rbinom(N_truth, 1, pM2(0, M1_for_M2star, C))

  # Y(1,M1(1),M2*(0)) - note: M1 is set to M1(1) (natural value under D=1), but M2 is drawn from D=0 distribution
  Y_1_M2star0 <- EY(1, M1_1, M2_star_0, C)

  OE <- mean(Y_1_natural - Y_0_natural)
  IDE <- mean(Y_1_M2star0 - Y_0_natural)
  IIE <- mean(Y_1_natural - Y_1_M2star0)

  c(OE = OE, IDE = IDE, IIE = IIE)
}

cat("Computing Monte Carlo truth for Interventional Effects (500K samples)...\n")
set.seed(12345)
truth_intv <- truth_intv_exp1(N_truth = 500000L)
cat("True Interventional Effects:\n")
print(round(truth_intv, 4))
cat("\n")

## ----------------------------
## 14) IPW for Interventional Effects (using ipwvent.R)
## ----------------------------
cat("Running IPW for Interventional Effects...\n")
ipw_intv_fit <- tryCatch({
  ipwvent_inner(
    data = df,
    D = "D",
    M = "M2",
    Y = "Y",
    L = "M1",
    D_formula = as.formula("D ~ C"),
    L_formula = as.formula("M1 ~ D + C"),
    M_formula = as.formula("M2 ~ D + M1 + C"),
    stabilize = TRUE,
    censor = TRUE
  )
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(ipw_intv_fit)) {
  ipw_intv_results <- c(
    OE = ipw_intv_fit$OE,
    IDE = ipw_intv_fit$IDE,
    IIE = ipw_intv_fit$IIE
  )
  cat("IPW Interventional estimates:\n")
  print(round(ipw_intv_results, 4))
} else {
  ipw_intv_results <- c(OE = NA, IDE = NA, IIE = NA)
}
cat("\n")

## ----------------------------
## 15) RWR (Regression-with-Residuals) for Interventional Effects
## ----------------------------
# Source rwrlite.R (requires rwrmed package)
cat("Running RWR for Interventional Effects...\n")

# Check if rwrmed is installed
rwr_intv_results <- c(OE = NA, IDE = NA, IIE = NA)
if (requireNamespace("rwrmed", quietly = TRUE)) {
  source("rwrlite.R")

  rwr_fit <- tryCatch({
    # For RWR: M_formula should NOT include L (the focal mediator is modeled without L)
    # L_formula_list contains the model for L (M1)
    rwrlite(
      data = df,
      D = "D",
      C = "C",
      Y_formula = as.formula("Y ~ D + M1 + M2 + C + D:M2"),
      M_formula = as.formula("M2 ~ D + C"),  # M2 without M1
      L_formula_list = list(as.formula("M1 ~ D + C")),
      boot = FALSE
    )
  }, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

  if (!is.null(rwr_fit)) {
    rwr_intv_results <- c(
      OE = rwr_fit$OE,
      IDE = rwr_fit$IDE,
      IIE = rwr_fit$IIE
    )
    cat("RWR Interventional estimates:\n")
    print(round(rwr_intv_results, 4))
  }
} else {
  cat("  Note: rwrmed package not installed. Skipping RWR estimator.\n")
  cat("  Install with: devtools::install_github('xiangzhou09/rwrmed')\n")
}
cat("\n")

## ----------------------------
## 16) MedSim for Interventional Effects
## ----------------------------
cat("Running MedSim for Interventional Effects...\n")
medsim_intv_fit <- tryCatch({
  # For interventional effects, the M2 model should NOT include M1
  spec_intv <- list(
    list(func = "glm", formula = M1 ~ D + C, args = list(family = binomial())),
    list(func = "glm", formula = M2 ~ D + C, args = list(family = binomial())),  # M2 without M1
    list(func = "lm", formula = Y ~ D + M1 + M2 + C)
  )
  medsim(
    data = df,
    num_sim = 5000,
    treatment = "D",
    intv_med = "M2",  # specify M2 as the focal mediator for interventional effects
    model_spec = spec_intv,
    seed = 999,
    boot = FALSE
  )
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(medsim_intv_fit)) {
  # Extract estimates
  nm <- names(medsim_intv_fit)
  is_scalar <- vapply(medsim_intv_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))
  scalars <- unlist(medsim_intv_fit[is_scalar])

  grab <- function(pattern) {
    idx <- grep(pattern, names(scalars), fixed = FALSE)
    if (length(idx) == 0) return(NA_real_)
    scalars[idx[1]]
  }

  medsim_intv_results <- c(
    OE = grab("OE\\("),
    IDE = grab("IDE\\("),
    IIE = grab("IIE\\(")
  )
  cat("MedSim Interventional estimates:\n")
  print(round(medsim_intv_results, 4))
} else {
  medsim_intv_results <- c(OE = NA, IDE = NA, IIE = NA)
}
cat("\n")

## ----------------------------
## 17) Summary of Interventional Effects
## ----------------------------
cat("=" , rep("=", 60), "\n", sep = "")
cat("SUMMARY: Interventional Effects Comparison\n")
cat("=" , rep("=", 60), "\n")

intv_results_df <- data.frame(
  Estimand = c("OE", "IDE", "IIE"),
  Truth = round(truth_intv, 4),
  IPW = round(ipw_intv_results, 4),
  RWR = round(rwr_intv_results, 4),
  MedSim = round(medsim_intv_results, 4)
)
print(intv_results_df, row.names = FALSE)

cat("\nBias relative to truth (Interventional Effects):\n")
intv_bias_df <- data.frame(
  Estimand = c("OE", "IDE", "IIE"),
  IPW = round(ipw_intv_results - truth_intv, 4),
  RWR = round(rwr_intv_results - truth_intv, 4),
  MedSim = round(medsim_intv_results - truth_intv, 4)
)
print(intv_bias_df, row.names = FALSE)
cat("\n")

##############################################################
## PART 3: Monte Carlo Simulation (100 reps, 10 cores)
##############################################################

cat("\n\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("PART 3: MONTE CARLO SIMULATION (100 replications, 10 cores)\n")
cat("=" , rep("=", 60), "\n\n")

library(parallel)

# Function to run one replication
run_one_rep <- function(rep_id, n = 2000, par = exp1_params) {
  # Generate data
  df <- gen_exp1(n, par)

  results <- list(rep_id = rep_id)

  # ---- PSE Estimators ----

  # 1. Linear
  lin_fit <- tryCatch({
    linpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C", boot = FALSE)
  }, error = function(e) NULL)

  if (!is.null(lin_fit)) {
    results$lin_ATE <- as.numeric(lin_fit$ATE)
    results$lin_DY <- as.numeric(lin_fit$PSE[["D->Y"]])
    results$lin_DM2Y <- as.numeric(lin_fit$PSE[["D->M2->Y"]])
    results$lin_DM1Y <- as.numeric(lin_fit$PSE[["D->M1~>Y"]])
  } else {
    results$lin_ATE <- results$lin_DY <- results$lin_DM2Y <- results$lin_DM1Y <- NA
  }

  # 2. IPW
  ipw_fit <- tryCatch({
    ipwpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C", boot = FALSE)
  }, error = function(e) NULL)

  if (!is.null(ipw_fit)) {
    results$ipw_ATE <- as.numeric(ipw_fit$ATE)
    results$ipw_DY <- as.numeric(ipw_fit$PSE[["D->Y"]])
    results$ipw_DM2Y <- as.numeric(ipw_fit$PSE[["D->M2->Y"]])
    results$ipw_DM1Y <- as.numeric(ipw_fit$PSE[["D->M1~>Y"]])
  } else {
    results$ipw_ATE <- results$ipw_DY <- results$ipw_DM2Y <- results$ipw_DM1Y <- NA
  }

  # 3. Regression Imputation (pathimp)
  reg_fit <- tryCatch({
    y0 <- glm(Y ~ D + C, data = df)
    y1 <- glm(Y ~ D + C + M1, data = df)
    y2 <- glm(Y ~ D + C + M1 + M2, data = df)
    pathimp(D = "D", Y = "Y", M = list("M1", "M2"),
            Y_models = list(y0, y1, y2), D_model = NULL,
            data = df, boot_reps = 2, out_ipw = FALSE)
  }, error = function(e) NULL)

  if (!is.null(reg_fit)) {
    pure_df <- reg_fit$org_obj$pure
    pure_t1 <- pure_df[pure_df$decomposition == "Type I", ]
    results$reg_ATE <- pure_t1$estimate[pure_t1$estimand == "total"]
    results$reg_DY <- pure_t1$estimate[pure_t1$estimand == "direct"]
    results$reg_DM2Y <- pure_t1$estimate[pure_t1$estimand == "via M2"]
    results$reg_DM1Y <- pure_t1$estimate[pure_t1$estimand == "via M1"]
  } else {
    results$reg_ATE <- results$reg_DY <- results$reg_DM2Y <- results$reg_DM1Y <- NA
  }

  # 4. MedSim (PSE)
  med_fit <- tryCatch({
    spec <- list(
      list(func = "glm", formula = M1 ~ D + C, args = list(family = binomial())),
      list(func = "glm", formula = M2 ~ D + M1 + C, args = list(family = binomial())),
      list(func = "lm", formula = Y ~ D + M1 + M2 + C)
    )
    medsim(data = df, num_sim = 2000, treatment = "D", intv_med = NULL,
           model_spec = spec, seed = rep_id, boot = FALSE)
  }, error = function(e) NULL)

  if (!is.null(med_fit)) {
    scalars <- unlist(med_fit[vapply(med_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))])
    grab <- function(pattern) {
      idx <- grep(pattern, names(scalars), fixed = FALSE)
      if (length(idx) == 0) NA_real_ else scalars[idx[1]]
    }
    results$med_ATE <- grab("TE\\(")
    results$med_DY <- grab("PSE.*D\\s*->\\s*Y\\s*\\)")
    results$med_DM2Y <- grab("PSE.*M2\\s*->\\s*Y")
    results$med_DM1Y <- grab("PSE.*M1\\s*~>\\s*Y")
  } else {
    results$med_ATE <- results$med_DY <- results$med_DM2Y <- results$med_DM1Y <- NA
  }

  # ---- Interventional Effects ----

  # IPW Interventional
  ipw_intv_fit <- tryCatch({
    ipwvent_inner(
      data = df, D = "D", M = "M2", Y = "Y", L = "M1",
      D_formula = as.formula("D ~ C"),
      L_formula = as.formula("M1 ~ D + C"),
      M_formula = as.formula("M2 ~ D + M1 + C"),
      stabilize = TRUE, censor = TRUE
    )
  }, error = function(e) NULL)

  if (!is.null(ipw_intv_fit)) {
    results$ipw_intv_OE <- ipw_intv_fit$OE
    results$ipw_intv_IDE <- ipw_intv_fit$IDE
    results$ipw_intv_IIE <- ipw_intv_fit$IIE
  } else {
    results$ipw_intv_OE <- results$ipw_intv_IDE <- results$ipw_intv_IIE <- NA
  }

  # RWR Interventional
  if (requireNamespace("rwrmed", quietly = TRUE)) {
    rwr_fit <- tryCatch({
      rwrlite(
        data = df, D = "D", C = "C",
        Y_formula = as.formula("Y ~ D + M1 + M2 + C + D:M2"),
        M_formula = as.formula("M2 ~ D + C"),
        L_formula_list = list(as.formula("M1 ~ D + C")),
        boot = FALSE
      )
    }, error = function(e) NULL)

    if (!is.null(rwr_fit)) {
      results$rwr_intv_OE <- rwr_fit$OE
      results$rwr_intv_IDE <- rwr_fit$IDE
      results$rwr_intv_IIE <- rwr_fit$IIE
    } else {
      results$rwr_intv_OE <- results$rwr_intv_IDE <- results$rwr_intv_IIE <- NA
    }
  } else {
    results$rwr_intv_OE <- results$rwr_intv_IDE <- results$rwr_intv_IIE <- NA
  }

  # MedSim Interventional
  medsim_intv_fit <- tryCatch({
    spec_intv <- list(
      list(func = "glm", formula = M1 ~ D + C, args = list(family = binomial())),
      list(func = "glm", formula = M2 ~ D + C, args = list(family = binomial())),
      list(func = "lm", formula = Y ~ D + M1 + M2 + C)
    )
    medsim(data = df, num_sim = 2000, treatment = "D", intv_med = "M2",
           model_spec = spec_intv, seed = rep_id + 1000, boot = FALSE)
  }, error = function(e) NULL)

  if (!is.null(medsim_intv_fit)) {
    scalars <- unlist(medsim_intv_fit[vapply(medsim_intv_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))])
    grab <- function(pattern) {
      idx <- grep(pattern, names(scalars), fixed = FALSE)
      if (length(idx) == 0) NA_real_ else scalars[idx[1]]
    }
    results$med_intv_OE <- grab("OE\\(")
    results$med_intv_IDE <- grab("IDE\\(")
    results$med_intv_IIE <- grab("IIE\\(")
  } else {
    results$med_intv_OE <- results$med_intv_IDE <- results$med_intv_IIE <- NA
  }

  return(as.data.frame(results))
}

# Monte Carlo settings
N_REPS <- 100
N_CORES <- 10

cat(sprintf("Running %d Monte Carlo replications on %d cores...\n", N_REPS, N_CORES))
cat("This may take a few minutes.\n\n")

# Set seed for reproducibility
set.seed(20240101)

# Run MC simulation in parallel
start_time <- Sys.time()
mc_results <- mclapply(1:N_REPS, run_one_rep, mc.cores = N_CORES)
end_time <- Sys.time()

cat(sprintf("Completed in %.1f minutes.\n\n", difftime(end_time, start_time, units = "mins")))

# Combine results
results_df <- do.call(rbind, mc_results)

# Recompute truth with new parameters
cat("Recomputing Monte Carlo truth with modified parameters...\n")
set.seed(99999)
truth_pse <- truth_exp1(N_truth = 1000000L, par = exp1_params)
truth_intv <- truth_intv_exp1(N_truth = 1000000L, par = exp1_params)

cat("\nTrue PSE values:\n")
print(round(truth_pse, 4))
cat("\nTrue Interventional values:\n")
print(round(truth_intv, 4))

# ---- Compute Bias and RMSE for PSE ----
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("MONTE CARLO RESULTS: PATH-SPECIFIC EFFECTS\n")
cat("=" , rep("=", 60), "\n\n")

pse_estimators <- c("lin", "ipw", "reg", "med")
pse_estimands <- list(
  ATE = "ATE",
  DY = "D->Y",
  DM2Y = "D->M2->Y",
  DM1Y = "D->M1~>Y"
)

pse_stats <- data.frame()
for (est in pse_estimators) {
  for (est_key in names(pse_estimands)) {
    col_name <- paste0(est, "_", est_key)
    true_val <- truth_pse[pse_estimands[[est_key]]]

    if (col_name %in% names(results_df)) {
      estimates <- results_df[[col_name]]
      n_valid <- sum(!is.na(estimates))
      mean_est <- mean(estimates, na.rm = TRUE)
      bias <- mean_est - true_val
      rmse <- sqrt(mean((estimates - true_val)^2, na.rm = TRUE))

      pse_stats <- rbind(pse_stats, data.frame(
        Estimator = toupper(est),
        Estimand = pse_estimands[[est_key]],
        Truth = round(true_val, 4),
        Mean = round(mean_est, 4),
        Bias = round(bias, 4),
        RMSE = round(rmse, 4),
        N = n_valid
      ))
    }
  }
}

print(pse_stats, row.names = FALSE)

# Bias summary table
cat("\n--- Bias Summary (PSE) ---\n")
bias_wide <- reshape(pse_stats[, c("Estimator", "Estimand", "Bias")],
                     idvar = "Estimand", timevar = "Estimator",
                     direction = "wide")
names(bias_wide) <- gsub("Bias\\.", "", names(bias_wide))
print(bias_wide, row.names = FALSE)

# RMSE summary table
cat("\n--- RMSE Summary (PSE) ---\n")
rmse_wide <- reshape(pse_stats[, c("Estimator", "Estimand", "RMSE")],
                     idvar = "Estimand", timevar = "Estimator",
                     direction = "wide")
names(rmse_wide) <- gsub("RMSE\\.", "", names(rmse_wide))
print(rmse_wide, row.names = FALSE)

# ---- Compute Bias and RMSE for Interventional Effects ----
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("MONTE CARLO RESULTS: INTERVENTIONAL EFFECTS\n")
cat("=" , rep("=", 60), "\n\n")

intv_estimators <- list(
  ipw_intv = "IPW",
  rwr_intv = "RWR",
  med_intv = "MedSim"
)
intv_estimands <- c("OE", "IDE", "IIE")

intv_stats <- data.frame()
for (est_key in names(intv_estimators)) {
  est_label <- intv_estimators[[est_key]]
  for (estimand in intv_estimands) {
    col_name <- paste0(est_key, "_", estimand)
    true_val <- truth_intv[estimand]

    if (col_name %in% names(results_df)) {
      estimates <- results_df[[col_name]]
      n_valid <- sum(!is.na(estimates))
      mean_est <- mean(estimates, na.rm = TRUE)
      bias <- mean_est - true_val
      rmse <- sqrt(mean((estimates - true_val)^2, na.rm = TRUE))

      intv_stats <- rbind(intv_stats, data.frame(
        Estimator = est_label,
        Estimand = estimand,
        Truth = round(true_val, 4),
        Mean = round(mean_est, 4),
        Bias = round(bias, 4),
        RMSE = round(rmse, 4),
        N = n_valid
      ))
    }
  }
}

print(intv_stats, row.names = FALSE)

# Bias summary table
cat("\n--- Bias Summary (Interventional) ---\n")
bias_wide_intv <- reshape(intv_stats[, c("Estimator", "Estimand", "Bias")],
                          idvar = "Estimand", timevar = "Estimator",
                          direction = "wide")
names(bias_wide_intv) <- gsub("Bias\\.", "", names(bias_wide_intv))
print(bias_wide_intv, row.names = FALSE)

# RMSE summary table
cat("\n--- RMSE Summary (Interventional) ---\n")
rmse_wide_intv <- reshape(intv_stats[, c("Estimator", "Estimand", "RMSE")],
                          idvar = "Estimand", timevar = "Estimator",
                          direction = "wide")
names(rmse_wide_intv) <- gsub("RMSE\\.", "", names(rmse_wide_intv))
print(rmse_wide_intv, row.names = FALSE)

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("MONTE CARLO SIMULATION COMPLETE\n")
cat("=" , rep("=", 60), "\n")
