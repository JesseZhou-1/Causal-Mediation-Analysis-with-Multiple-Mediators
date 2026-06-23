############################################################
## Run Experiment 2: Ordinal Treatment, Mediators, Outcome
## All variables ordinal (3 levels: 0, 1, 2)
## Only medsim (with correct polr models) should be unbiased
############################################################

## ----------------------------
## 0) Setup and load packages
## ----------------------------
suppressPackageStartupMessages({
  library(stats)
  library(MASS)
  library(nnet)
  library(dplyr)
  library(stringr)
  library(paths)
})

cat("Packages loaded successfully.\n\n")

## ----------------------------
## 1) Define helper functions
## ----------------------------
demean <- function(x, w = rep(1, length(x))) x - weighted.mean(x, w, na.rm = TRUE)

trimQ <- function(x, low = 0.01, high = 0.99) {
  min <- quantile(x, low)
  max <- quantile(x, high)
  x[x < min] <- min
  x[x > max] <- max
  x
}

invlogit <- function(x) 1 / (1 + exp(-x))

## ----------------------------
## 2) Ordinal DGP (ordered logit)
## ----------------------------
# Generate ordinal variable from ordered logit:
#   latent = linear_pred + logistic_noise
#   Y = 0 if latent < cut1, 1 if cut1 <= latent < cut2, 2 if latent >= cut2
rordlogit <- function(n, linear_pred, cuts) {
  # cuts = c(cut1, cut2) for 3-level ordinal
  # P(Y <= j) = invlogit(cuts[j] - linear_pred)
  u <- runif(n)
  p0 <- invlogit(cuts[1] - linear_pred)
  p1 <- invlogit(cuts[2] - linear_pred)
  ifelse(u < p0, 0L, ifelse(u < p1, 1L, 2L))
}

# DGP parameters for Experiment 2
exp2_params <- list(
  # Treatment model: D | C ~ OrderedLogit
  D_cuts = c(-0.5, 0.5),
  D_beta = 0.7,           # effect of C on D

  # M1 model: M1 | D, C ~ OrderedLogit
  M1_cuts = c(-0.8, 0.8),
  M1_beta_D = 0.9,
  M1_beta_C = 0.6,

  # M2 model: M2 | D, M1, C ~ OrderedLogit
  M2_cuts = c(-0.5, 1.0),
  M2_beta_D = 0.7,
  M2_beta_M1 = 1.0,
  M2_beta_C = 0.5,

  # Y model: Y | D, M1, M2, C ~ OrderedLogit
  Y_cuts = c(-0.3, 0.9),
  Y_beta_D = 0.5,
  Y_beta_M1 = 0.6,
  Y_beta_M2 = 0.7,
  Y_beta_C = 0.4
)

gen_exp2 <- function(n, par = exp2_params) {
  C <- rnorm(n, 0, 1)

  # D | C
  D <- rordlogit(n, par$D_beta * C, par$D_cuts)

  # M1 | D, C
  lp_M1 <- par$M1_beta_D * D + par$M1_beta_C * C
  M1 <- rordlogit(n, lp_M1, par$M1_cuts)

  # M2 | D, M1, C
  lp_M2 <- par$M2_beta_D * D + par$M2_beta_M1 * M1 + par$M2_beta_C * C
  M2 <- rordlogit(n, lp_M2, par$M2_cuts)

  # Y | D, M1, M2, C
  lp_Y <- par$Y_beta_D * D + par$Y_beta_M1 * M1 + par$Y_beta_M2 * M2 + par$Y_beta_C * C
  Y <- rordlogit(n, lp_Y, par$Y_cuts)

  data.frame(
    C = C,
    D = D, M1 = M1, M2 = M2, Y = Y
  )
}

## ----------------------------
## 3) Compute Monte Carlo truth
## ----------------------------
# For ordinal Y, E[Y] is just the mean of the numeric values (0, 1, 2)
# Truth: simulate potential outcomes with large N

truth_pse_exp2 <- function(N_truth = 1000000L, par = exp2_params, d = 2, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)

  # Helper: compute E[Y | d, m1, m2, c] for ordinal Y
  # P(Y=0) = invlogit(cut1 - lp), P(Y=1) = invlogit(cut2 - lp) - P(Y=0), P(Y=2) = 1 - invlogit(cut2 - lp)
  EY_ord <- function(d_val, m1_val, m2_val, c_val) {
    lp <- par$Y_beta_D * d_val + par$Y_beta_M1 * m1_val + par$Y_beta_M2 * m2_val + par$Y_beta_C * c_val
    p0 <- invlogit(par$Y_cuts[1] - lp)
    p01 <- invlogit(par$Y_cuts[2] - lp)
    p1 <- p01 - p0
    p2 <- 1 - p01
    0 * p0 + 1 * p1 + 2 * p2
  }

  # Draw M1 under D=d and D=dstar
  lp_M1_d <- par$M1_beta_D * d + par$M1_beta_C * C
  M1_d <- rordlogit(N_truth, lp_M1_d, par$M1_cuts)
  lp_M1_dstar <- par$M1_beta_D * dstar + par$M1_beta_C * C
  M1_dstar <- rordlogit(N_truth, lp_M1_dstar, par$M1_cuts)

  # Draw M2 under various conditions
  M2_d_M1d <- rordlogit(N_truth,
    par$M2_beta_D * d + par$M2_beta_M1 * M1_d + par$M2_beta_C * C, par$M2_cuts)
  M2_dstar_M1dstar <- rordlogit(N_truth,
    par$M2_beta_D * dstar + par$M2_beta_M1 * M1_dstar + par$M2_beta_C * C, par$M2_cuts)
  M2_d_M1dstar <- rordlogit(N_truth,
    par$M2_beta_D * d + par$M2_beta_M1 * M1_dstar + par$M2_beta_C * C, par$M2_cuts)

  # PSE truth via simulation (draw Y, take mean)
  Y_d_M1d_M2d <- rordlogit(N_truth,
    par$Y_beta_D * d + par$Y_beta_M1 * M1_d + par$Y_beta_M2 * M2_d_M1d + par$Y_beta_C * C, par$Y_cuts)
  Y_dstar_M1dstar_M2dstar <- rordlogit(N_truth,
    par$Y_beta_D * dstar + par$Y_beta_M1 * M1_dstar + par$Y_beta_M2 * M2_dstar_M1dstar + par$Y_beta_C * C, par$Y_cuts)
  Y_d_M1dstar_M2dstar <- rordlogit(N_truth,
    par$Y_beta_D * d + par$Y_beta_M1 * M1_dstar + par$Y_beta_M2 * M2_dstar_M1dstar + par$Y_beta_C * C, par$Y_cuts)
  Y_d_M1dstar_M2d <- rordlogit(N_truth,
    par$Y_beta_D * d + par$Y_beta_M1 * M1_dstar + par$Y_beta_M2 * M2_d_M1dstar + par$Y_beta_C * C, par$Y_cuts)

  c(
    ATE       = mean(Y_d_M1d_M2d) - mean(Y_dstar_M1dstar_M2dstar),
    `D->Y`    = mean(Y_d_M1dstar_M2dstar) - mean(Y_dstar_M1dstar_M2dstar),
    `D->M2->Y` = mean(Y_d_M1dstar_M2d) - mean(Y_d_M1dstar_M2dstar),
    `D->M1~>Y` = mean(Y_d_M1d_M2d) - mean(Y_d_M1dstar_M2d)
  )
}

truth_intv_exp2 <- function(N_truth = 1000000L, par = exp2_params, d = 2, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)

  # M1 under d and dstar
  M1_d <- rordlogit(N_truth, par$M1_beta_D * d + par$M1_beta_C * C, par$M1_cuts)
  M1_dstar <- rordlogit(N_truth, par$M1_beta_D * dstar + par$M1_beta_C * C, par$M1_cuts)

  # Natural M2
  M2_d_natural <- rordlogit(N_truth,
    par$M2_beta_D * d + par$M2_beta_M1 * M1_d + par$M2_beta_C * C, par$M2_cuts)
  M2_dstar_natural <- rordlogit(N_truth,
    par$M2_beta_D * dstar + par$M2_beta_M1 * M1_dstar + par$M2_beta_C * C, par$M2_cuts)

  # Natural Y
  Y_d_natural <- rordlogit(N_truth,
    par$Y_beta_D * d + par$Y_beta_M1 * M1_d + par$Y_beta_M2 * M2_d_natural + par$Y_beta_C * C, par$Y_cuts)
  Y_dstar_natural <- rordlogit(N_truth,
    par$Y_beta_D * dstar + par$Y_beta_M1 * M1_dstar + par$Y_beta_M2 * M2_dstar_natural + par$Y_beta_C * C, par$Y_cuts)

  # M2*(dstar): marginalize M2 over M1|D=dstar
  M1_for_M2star <- rordlogit(N_truth, par$M1_beta_D * dstar + par$M1_beta_C * C, par$M1_cuts)
  M2_star_dstar <- rordlogit(N_truth,
    par$M2_beta_D * dstar + par$M2_beta_M1 * M1_for_M2star + par$M2_beta_C * C, par$M2_cuts)

  # Y(d, M1(d), M2*(dstar))
  Y_d_M2star <- rordlogit(N_truth,
    par$Y_beta_D * d + par$Y_beta_M1 * M1_d + par$Y_beta_M2 * M2_star_dstar + par$Y_beta_C * C, par$Y_cuts)

  c(
    OE  = mean(Y_d_natural) - mean(Y_dstar_natural),
    IDE = mean(Y_d_M2star) - mean(Y_dstar_natural),
    IIE = mean(Y_d_natural) - mean(Y_d_M2star)
  )
}

## ============================================================
## 4) Note: Regression Imputation uses pathimp (paths package)
## ============================================================
# pathimp() wraps paths::paths() with d/dstar support.
# For exp2, we pass misspecified lm models (treating ordinal Y as numeric).

## ============================================================
## 5) linmed_inner for linpath (same as exp1)
## ============================================================
linmed_inner <- function(
    data, D, M, Y, C = NULL, d = 1, dstar = 0,
    m = rep(0, length(M)),
    interaction_DM = FALSE, interaction_DC = FALSE, interaction_MC = FALSE,
    weights_name = NULL, minimal = FALSE
) {
  df <- data
  key_vars <- c(D, M, Y, C)
  if (!minimal) {
    miss_summary <- sapply(key_vars, FUN = function(v) c(nmiss = sum(!is.na(df[[v]])), miss = sum(is.na(df[[v]])))) |> t() |> as.data.frame()
  }
  if (is.null(weights_name)) { weights <- rep(1, nrow(df)) } else { weights <- df[[weights_name]] }
  for (covariate in C) df[[covariate]] <- demean(df[[covariate]], w = weights)
  m_preds <- paste(c(C, D), collapse = " + ")
  if (interaction_DC) m_preds <- paste(m_preds, "+", paste(D, C, sep = ":", collapse = " + "))
  y_preds <- paste(c(C, D, M), collapse = " + ")
  if (interaction_DM) y_preds <- paste(y_preds, "+", paste(D, M, sep = ":", collapse = " + "))
  if (interaction_DC) y_preds <- paste(y_preds, "+", paste(D, C, sep = ":", collapse = " + "))
  if (interaction_MC) y_preds <- paste(y_preds, "+", paste(outer(M, C, FUN = "paste", sep = ":"), collapse = " + "))
  m_forms <- lapply(M, function(x) paste(x, "~", m_preds))
  y_form <- paste(Y, "~", y_preds)
  m_models <- lapply(m_forms, function(x) lm(as.formula(x), data = df, weights = weights))
  names(m_models) <- M
  y_model <- lm(as.formula(y_form), data = df, weights = weights)

  if (interaction_DM) {
    NDE_part <- mapply(function(M_k, M_model_k) y_model$coef[[paste0(D, ":", M_k)]] * (M_model_k$coef[["(Intercept)"]] + M_model_k$coef[[D]] * dstar), M, m_models)
    NIE_part <- mapply(function(M_k, M_model_k) M_model_k$coef[[D]] * (y_model$coef[[M_k]] + y_model$coef[[paste0(D, ":", M_k)]] * d), M, m_models)
    CDE_part <- mapply(function(M_k, m_k) y_model$coef[[paste0(D, ":", M_k)]] * m_k, M, m)
  } else {
    NDE_part <- 0
    NIE_part <- mapply(function(M_k, M_model_k) M_model_k$coef[[D]] * y_model$coef[[M_k]], M, m_models)
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

## ============================================================
## 6) Source estimator functions (ipwmed, ipwpath, ipwvent, etc.)
## ============================================================
source("ipwmed.R")
source("ipwpath.R")
source("ipwvent.R")
source("linpath.R")
source("medsim.R")
source("rwrlite.R")
source("pathimp.R")
library(paths)

cat("All functions loaded.\n\n")

## ============================================================
## 7) Compute truth
## ============================================================
cat("Computing Monte Carlo truth (1M samples)...\n")
truth_pse <- truth_pse_exp2()
truth_intv <- truth_intv_exp2()

cat("True PSEs (d=2 vs dstar=0):\n")
print(round(truth_pse, 4))
cat("\nTrue Interventional Effects:\n")
print(round(truth_intv, 4))
cat("\n")

## ============================================================
## 8) Single-dataset demonstration
## ============================================================
set.seed(2024)
n <- 1000
df <- gen_exp2(n)
cat("Generated dataset with n =", n, "observations.\n")
cat("D distribution:", table(df$D), "\n")
cat("M1 distribution:", table(df$M1), "\n")
cat("M2 distribution:", table(df$M2), "\n")
cat("Y distribution:", table(df$Y), "\n\n")

# Treatment contrast: d=2, dstar=0
d_val <- 2
dstar_val <- 0

## -- Linear estimator --
cat("Running Linear estimator (d=2 vs dstar=0)...\n")
# linpath needs numeric D, M, Y
df_lin <- df
df_lin$D <- as.numeric(df_lin$D)
df_lin$M1 <- as.numeric(df_lin$M1)
df_lin$M2 <- as.numeric(df_lin$M2)
df_lin$Y <- as.numeric(df_lin$Y)

lin_fit <- tryCatch(
  linpath(data = df_lin, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
          d = d_val, dstar = dstar_val, boot = FALSE),
  error = function(e) { cat("  Error:", e$message, "\n"); NULL }
)

if (!is.null(lin_fit)) {
  lin_results <- c(ATE = as.numeric(lin_fit$ATE),
    `D->Y` = as.numeric(lin_fit$PSE[["D->Y"]]),
    `D->M2->Y` = as.numeric(lin_fit$PSE[["D->M2->Y"]]),
    `D->M1~>Y` = as.numeric(lin_fit$PSE[["D->M1~>Y"]]))
  cat("Linear:\n"); print(round(lin_results, 4))
} else {
  lin_results <- c(ATE = NA, `D->Y` = NA, `D->M2->Y` = NA, `D->M1~>Y` = NA)
}
cat("\n")

## -- IPW estimator (ordinal) --
cat("Running IPW estimator (ordinal)...\n")
ipw_fit <- tryCatch(
  ipwpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
          d = d_val, dstar = dstar_val, boot = FALSE),
  error = function(e) { cat("  Error:", e$message, "\n"); NULL }
)

if (!is.null(ipw_fit)) {
  ipw_results <- c(ATE = as.numeric(ipw_fit$ATE),
    `D->Y` = as.numeric(ipw_fit$PSE[["D->Y"]]),
    `D->M2->Y` = as.numeric(ipw_fit$PSE[["D->M2->Y"]]),
    `D->M1~>Y` = as.numeric(ipw_fit$PSE[["D->M1~>Y"]]))
  cat("IPW:\n"); print(round(ipw_results, 4))
} else {
  ipw_results <- c(ATE = NA, `D->Y` = NA, `D->M2->Y` = NA, `D->M1~>Y` = NA)
}
cat("\n")

## -- RegImp estimator (pathimp with misspecified lm) --
cat("Running RegImp estimator (pathimp, lm misspecified)...\n")
reg_fit <- tryCatch({
  # Fit K+1 nested lm models (intentionally misspecified for ordinal Y)
  df_reg <- df
  df_reg$D <- as.numeric(df_reg$D)
  df_reg$M1 <- as.numeric(df_reg$M1)
  df_reg$M2 <- as.numeric(df_reg$M2)
  df_reg$Y <- as.numeric(df_reg$Y)
  lm_m0 <- lm(Y ~ D + C, data = df_reg)
  lm_m1 <- lm(Y ~ D + C + M1, data = df_reg)
  lm_m2 <- lm(Y ~ D + C + M1 + M2, data = df_reg)
  pathimp(D = "D", Y = "Y", M = list("M1", "M2"),
          Y_models = list(lm_m0, lm_m1, lm_m2),
          data = df_reg, d = d_val, dstar = dstar_val,
          boot_reps = 2, out_ipw = FALSE)
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(reg_fit)) {
  # Extract point estimates from paths object (Type I pure imputation)
  pure_df <- reg_fit$org_obj$pure
  pure_t1 <- pure_df[pure_df$decomposition == "Type I", ]
  reg_results <- c(
    ATE = pure_t1$estimate[pure_t1$estimand == "total"],
    `D->Y` = pure_t1$estimate[pure_t1$estimand == "direct"],
    `D->M2->Y` = pure_t1$estimate[pure_t1$estimand == "via M2"],
    `D->M1~>Y` = pure_t1$estimate[pure_t1$estimand == "via M1"])
  cat("RegImp:\n"); print(round(reg_results, 4))
} else {
  reg_results <- c(ATE = NA, `D->Y` = NA, `D->M2->Y` = NA, `D->M1~>Y` = NA)
}
cat("\n")

## -- MedSim estimator (correctly specified polr) --
cat("Running MedSim estimator (polr models)...\n")
med_fit <- tryCatch({
  # Convert to factors for polr
  df_med <- df
  df_med$D <- factor(df_med$D, ordered = FALSE)
  df_med$M1 <- factor(df_med$M1, ordered = FALSE)
  df_med$M2 <- factor(df_med$M2, ordered = FALSE)
  df_med$Y <- factor(df_med$Y, ordered = FALSE)

  spec <- list(
    list(func = "polr", formula = M1 ~ D + C, args = list(method = "logistic")),
    list(func = "polr", formula = M2 ~ D + M1 + C, args = list(method = "logistic")),
    list(func = "polr", formula = Y ~ D + M1 + M2 + C, args = list(method = "logistic"))
  )
  medsim(data = df_med, num_sim = 2000, cat_list = c("0", "2"),
         treatment = "D", intv_med = NULL, model_spec = spec,
         seed = 999, boot = FALSE)
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(med_fit)) {
  nm <- names(med_fit)
  is_scalar <- vapply(med_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))
  scalars <- unlist(med_fit[is_scalar])
  grab <- function(pattern) {
    idx <- grep(pattern, names(scalars), fixed = FALSE)
    if (length(idx) == 0) return(NA_real_)
    scalars[idx[1]]
  }
  med_results <- c(ATE = grab("TE\\("),
    `D->Y` = grab("PSE.*D\\s*->\\s*Y\\s*\\)"),
    `D->M2->Y` = grab("PSE.*M2\\s*->\\s*Y"),
    `D->M1~>Y` = grab("PSE.*M1\\s*~>\\s*Y"))
  cat("MedSim:\n"); print(round(med_results, 4))
} else {
  med_results <- c(ATE = NA, `D->Y` = NA, `D->M2->Y` = NA, `D->M1~>Y` = NA)
}
cat("\n")

## -- Summary PSE --
cat("=" , rep("=", 60), "\n", sep = "")
cat("SUMMARY: Path-Specific Effects (d=2 vs dstar=0)\n")
cat("=" , rep("=", 60), "\n")
results_df <- data.frame(
  Estimand = c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y"),
  Truth = round(truth_pse, 4),
  Linear = round(lin_results, 4),
  IPW = round(ipw_results, 4),
  RegImp = round(reg_results, 4),
  MedSim = round(med_results, 4)
)
print(results_df, row.names = FALSE)
cat("\nBias:\n")
bias_df <- data.frame(
  Estimand = c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y"),
  Linear = round(lin_results - truth_pse, 4),
  IPW = round(ipw_results - truth_pse, 4),
  RegImp = round(reg_results - truth_pse, 4),
  MedSim = round(med_results - truth_pse, 4)
)
print(bias_df, row.names = FALSE)

## ============================================================
## PART 2: Interventional Effects
## ============================================================
cat("\n\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("INTERVENTIONAL EFFECTS (d=2 vs dstar=0, intervening on M2)\n")
cat("=" , rep("=", 60), "\n\n")

## -- IPW Interventional (ordinal) --
cat("Running IPW Interventional (ordinal)...\n")
ipw_intv_fit <- tryCatch({
  ipwvent_inner(
    data = df, D = "D", M = "M2", Y = "Y", L = "M1",
    D_formula = as.formula("D ~ C"),
    L_formula = as.formula("M1 ~ D + C"),
    M_formula = as.formula("M2 ~ D + M1 + C"),
    d = d_val, dstar = dstar_val
  )
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(ipw_intv_fit)) {
  ipw_intv_results <- c(OE = ipw_intv_fit$OE, IDE = ipw_intv_fit$IDE, IIE = ipw_intv_fit$IIE)
  cat("IPW Interventional:\n"); print(round(ipw_intv_results, 4))
} else {
  ipw_intv_results <- c(OE = NA, IDE = NA, IIE = NA)
}
cat("\n")

## -- MedSim Interventional --
cat("Running MedSim Interventional (polr)...\n")
medsim_intv_fit <- tryCatch({
  df_med <- df
  df_med$D <- factor(df_med$D, ordered = FALSE)
  df_med$M1 <- factor(df_med$M1, ordered = FALSE)
  df_med$M2 <- factor(df_med$M2, ordered = FALSE)
  df_med$Y <- factor(df_med$Y, ordered = FALSE)

  spec_intv <- list(
    list(func = "polr", formula = M1 ~ D + C, args = list(method = "logistic")),
    list(func = "polr", formula = M2 ~ D + C, args = list(method = "logistic")),
    list(func = "polr", formula = Y ~ D + M1 + M2 + C, args = list(method = "logistic"))
  )
  medsim(data = df_med, num_sim = 2000, cat_list = c("0", "2"),
         treatment = "D", intv_med = "M2", model_spec = spec_intv,
         seed = 999, boot = FALSE)
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(medsim_intv_fit)) {
  scalars <- unlist(medsim_intv_fit[vapply(medsim_intv_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))])
  grab <- function(pattern) {
    idx <- grep(pattern, names(scalars), fixed = FALSE)
    if (length(idx) == 0) return(NA_real_)
    scalars[idx[1]]
  }
  medsim_intv_results <- c(OE = grab("OE\\("), IDE = grab("IDE\\("), IIE = grab("IIE\\("))
  cat("MedSim Interventional:\n"); print(round(medsim_intv_results, 4))
} else {
  medsim_intv_results <- c(OE = NA, IDE = NA, IIE = NA)
}
cat("\n")

## -- RWR Interventional --
cat("Running RWR Interventional...\n")
rwr_intv_fit <- tryCatch({
  rwrlite(
    data = df, D = "D", C = "C",
    d = d_val, dstar = dstar_val,
    Y_formula = as.formula("Y ~ D + M2 + M1 + C"),
    M_formula = as.formula("M2 ~ D + C"),
    L_formula_list = list(as.formula("M1 ~ D + C")),
    boot = FALSE
  )
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(rwr_intv_fit)) {
  rwr_intv_results <- c(OE = rwr_intv_fit$OE, IDE = rwr_intv_fit$IDE, IIE = rwr_intv_fit$IIE)
  cat("RWR Interventional:\n"); print(round(rwr_intv_results, 4))
} else {
  rwr_intv_results <- c(OE = NA, IDE = NA, IIE = NA)
}
cat("\n")

## -- Summary Interventional --
cat("=" , rep("=", 60), "\n", sep = "")
cat("SUMMARY: Interventional Effects\n")
cat("=" , rep("=", 60), "\n")
intv_results_df <- data.frame(
  Estimand = c("OE", "IDE", "IIE"),
  Truth = round(truth_intv, 4),
  IPW = round(ipw_intv_results, 4),
  RWR = round(rwr_intv_results, 4),
  MedSim = round(medsim_intv_results, 4)
)
print(intv_results_df, row.names = FALSE)
cat("\n")

##############################################################
## PART 3: Monte Carlo Simulation (100 reps, 10 cores)
##############################################################

cat("\n\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("PART 3: MONTE CARLO SIMULATION (100 replications, 10 cores)\n")
cat("=" , rep("=", 60), "\n\n")

library(parallel)

run_one_rep_exp2 <- function(rep_id, n = 1000, par = exp2_params, d_val = 2, dstar_val = 0) {
  # Set seed for reproducibility
  set.seed(rep_id)

  df <- gen_exp2(n, par)
  results <- list(rep_id = rep_id)

  # ---- PSE Estimators ----

  # 1. Linear
  lin_fit <- tryCatch({
    df_lin <- df
    df_lin$D <- as.numeric(df_lin$D)
    df_lin$M1 <- as.numeric(df_lin$M1)
    df_lin$M2 <- as.numeric(df_lin$M2)
    df_lin$Y <- as.numeric(df_lin$Y)
    linpath(data = df_lin, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
            d = d_val, dstar = dstar_val, boot = FALSE)
  }, error = function(e) NULL)

  if (!is.null(lin_fit)) {
    results$lin_ATE <- as.numeric(lin_fit$ATE)
    results$lin_DY <- as.numeric(lin_fit$PSE[["D->Y"]])
    results$lin_DM2Y <- as.numeric(lin_fit$PSE[["D->M2->Y"]])
    results$lin_DM1Y <- as.numeric(lin_fit$PSE[["D->M1~>Y"]])
  } else {
    results$lin_ATE <- results$lin_DY <- results$lin_DM2Y <- results$lin_DM1Y <- NA
  }

  # 2. IPW (ordinal)
  ipw_fit <- tryCatch({
    ipwpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
            d = d_val, dstar = dstar_val, boot = FALSE)
  }, error = function(e) NULL)

  if (!is.null(ipw_fit)) {
    results$ipw_ATE <- as.numeric(ipw_fit$ATE)
    results$ipw_DY <- as.numeric(ipw_fit$PSE[["D->Y"]])
    results$ipw_DM2Y <- as.numeric(ipw_fit$PSE[["D->M2->Y"]])
    results$ipw_DM1Y <- as.numeric(ipw_fit$PSE[["D->M1~>Y"]])
  } else {
    results$ipw_ATE <- results$ipw_DY <- results$ipw_DM2Y <- results$ipw_DM1Y <- NA
  }

  # 3. RegImp (pathimp with misspecified lm)
  reg_fit <- tryCatch({
    df_reg <- df
    df_reg$D <- as.numeric(df_reg$D); df_reg$M1 <- as.numeric(df_reg$M1)
    df_reg$M2 <- as.numeric(df_reg$M2); df_reg$Y <- as.numeric(df_reg$Y)
    lm_m0 <- lm(Y ~ D + C, data = df_reg)
    lm_m1 <- lm(Y ~ D + C + M1, data = df_reg)
    lm_m2 <- lm(Y ~ D + C + M1 + M2, data = df_reg)
    pathimp(D = "D", Y = "Y", M = list("M1", "M2"),
            Y_models = list(lm_m0, lm_m1, lm_m2),
            data = df_reg, d = d_val, dstar = dstar_val,
            boot_reps = 2, out_ipw = FALSE)
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

  # 4. MedSim (PSE, polr)
  med_fit <- tryCatch({
    df_med <- df
    df_med$D <- factor(df_med$D, ordered = FALSE)
    df_med$M1 <- factor(df_med$M1, ordered = FALSE)
    df_med$M2 <- factor(df_med$M2, ordered = FALSE)
    df_med$Y <- factor(df_med$Y, ordered = FALSE)
    spec <- list(
      list(func = "polr", formula = M1 ~ D + C, args = list(method = "logistic")),
      list(func = "polr", formula = M2 ~ D + M1 + C, args = list(method = "logistic")),
      list(func = "polr", formula = Y ~ D + M1 + M2 + C, args = list(method = "logistic"))
    )
    medsim(data = df_med, num_sim = 2000, cat_list = c("0", "2"),
           treatment = "D", intv_med = NULL, model_spec = spec,
           seed = rep_id, boot = FALSE)
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

  # IPW Interventional (ordinal)
  ipw_intv_fit <- tryCatch({
    ipwvent_inner(
      data = df, D = "D", M = "M2", Y = "Y", L = "M1",
      D_formula = as.formula("D ~ C"),
      L_formula = as.formula("M1 ~ D + C"),
      M_formula = as.formula("M2 ~ D + M1 + C"),
      d = d_val, dstar = dstar_val
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
  rwr_intv_fit <- tryCatch({
    rwrlite(
      data = df, D = "D", C = "C",
      d = d_val, dstar = dstar_val,
      Y_formula = as.formula("Y ~ D + M2 + M1 + C"),
      M_formula = as.formula("M2 ~ D + C"),
      L_formula_list = list(as.formula("M1 ~ D + C")),
      boot = FALSE
    )
  }, error = function(e) NULL)

  if (!is.null(rwr_intv_fit)) {
    results$rwr_intv_OE <- rwr_intv_fit$OE
    results$rwr_intv_IDE <- rwr_intv_fit$IDE
    results$rwr_intv_IIE <- rwr_intv_fit$IIE
  } else {
    results$rwr_intv_OE <- results$rwr_intv_IDE <- results$rwr_intv_IIE <- NA
  }

  # MedSim Interventional (polr)
  medsim_intv_fit <- tryCatch({
    df_med <- df
    df_med$D <- factor(df_med$D, ordered = FALSE)
    df_med$M1 <- factor(df_med$M1, ordered = FALSE)
    df_med$M2 <- factor(df_med$M2, ordered = FALSE)
    df_med$Y <- factor(df_med$Y, ordered = FALSE)
    spec_intv <- list(
      list(func = "polr", formula = M1 ~ D + C, args = list(method = "logistic")),
      list(func = "polr", formula = M2 ~ D + C, args = list(method = "logistic")),
      list(func = "polr", formula = Y ~ D + M1 + M2 + C, args = list(method = "logistic"))
    )
    medsim(data = df_med, num_sim = 2000, cat_list = c("0", "2"),
           treatment = "D", intv_med = "M2", model_spec = spec_intv,
           seed = rep_id + 1000, boot = FALSE)
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

# MC settings
N_REPS <- 100
N_CORES <- 10

cat(sprintf("Running %d Monte Carlo replications on %d cores...\n", N_REPS, N_CORES))
cat("This may take several minutes.\n\n")

set.seed(20240201)
start_time <- Sys.time()
mc_results <- mclapply(1:N_REPS, run_one_rep_exp2, mc.cores = N_CORES)
end_time <- Sys.time()
cat(sprintf("Completed in %.1f minutes.\n\n", difftime(end_time, start_time, units = "mins")))

# Combine results
results_all <- do.call(rbind, mc_results)

# ---- Compute Bias and RMSE for PSE ----
cat("=" , rep("=", 60), "\n", sep = "")
cat("MONTE CARLO RESULTS: PATH-SPECIFIC EFFECTS\n")
cat("=" , rep("=", 60), "\n\n")

pse_estimators <- c("lin", "ipw", "reg", "med")
pse_estimands <- list(ATE = "ATE", DY = "D->Y", DM2Y = "D->M2->Y", DM1Y = "D->M1~>Y")

pse_stats <- data.frame()
for (est in pse_estimators) {
  for (est_key in names(pse_estimands)) {
    col_name <- paste0(est, "_", est_key)
    true_val <- truth_pse[pse_estimands[[est_key]]]
    if (col_name %in% names(results_all)) {
      estimates <- results_all[[col_name]]
      n_valid <- sum(!is.na(estimates))
      mean_est <- mean(estimates, na.rm = TRUE)
      bias <- mean_est - true_val
      rmse <- sqrt(mean((estimates - true_val)^2, na.rm = TRUE))
      pse_stats <- rbind(pse_stats, data.frame(
        Estimator = toupper(est), Estimand = pse_estimands[[est_key]],
        Truth = round(true_val, 4), Mean = round(mean_est, 4),
        Bias = round(bias, 4), RMSE = round(rmse, 4), N = n_valid))
    }
  }
}
print(pse_stats, row.names = FALSE)

cat("\n--- Bias Summary (PSE) ---\n")
bias_wide <- reshape(pse_stats[, c("Estimator", "Estimand", "Bias")],
                     idvar = "Estimand", timevar = "Estimator", direction = "wide")
names(bias_wide) <- gsub("Bias\\.", "", names(bias_wide))
print(bias_wide, row.names = FALSE)

cat("\n--- RMSE Summary (PSE) ---\n")
rmse_wide <- reshape(pse_stats[, c("Estimator", "Estimand", "RMSE")],
                     idvar = "Estimand", timevar = "Estimator", direction = "wide")
names(rmse_wide) <- gsub("RMSE\\.", "", names(rmse_wide))
print(rmse_wide, row.names = FALSE)

# ---- Compute Bias and RMSE for Interventional Effects ----
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("MONTE CARLO RESULTS: INTERVENTIONAL EFFECTS\n")
cat("=" , rep("=", 60), "\n\n")

intv_estimators <- list(ipw_intv = "IPW", rwr_intv = "RWR", med_intv = "MedSim")
intv_estimands <- c("OE", "IDE", "IIE")

intv_stats <- data.frame()
for (est_key in names(intv_estimators)) {
  est_label <- intv_estimators[[est_key]]
  for (estimand in intv_estimands) {
    col_name <- paste0(est_key, "_", estimand)
    true_val <- truth_intv[estimand]
    if (col_name %in% names(results_all)) {
      estimates <- results_all[[col_name]]
      n_valid <- sum(!is.na(estimates))
      mean_est <- mean(estimates, na.rm = TRUE)
      bias <- mean_est - true_val
      rmse <- sqrt(mean((estimates - true_val)^2, na.rm = TRUE))
      intv_stats <- rbind(intv_stats, data.frame(
        Estimator = est_label, Estimand = estimand,
        Truth = round(true_val, 4), Mean = round(mean_est, 4),
        Bias = round(bias, 4), RMSE = round(rmse, 4), N = n_valid))
    }
  }
}
print(intv_stats, row.names = FALSE)

cat("\n--- Bias Summary (Interventional) ---\n")
bias_wide_intv <- reshape(intv_stats[, c("Estimator", "Estimand", "Bias")],
                          idvar = "Estimand", timevar = "Estimator", direction = "wide")
names(bias_wide_intv) <- gsub("Bias\\.", "", names(bias_wide_intv))
print(bias_wide_intv, row.names = FALSE)

cat("\n--- RMSE Summary (Interventional) ---\n")
rmse_wide_intv <- reshape(intv_stats[, c("Estimator", "Estimand", "RMSE")],
                          idvar = "Estimand", timevar = "Estimator", direction = "wide")
names(rmse_wide_intv) <- gsub("RMSE\\.", "", names(rmse_wide_intv))
print(rmse_wide_intv, row.names = FALSE)

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("MONTE CARLO SIMULATION COMPLETE\n")
cat("=" , rep("=", 60), "\n")
