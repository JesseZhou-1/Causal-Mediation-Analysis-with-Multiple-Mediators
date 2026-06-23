############################################################
## Run Experiment 4: Nonlinear DGP with Mixed Variable Types
## D ~ Poisson, M1 ~ Binary, M2 ~ Normal, Y ~ OrdLogit
## Same types as exp3 but with nonlinear terms & interactions
## All parametric methods biased; medflow should perform best
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
  library(SuperLearner)
  library(ranger)
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

# Ordinal logit simulator (3-level: 0, 1, 2)
rordlogit <- function(n, linear_pred, cuts) {
  u <- runif(n)
  p0 <- invlogit(cuts[1] - linear_pred)
  p1 <- invlogit(cuts[2] - linear_pred)
  ifelse(u < p0, 0L, ifelse(u < p1, 1L, 2L))
}

## ----------------------------
## 2) DGP for Experiment 4
## ----------------------------
exp4_params <- list(
  # Treatment model: D | C ~ Poisson(exp(alpha0 + alpha1*C + alpha2*C^2))
  D_alpha0 = 0.7,
  D_alpha1 = 0.3,
  D_alpha2 = 0.15,   # nonlinear: C^2

  # M1 model: M1 | D, C ~ Bernoulli(logistic(beta0 + beta_D*D + beta_C*C + beta_DC*D*C))
  M1_beta0 = -1.5,
  M1_beta_D = 0.5,
  M1_beta_C = 0.4,
  M1_beta_DC = 0.25,  # interaction: D*C

  # M2 model: M2 | D, M1, C ~ Normal(mu, sigma=1)
  #   mu = gamma_D*D + gamma_M1*M1 + gamma_C*C + gamma_D2*D^2 + gamma_DM1*D*M1
  M2_gamma_D = 0.4,
  M2_gamma_M1 = 0.6,
  M2_gamma_C = 0.3,
  M2_gamma_D2 = 0.15,    # nonlinear: D^2
  M2_gamma_DM1 = -0.1,   # interaction: D*M1
  M2_sigma = 1.0,

  # Y model: Y | D, M1, M2, C ~ OrdLogit(lp_Y, cuts)
  #   lp_Y = delta_D*D + delta_M1*M1 + delta_M2*M2 + delta_C*C
  #        + delta_M1M2*M1*M2 + delta_sinM2*sin(pi*M2/2) + delta_DC*D*C
  Y_delta_D = 0.3,
  Y_delta_M1 = 0.4,
  Y_delta_M2 = 0.5,
  Y_delta_C = 0.3,
  Y_delta_M1M2 = 0.2,     # interaction: M1*M2
  Y_delta_sinM2 = 0.15,   # nonlinear: sin(pi*M2/2) — cannot be captured by polr
  Y_delta_DC = 0.1,        # interaction: D*C
  Y_cuts = c(-0.5, 0.8)
)

gen_exp4 <- function(n, par = exp4_params) {
  C <- rnorm(n)

  # D | C ~ Poisson (with C^2)
  D <- rpois(n, exp(par$D_alpha0 + par$D_alpha1 * C + par$D_alpha2 * C^2))

  # M1 | D, C ~ Bernoulli (with D*C interaction)
  M1 <- rbinom(n, 1, plogis(par$M1_beta0 + par$M1_beta_D * D +
         par$M1_beta_C * C + par$M1_beta_DC * D * C))

  # M2 | D, M1, C ~ Normal (with D^2 and D*M1)
  mu_M2 <- par$M2_gamma_D * D + par$M2_gamma_M1 * M1 + par$M2_gamma_C * C +
            par$M2_gamma_D2 * D^2 + par$M2_gamma_DM1 * D * M1
  M2 <- rnorm(n, mu_M2, par$M2_sigma)

  # Y | D, M1, M2, C ~ OrdLogit (with M1*M2, sin(pi*M2/2), D*C)
  lp_Y <- par$Y_delta_D * D + par$Y_delta_M1 * M1 + par$Y_delta_M2 * M2 +
           par$Y_delta_C * C + par$Y_delta_M1M2 * M1 * M2 +
           par$Y_delta_sinM2 * sin(pi * M2 / 2) + par$Y_delta_DC * D * C
  Y <- rordlogit(n, lp_Y, par$Y_cuts)

  data.frame(C = C, D = D, M1 = M1, M2 = M2, Y = Y)
}

## ----------------------------
## 3) Compute Monte Carlo truth
## ----------------------------

# Helper: M2 mean with nonlinear terms
mu_M2_fn <- function(D_val, M1_val, C_val, par = exp4_params) {
  par$M2_gamma_D * D_val + par$M2_gamma_M1 * M1_val + par$M2_gamma_C * C_val +
  par$M2_gamma_D2 * D_val^2 + par$M2_gamma_DM1 * D_val * M1_val
}

# Helper: Y linear predictor with nonlinear terms
lp_Y_fn <- function(D_val, M1_val, M2_val, C_val, par = exp4_params) {
  par$Y_delta_D * D_val + par$Y_delta_M1 * M1_val + par$Y_delta_M2 * M2_val +
  par$Y_delta_C * C_val + par$Y_delta_M1M2 * M1_val * M2_val +
  par$Y_delta_sinM2 * sin(pi * M2_val / 2) + par$Y_delta_DC * D_val * C_val
}

truth_pse_exp4 <- function(N_truth = 2000000L, par = exp4_params, d = 5, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)

  # Draw M1 under D=d and D=dstar (with D*C interaction)
  M1_d     <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * d +
               par$M1_beta_C * C + par$M1_beta_DC * d * C))
  M1_dstar <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
               par$M1_beta_C * C + par$M1_beta_DC * dstar * C))

  # Draw M2 under various conditions (with D^2 and D*M1)
  M2_d_M1d      <- rnorm(N_truth, mu_M2_fn(d, M1_d, C, par), par$M2_sigma)
  M2_dstar_M1ds <- rnorm(N_truth, mu_M2_fn(dstar, M1_dstar, C, par), par$M2_sigma)
  M2_d_M1ds     <- rnorm(N_truth, mu_M2_fn(d, M1_dstar, C, par), par$M2_sigma)

  # Draw Y under various conditions (with M1*M2, sin, D*C)
  Y_d_M1d_M2d   <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_d_M1d, C, par), par$Y_cuts)
  Y_ds_M1ds_M2ds <- rordlogit(N_truth, lp_Y_fn(dstar, M1_dstar, M2_dstar_M1ds, C, par), par$Y_cuts)
  Y_d_M1ds_M2ds <- rordlogit(N_truth, lp_Y_fn(d, M1_dstar, M2_dstar_M1ds, C, par), par$Y_cuts)
  Y_d_M1ds_M2d  <- rordlogit(N_truth, lp_Y_fn(d, M1_dstar, M2_d_M1ds, C, par), par$Y_cuts)

  c(
    ATE        = mean(Y_d_M1d_M2d) - mean(Y_ds_M1ds_M2ds),
    `D->Y`     = mean(Y_d_M1ds_M2ds) - mean(Y_ds_M1ds_M2ds),
    `D->M2->Y` = mean(Y_d_M1ds_M2d) - mean(Y_d_M1ds_M2ds),
    `D->M1~>Y` = mean(Y_d_M1d_M2d) - mean(Y_d_M1ds_M2d)
  )
}

truth_intv_exp4 <- function(N_truth = 2000000L, par = exp4_params, d = 5, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)

  # M1 under d and dstar (with D*C interaction)
  M1_d     <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * d +
               par$M1_beta_C * C + par$M1_beta_DC * d * C))
  M1_dstar <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
               par$M1_beta_C * C + par$M1_beta_DC * dstar * C))

  # Natural M2
  M2_d_natural     <- rnorm(N_truth, mu_M2_fn(d, M1_d, C, par), par$M2_sigma)
  M2_dstar_natural <- rnorm(N_truth, mu_M2_fn(dstar, M1_dstar, C, par), par$M2_sigma)

  # Natural Y
  Y_d_natural     <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_d_natural, C, par), par$Y_cuts)
  Y_dstar_natural <- rordlogit(N_truth, lp_Y_fn(dstar, M1_dstar, M2_dstar_natural, C, par), par$Y_cuts)

  # M2*(dstar): marginalize M2 over M1|D=dstar (re-draw M1 for independence)
  M1_for_M2star <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
                    par$M1_beta_C * C + par$M1_beta_DC * dstar * C))
  M2_star_dstar <- rnorm(N_truth, mu_M2_fn(dstar, M1_for_M2star, C, par), par$M2_sigma)

  # Y(d, M1(d), M2*(dstar))
  Y_d_M2star <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_star_dstar, C, par), par$Y_cuts)

  c(
    OE  = mean(Y_d_natural) - mean(Y_dstar_natural),
    IDE = mean(Y_d_M2star) - mean(Y_dstar_natural),
    IIE = mean(Y_d_natural) - mean(Y_d_M2star)
  )
}

## ============================================================
## 4) linmed_inner for linpath
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
## 5) Source estimator functions
## ============================================================
source("ipwmed.R")
source("ipwpath.R")
source("linpath.R")
source("medsim.R")
source("rwrlite.R")
source("pathimp.R")
source("dmlpath_custom.R")
library(paths)

cat("All functions loaded.\n\n")

## ============================================================
## 6) Compute truth
## ============================================================
cat("Computing Monte Carlo truth (2M samples)...\n")
truth_pse <- truth_pse_exp4()
truth_intv <- truth_intv_exp4()

cat("True PSEs (d=5 vs dstar=0):\n")
print(round(truth_pse, 4))
cat("\nTrue Interventional Effects:\n")
print(round(truth_intv, 4))
cat("\n")

## ============================================================
## 7) Single-dataset demonstration
## ============================================================
set.seed(2024)
n <- 20000
df <- gen_exp4(n)
cat("Generated dataset with n =", n, "observations.\n")
cat("D distribution:", table(df$D), "\n")
cat("D unique values:", length(unique(df$D)), "\n")
cat("M1 distribution:", table(df$M1), "\n")
cat("M2 range:", round(range(df$M2), 2), "\n")
cat("Y distribution:", table(df$Y), "\n\n")

# Treatment contrast: d=5, dstar=0
d_val <- 5
dstar_val <- 0

## -- Linear estimator --
cat("Running Linear estimator (d=5 vs dstar=0)...\n")
lin_fit <- tryCatch(
  linpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
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

## -- IPW estimator (Poisson auto-detected) --
cat("Running IPW estimator (count treatment)...\n")
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

## -- RegImp estimator (pathimp with simple additive lm models) --
cat("Running RegImp estimator (pathimp, additive lm)...\n")
reg_fit <- tryCatch({
  lm_m0 <- lm(Y ~ D + C, data = df)
  lm_m1 <- lm(Y ~ D + C + M1, data = df)
  lm_m2 <- lm(Y ~ D + C + M1 + M2, data = df)
  pathimp(D = "D", Y = "Y", M = list("M1", "M2"),
          Y_models = list(lm_m0, lm_m1, lm_m2),
          data = df, d = d_val, dstar = dstar_val,
          boot_reps = 2, out_ipw = FALSE)
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(reg_fit)) {
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

## -- MedSim estimator (simple additive parametric specs) --
cat("Running MedSim estimator (glm+lm+polr additive)...\n")
med_fit <- tryCatch({
  df_med <- df
  df_med$M1 <- factor(df_med$M1, ordered = FALSE)
  df_med$Y  <- factor(df_med$Y, ordered = FALSE)

  spec <- list(
    list(func = "glm",  formula = M1 ~ D + C,
         args = list(family = binomial())),
    list(func = "lm",   formula = M2 ~ D + M1 + C),
    list(func = "polr", formula = Y ~ D + M1 + M2 + C,
         args = list(method = "logistic"))
  )
  medsim(data = df_med, num_sim = 5000, cat_list = c("0", "5"),
         treatment = "D", intv_med = NULL, model_spec = spec,
         seed = 999, boot = FALSE)
}, error = function(e) { cat("  Error:", e$message, "\n"); NULL })

if (!is.null(med_fit)) {
  scalars <- unlist(med_fit[vapply(med_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))])
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

sl_lib_rf <- c("SL.ranger", "SL.glm")

## -- DML MR2 estimator (SL.ranger + SL.glm SuperLearner, cross-fitted) --
cat("Running DML MR2 estimator (SL.ranger + SL.glm, 5-fold cross-fitting)...\n")
dml_fit <- tryCatch(
  dmlpath_custom(
    data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
    d = d_val, dstar = dstar_val, n_folds = 5,
    SL_library = sl_lib_rf,
    seed = 2024
  ),
  error = function(e) { cat("  Error:", e$message, "\n"); NULL }
)

if (!is.null(dml_fit)) {
  dml_results <- c(ATE = dml_fit$ATE,
    `D->Y` = dml_fit$PSE[["D->Y"]],
    `D->M2->Y` = dml_fit$PSE[["D->M2->Y"]],
    `D->M1~>Y` = dml_fit$PSE[["D->M1~>Y"]])
  cat("DML MR2:\n"); print(round(dml_results, 4))
  cat("  PSE sum check:", round(sum(dml_fit$PSE) - dml_fit$ATE, 8), "\n")
} else {
  dml_results <- c(ATE = NA, `D->Y` = NA, `D->M2->Y` = NA, `D->M1~>Y` = NA)
}
cat("\n")

## -- IPW-ML estimator (SL.ranger + SL.glm propensities, no regression) --
cat("Running IPW-ML estimator (SL.ranger + SL.glm propensities)...\n")
ipw_sl_fit <- tryCatch(
  dmlpath_custom(
    data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
    d = d_val, dstar = dstar_val, n_folds = 5,
    SL_library = sl_lib_rf,
    seed = 2024, method = "ipw"
  ),
  error = function(e) { cat("  Error:", e$message, "\n"); NULL }
)

if (!is.null(ipw_sl_fit)) {
  ipw_sl_results <- c(ATE = ipw_sl_fit$ATE,
    `D->Y` = ipw_sl_fit$PSE[["D->Y"]],
    `D->M2->Y` = ipw_sl_fit$PSE[["D->M2->Y"]],
    `D->M1~>Y` = ipw_sl_fit$PSE[["D->M1~>Y"]])
  cat("IPW-SL:\n"); print(round(ipw_sl_results, 4))
  cat("  PSE sum check:", round(sum(ipw_sl_fit$PSE) - ipw_sl_fit$ATE, 8), "\n")
} else {
  ipw_sl_results <- c(ATE = NA, `D->Y` = NA, `D->M2->Y` = NA, `D->M1~>Y` = NA)
}
cat("\n")

## -- REG-ML estimator (SL.ranger + SL.glm regressions, no IPW) --
cat("Running REG-ML estimator (SL.ranger + SL.glm regressions)...\n")
reg_sl_fit <- tryCatch(
  dmlpath_custom(
    data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
    d = d_val, dstar = dstar_val, n_folds = 5,
    SL_library = sl_lib_rf,
    seed = 2024, method = "reg"
  ),
  error = function(e) { cat("  Error:", e$message, "\n"); NULL }
)

if (!is.null(reg_sl_fit)) {
  reg_sl_results <- c(ATE = reg_sl_fit$ATE,
    `D->Y` = reg_sl_fit$PSE[["D->Y"]],
    `D->M2->Y` = reg_sl_fit$PSE[["D->M2->Y"]],
    `D->M1~>Y` = reg_sl_fit$PSE[["D->M1~>Y"]])
  cat("REG-SL:\n"); print(round(reg_sl_results, 4))
  cat("  PSE sum check:", round(sum(reg_sl_fit$PSE) - reg_sl_fit$ATE, 8), "\n")
} else {
  reg_sl_results <- c(ATE = NA, `D->Y` = NA, `D->M2->Y` = NA, `D->M1~>Y` = NA)
}
cat("\n")

## -- Summary PSE --
cat("=" , rep("=", 60), "\n", sep = "")
cat("SUMMARY: Path-Specific Effects (d=5 vs dstar=0)\n")
cat("=" , rep("=", 60), "\n")
results_df <- data.frame(
  Estimand = c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y"),
  Truth = round(truth_pse, 4),
  Linear = round(lin_results, 4),
  IPW = round(ipw_results, 4),
  RegImp = round(reg_results, 4),
  MedSim = round(med_results, 4),
  DML = round(dml_results, 4),
  IPW_SL = round(ipw_sl_results, 4),
  REG_SL = round(reg_sl_results, 4)
)
print(results_df, row.names = FALSE)
cat("\nBias:\n")
bias_df <- data.frame(
  Estimand = c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y"),
  Linear = round(lin_results - truth_pse, 4),
  IPW = round(ipw_results - truth_pse, 4),
  RegImp = round(reg_results - truth_pse, 4),
  MedSim = round(med_results - truth_pse, 4),
  DML = round(dml_results - truth_pse, 4),
  IPW_SL = round(ipw_sl_results - truth_pse, 4),
  REG_SL = round(reg_sl_results - truth_pse, 4)
)
print(bias_df, row.names = FALSE)

## ============================================================
## PART 2: Interventional Effects (skip IPW interventional)
## ============================================================
cat("\n\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("INTERVENTIONAL EFFECTS (d=5 vs dstar=0, intervening on M2)\n")
cat("=" , rep("=", 60), "\n\n")

## -- RWR Interventional (simple additive parametric formulas) --
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

## -- MedSim Interventional (simple additive parametric formulas) --
cat("Running MedSim Interventional (glm+lm+polr additive)...\n")
medsim_intv_fit <- tryCatch({
  df_med <- df
  df_med$M1 <- factor(df_med$M1, ordered = FALSE)
  df_med$Y  <- factor(df_med$Y, ordered = FALSE)

  spec_intv <- list(
    list(func = "glm",  formula = M1 ~ D + C,
         args = list(family = binomial())),
    list(func = "lm",   formula = M2 ~ D + C),  # M2 without M1 for interventional
    list(func = "polr", formula = Y ~ D + M1 + M2 + C,
         args = list(method = "logistic"))
  )
  medsim(data = df_med, num_sim = 5000, cat_list = c("0", "5"),
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

## -- Summary Interventional --
cat("=" , rep("=", 60), "\n", sep = "")
cat("SUMMARY: Interventional Effects\n")
cat("=" , rep("=", 60), "\n")
intv_results_df <- data.frame(
  Estimand = c("OE", "IDE", "IIE"),
  Truth = round(truth_intv, 4),
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

run_one_rep_exp4 <- function(rep_id, n = 20000, par = exp4_params, d_val = 5, dstar_val = 0) {
  set.seed(rep_id)

  df <- gen_exp4(n, par)
  results <- list(rep_id = rep_id)

  # ---- PSE Estimators ----

  # 1. Linear
  lin_fit <- tryCatch({
    linpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
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

  # 2. IPW (Poisson auto-detected)
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

  # 3. RegImp (pathimp with simple additive lm models)
  reg_fit <- tryCatch({
    lm_m0 <- lm(Y ~ D + C, data = df)
    lm_m1 <- lm(Y ~ D + C + M1, data = df)
    lm_m2 <- lm(Y ~ D + C + M1 + M2, data = df)
    pathimp(D = "D", Y = "Y", M = list("M1", "M2"),
            Y_models = list(lm_m0, lm_m1, lm_m2),
            data = df, d = d_val, dstar = dstar_val,
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

  # 4. MedSim (PSE, simple additive parametric specs)
  med_fit <- tryCatch({
    df_med <- df
    df_med$M1 <- factor(df_med$M1, ordered = FALSE)
    df_med$Y  <- factor(df_med$Y, ordered = FALSE)
    spec <- list(
      list(func = "glm",  formula = M1 ~ D + C,
           args = list(family = binomial())),
      list(func = "lm",   formula = M2 ~ D + M1 + C),
      list(func = "polr", formula = Y ~ D + M1 + M2 + C,
           args = list(method = "logistic"))
    )
    medsim(data = df_med, num_sim = 5000, cat_list = c("0", "5"),
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

  sl_lib_rf <- c("SL.ranger", "SL.glm")

  # 5. DML MR2 (SL.ranger + SL.glm SuperLearner, 5-fold cross-fitting)
  dml_fit <- tryCatch({
    dmlpath_custom(
      data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
      d = d_val, dstar = dstar_val, n_folds = 5,
      SL_library = sl_lib_rf,
      seed = rep_id + 50000
    )
  }, error = function(e) NULL)

  if (!is.null(dml_fit)) {
    results$dml_ATE <- dml_fit$ATE
    results$dml_DY <- dml_fit$PSE[["D->Y"]]
    results$dml_DM2Y <- dml_fit$PSE[["D->M2->Y"]]
    results$dml_DM1Y <- dml_fit$PSE[["D->M1~>Y"]]
  } else {
    results$dml_ATE <- results$dml_DY <- results$dml_DM2Y <- results$dml_DM1Y <- NA
  }

  # 6. IPW-ML (SL.ranger + SL.glm propensities, no regression)
  ipw_sl_fit <- tryCatch({
    dmlpath_custom(
      data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
      d = d_val, dstar = dstar_val, n_folds = 5,
      SL_library = sl_lib_rf,
      seed = rep_id + 60000, method = "ipw"
    )
  }, error = function(e) NULL)

  if (!is.null(ipw_sl_fit)) {
    results$ipw_sl_ATE <- ipw_sl_fit$ATE
    results$ipw_sl_DY <- ipw_sl_fit$PSE[["D->Y"]]
    results$ipw_sl_DM2Y <- ipw_sl_fit$PSE[["D->M2->Y"]]
    results$ipw_sl_DM1Y <- ipw_sl_fit$PSE[["D->M1~>Y"]]
  } else {
    results$ipw_sl_ATE <- results$ipw_sl_DY <- results$ipw_sl_DM2Y <- results$ipw_sl_DM1Y <- NA
  }

  # 7. REG-ML (SL.ranger + SL.glm regressions, no IPW)
  reg_sl_fit <- tryCatch({
    dmlpath_custom(
      data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
      d = d_val, dstar = dstar_val, n_folds = 5,
      SL_library = sl_lib_rf,
      seed = rep_id + 70000, method = "reg"
    )
  }, error = function(e) NULL)

  if (!is.null(reg_sl_fit)) {
    results$reg_sl_ATE <- reg_sl_fit$ATE
    results$reg_sl_DY <- reg_sl_fit$PSE[["D->Y"]]
    results$reg_sl_DM2Y <- reg_sl_fit$PSE[["D->M2->Y"]]
    results$reg_sl_DM1Y <- reg_sl_fit$PSE[["D->M1~>Y"]]
  } else {
    results$reg_sl_ATE <- results$reg_sl_DY <- results$reg_sl_DM2Y <- results$reg_sl_DM1Y <- NA
  }

  # ---- Interventional Effects (no IPW interventional) ----

  # RWR Interventional (simple additive parametric formulas)
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

  # MedSim Interventional (simple additive parametric formulas)
  medsim_intv_fit <- tryCatch({
    df_med <- df
    df_med$M1 <- factor(df_med$M1, ordered = FALSE)
    df_med$Y  <- factor(df_med$Y, ordered = FALSE)
    spec_intv <- list(
      list(func = "glm",  formula = M1 ~ D + C,
           args = list(family = binomial())),
      list(func = "lm",   formula = M2 ~ D + C),
      list(func = "polr", formula = Y ~ D + M1 + M2 + C,
           args = list(method = "logistic"))
    )
    medsim(data = df_med, num_sim = 5000, cat_list = c("0", "5"),
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

set.seed(20240401)
start_time <- Sys.time()
mc_results <- mclapply(1:N_REPS, run_one_rep_exp4, mc.cores = N_CORES)
end_time <- Sys.time()
cat(sprintf("Completed in %.1f minutes.\n\n", difftime(end_time, start_time, units = "mins")))

# Combine results
results_all <- do.call(rbind, mc_results)

# ---- Compute Bias and RMSE for PSE ----
cat("=" , rep("=", 60), "\n", sep = "")
cat("MONTE CARLO RESULTS: PATH-SPECIFIC EFFECTS\n")
cat("=" , rep("=", 60), "\n\n")

pse_estimators <- c("lin", "ipw", "reg", "med", "dml", "ipw_sl", "reg_sl")
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

intv_estimators <- list(rwr_intv = "RWR", med_intv = "MedSim")
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
