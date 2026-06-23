#!/usr/bin/env Rscript
# ============================================================
# Monte Carlo Worker Script for Experiment 2
# Ordinal Treatment, Mediators, Outcome (all 3 levels: 0,1,2)
# ============================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript mc_exp2_worker.R <task_id> <n_reps> <n_nodes> <n_cores> <output_dir>")
}

task_id    <- as.integer(args[1])
n_reps     <- as.integer(args[2])
n_nodes    <- as.integer(args[3])
n_cores    <- as.integer(args[4])
output_dir <- args[5]

cat(sprintf("Task %d started at %s\n", task_id, Sys.time()))
cat(sprintf("Total replications: %d, Nodes: %d, Cores per node: %d\n", n_reps, n_nodes, n_cores))

# Calculate which replications this task handles
reps_per_node <- ceiling(n_reps / n_nodes)
start_rep <- (task_id - 1) * reps_per_node + 1
end_rep <- min(task_id * reps_per_node, n_reps)
my_reps <- start_rep:end_rep

cat(sprintf("Task %d handling replications %d to %d (%d total)\n",
            task_id, start_rep, end_rep, length(my_reps)))

# ============================================================
# Load packages
# ============================================================
suppressPackageStartupMessages({
  library(stats)
  library(MASS)
  library(nnet)
  library(dplyr)
  library(stringr)
  library(paths)
  library(parallel)
})

# Check for rwrmed
has_rwrmed <- requireNamespace("rwrmed", quietly = TRUE)
if (!has_rwrmed) {
  cat("Warning: rwrmed package not available. RWR estimator will be skipped.\n")
}

# ============================================================
# Helper functions
# ============================================================
demean <- function(x, w = rep(1, length(x))) x - weighted.mean(x, w, na.rm = TRUE)

trimQ <- function(x, low = 0.01, high = 0.99) {
  min_val <- quantile(x, low)
  max_val <- quantile(x, high)
  x[x < min_val] <- min_val
  x[x > max_val] <- max_val
  x
}

invlogit <- function(x) 1 / (1 + exp(-x))

# ============================================================
# DGP: Ordinal (3 levels: 0, 1, 2)
# ============================================================
rordlogit <- function(n, linear_pred, cuts) {
  u <- runif(n)
  p0 <- invlogit(cuts[1] - linear_pred)
  p1 <- invlogit(cuts[2] - linear_pred)
  ifelse(u < p0, 0L, ifelse(u < p1, 1L, 2L))
}

exp2_params <- list(
  D_cuts = c(-0.5, 0.5), D_beta = 0.7,
  M1_cuts = c(-0.8, 0.8), M1_beta_D = 0.9, M1_beta_C = 0.6,
  M2_cuts = c(-0.5, 1.0), M2_beta_D = 0.7, M2_beta_M1 = 1.0, M2_beta_C = 0.5,
  Y_cuts = c(-0.3, 0.9), Y_beta_D = 0.5, Y_beta_M1 = 0.6, Y_beta_M2 = 0.7, Y_beta_C = 0.4
)

gen_exp2 <- function(n, par = exp2_params) {
  C <- rnorm(n)
  D <- rordlogit(n, par$D_beta * C, par$D_cuts)
  M1 <- rordlogit(n, par$M1_beta_D * D + par$M1_beta_C * C, par$M1_cuts)
  M2 <- rordlogit(n, par$M2_beta_D * D + par$M2_beta_M1 * M1 + par$M2_beta_C * C, par$M2_cuts)
  Y <- rordlogit(n, par$Y_beta_D * D + par$Y_beta_M1 * M1 + par$Y_beta_M2 * M2 + par$Y_beta_C * C, par$Y_cuts)
  data.frame(C = C, D = D, M1 = M1, M2 = M2, Y = Y)
}

# ============================================================
# linmed_inner for linpath
# ============================================================
linmed_inner <- function(data, D, M, Y, C = NULL, d = 1, dstar = 0,
                          m = rep(0, length(M)), interaction_DM = FALSE,
                          interaction_DC = FALSE, interaction_MC = FALSE,
                          weights_name = NULL, minimal = FALSE) {
  df <- data; key_vars <- c(D, M, Y, C)
  if (!minimal) {
    miss_summary <- sapply(key_vars, FUN = function(v) c(nmiss = sum(!is.na(df[[v]])), miss = sum(is.na(df[[v]])))) |> t() |> as.data.frame()
  }
  if (is.null(weights_name)) weights <- rep(1, nrow(df)) else weights <- df[[weights_name]]
  for (covariate in C) df[[covariate]] <- demean(df[[covariate]], w = weights)
  m_preds <- paste(c(C, D), collapse = " + ")
  y_preds <- paste(c(C, D, M), collapse = " + ")
  m_forms <- lapply(M, function(x) paste(x, "~", m_preds))
  y_form <- paste(Y, "~", y_preds)
  m_models <- lapply(m_forms, function(x) lm(as.formula(x), data = df, weights = weights))
  names(m_models) <- M
  y_model <- lm(as.formula(y_form), data = df, weights = weights)
  NIE_part <- mapply(function(M_k, M_model_k) M_model_k$coef[[D]] * y_model$coef[[M_k]], M, m_models)
  NDE <- y_model$coef[[D]] * (d - dstar)
  NIE <- sum(NIE_part) * (d - dstar)
  ATE <- NDE + NIE
  if (minimal) list(ATE = ATE, NDE = NDE, NIE = NIE, CDE = NDE)
  else list(ATE = ATE, NDE = NDE, NIE = NIE, CDE = NDE, model_m = m_models, model_y = y_model, miss_summary = miss_summary)
}

# ============================================================
# Source main estimator functions
# ============================================================
# These should be in the same directory (WORK_DIR from sbatch)
source("ipwmed.R")
source("ipwpath.R")
source("ipwvent.R")
source("linpath.R")
source("medsim.R")
source("pathimp.R")

# Source rwrlite if rwrmed is available
if (has_rwrmed) {
  source("rwrlite.R")
}

cat("All functions loaded.\n")

# ============================================================
# Single replication function
# ============================================================
d_val <- 2
dstar_val <- 0

run_one_rep <- function(rep_id, n = 2000, par = exp2_params) {
  # Set seed for reproducibility
  set.seed(rep_id)

  df <- gen_exp2(n, par)
  results <- list(rep_id = rep_id)

  # ---- PSE Estimators ----

  # 1. Linear
  lin_fit <- tryCatch({
    df_lin <- df
    df_lin$D <- as.numeric(df_lin$D); df_lin$M1 <- as.numeric(df_lin$M1)
    df_lin$M2 <- as.numeric(df_lin$M2); df_lin$Y <- as.numeric(df_lin$Y)
    linpath(data = df_lin, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
            d = d_val, dstar = dstar_val, boot = FALSE)
  }, error = function(e) NULL)
  if (!is.null(lin_fit)) {
    results$lin_ATE <- as.numeric(lin_fit$ATE); results$lin_DY <- as.numeric(lin_fit$PSE[["D->Y"]])
    results$lin_DM2Y <- as.numeric(lin_fit$PSE[["D->M2->Y"]]); results$lin_DM1Y <- as.numeric(lin_fit$PSE[["D->M1~>Y"]])
  } else { results$lin_ATE <- results$lin_DY <- results$lin_DM2Y <- results$lin_DM1Y <- NA }

  # 2. IPW (ordinal)
  ipw_fit <- tryCatch({
    ipwpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C", d = d_val, dstar = dstar_val, boot = FALSE)
  }, error = function(e) NULL)
  if (!is.null(ipw_fit)) {
    results$ipw_ATE <- as.numeric(ipw_fit$ATE); results$ipw_DY <- as.numeric(ipw_fit$PSE[["D->Y"]])
    results$ipw_DM2Y <- as.numeric(ipw_fit$PSE[["D->M2->Y"]]); results$ipw_DM1Y <- as.numeric(ipw_fit$PSE[["D->M1~>Y"]])
  } else { results$ipw_ATE <- results$ipw_DY <- results$ipw_DM2Y <- results$ipw_DM1Y <- NA }

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
  } else { results$reg_ATE <- results$reg_DY <- results$reg_DM2Y <- results$reg_DM1Y <- NA }

  # 4. MedSim (polr, PSE)
  med_fit <- tryCatch({
    df_med <- df
    df_med$D <- factor(df_med$D, ordered = FALSE); df_med$M1 <- factor(df_med$M1, ordered = FALSE)
    df_med$M2 <- factor(df_med$M2, ordered = FALSE); df_med$Y <- factor(df_med$Y, ordered = FALSE)
    spec <- list(
      list(func = "polr", formula = M1 ~ D + C, args = list(method = "logistic")),
      list(func = "polr", formula = M2 ~ D + M1 + C, args = list(method = "logistic")),
      list(func = "polr", formula = Y ~ D + M1 + M2 + C, args = list(method = "logistic"))
    )
    medsim(data = df_med, num_sim = 5000, cat_list = c("0", "2"),
           treatment = "D", intv_med = NULL, model_spec = spec,
           seed = rep_id, boot = FALSE)
  }, error = function(e) NULL)
  if (!is.null(med_fit)) {
    scalars <- unlist(med_fit[vapply(med_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))])
    grab <- function(pattern) { idx <- grep(pattern, names(scalars)); if (length(idx) == 0) NA_real_ else scalars[idx[1]] }
    results$med_ATE <- grab("TE\\("); results$med_DY <- grab("PSE.*D\\s*->\\s*Y\\s*\\)")
    results$med_DM2Y <- grab("PSE.*M2\\s*->\\s*Y"); results$med_DM1Y <- grab("PSE.*M1\\s*~>\\s*Y")
  } else { results$med_ATE <- results$med_DY <- results$med_DM2Y <- results$med_DM1Y <- NA }

  # ---- Interventional Effects ----

  # IPW Interventional (ordinal)
  ipw_intv_fit <- tryCatch({
    ipwvent_inner(data = df, D = "D", M = "M2", Y = "Y", L = "M1",
      D_formula = as.formula("D ~ C"), L_formula = as.formula("M1 ~ D + C"),
      M_formula = as.formula("M2 ~ D + M1 + C"), d = d_val, dstar = dstar_val)
  }, error = function(e) NULL)
  if (!is.null(ipw_intv_fit)) {
    results$ipw_intv_OE <- ipw_intv_fit$OE; results$ipw_intv_IDE <- ipw_intv_fit$IDE; results$ipw_intv_IIE <- ipw_intv_fit$IIE
  } else { results$ipw_intv_OE <- results$ipw_intv_IDE <- results$ipw_intv_IIE <- NA }

  # RWR Interventional
  if (has_rwrmed) {
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
      results$rwr_intv_OE <- rwr_intv_fit$OE; results$rwr_intv_IDE <- rwr_intv_fit$IDE; results$rwr_intv_IIE <- rwr_intv_fit$IIE
    } else { results$rwr_intv_OE <- results$rwr_intv_IDE <- results$rwr_intv_IIE <- NA }
  } else { results$rwr_intv_OE <- results$rwr_intv_IDE <- results$rwr_intv_IIE <- NA }

  # MedSim Interventional (polr)
  medsim_intv_fit <- tryCatch({
    df_med <- df
    df_med$D <- factor(df_med$D, ordered = FALSE); df_med$M1 <- factor(df_med$M1, ordered = FALSE)
    df_med$M2 <- factor(df_med$M2, ordered = FALSE); df_med$Y <- factor(df_med$Y, ordered = FALSE)
    spec_intv <- list(
      list(func = "polr", formula = M1 ~ D + C, args = list(method = "logistic")),
      list(func = "polr", formula = M2 ~ D + C, args = list(method = "logistic")),
      list(func = "polr", formula = Y ~ D + M1 + M2 + C, args = list(method = "logistic"))
    )
    medsim(data = df_med, num_sim = 5000, cat_list = c("0", "2"),
           treatment = "D", intv_med = "M2", model_spec = spec_intv,
           seed = rep_id + 10000, boot = FALSE)
  }, error = function(e) NULL)
  if (!is.null(medsim_intv_fit)) {
    scalars <- unlist(medsim_intv_fit[vapply(medsim_intv_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))])
    grab <- function(pattern) { idx <- grep(pattern, names(scalars)); if (length(idx) == 0) NA_real_ else scalars[idx[1]] }
    results$med_intv_OE <- grab("OE\\("); results$med_intv_IDE <- grab("IDE\\("); results$med_intv_IIE <- grab("IIE\\(")
  } else { results$med_intv_OE <- results$med_intv_IDE <- results$med_intv_IIE <- NA }

  return(as.data.frame(results))
}

# ============================================================
# Run replications in parallel
# ============================================================
cat(sprintf("Starting %d replications using %d cores...\n", length(my_reps), n_cores))

results_list <- mclapply(
  my_reps,
  run_one_rep,
  mc.cores = n_cores,
  mc.preschedule = TRUE
)

# Combine results
results_df <- bind_rows(results_list)

# Save results
output_file <- file.path(output_dir, sprintf("results_task_%d.rds", task_id))
saveRDS(results_df, output_file)

cat(sprintf("Task %d completed at %s\n", task_id, Sys.time()))
cat(sprintf("Results saved to: %s\n", output_file))
cat(sprintf("Completed %d replications\n", nrow(results_df)))
