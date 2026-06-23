#!/usr/bin/env Rscript
# ============================================================
# Monte Carlo Worker Script for Experiment 1
# Estimates Path-Specific Effects and Interventional Effects
# ============================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript mc_exp1_worker.R <task_id> <n_reps> <n_nodes> <n_cores> <output_dir>")
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
  library(tidyr)
  library(paths)
  library(stringr)
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

comb_list_vec <- function(...) {
  mapply(c, ..., SIMPLIFY = FALSE)
}

invlogit <- function(x) 1 / (1 + exp(-x))

# ============================================================
# DGP Parameters
# Modified: shift intercepts to create more extreme probabilities
# where logistic curves are more nonlinear (far from 0.5)
# b10 = -1.5: baseline P(M1=1|D=0,C=0) ≈ 0.18
# b20 = 1.2: baseline P(M2=1|D=0,M1=0,C=0) ≈ 0.77
# ============================================================
exp1_params <- list(
  a0 = -0.2, a1 = 0.7,                           # Treatment model (unchanged)
  b10 = -1.5, b11 = 1.2, b12 = 0.8,              # M1: low baseline, strong D effect
  b20 = 1.2, b21 = 0.8, b22 = 1.5, b23 = 0.6,   # M2: high baseline, strong M1 effect
  g0 = 0.0, g1 = 0.5, g2 = 0.6, g3 = 0.7, g4 = 0.4,  # Outcome model (unchanged)
  sigma = 1.0
)

# ============================================================
# Data generation function
# ============================================================
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

# ============================================================
# linmed_inner function (for linpath)
# ============================================================
linmed_inner <- function(
    data, D, M, Y, C = NULL, d = 1, dstar = 0,
    m = rep(0, length(M)), interaction_DM = FALSE,
    interaction_DC = FALSE, interaction_MC = FALSE,
    weights_name = NULL, minimal = FALSE
) {
  df <- data
  key_vars <- c(D, M, Y, C)
  if (!minimal) {
    miss_summary <- sapply(
      key_vars,
      FUN = function(v) c(nmiss = sum(!is.na(df[[v]])), miss = sum(is.na(df[[v]])))
    ) |> t() |> as.data.frame()
  }

  if (is.null(weights_name)) {
    weights <- rep(1, nrow(df))
  } else {
    weights <- df[[weights_name]]
  }

  for(covariate in C) df[[covariate]] <- demean(df[[covariate]], w = weights)

  m_preds <- paste(c(C, D), collapse = " + ")
  if (interaction_DC) {
    m_preds <- paste(m_preds, "+", paste(D, C, sep = ":", collapse = " + "))
  }

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

  m_forms <- lapply(M, function(x) paste(x, "~", m_preds))
  y_form <- paste(Y, "~", y_preds)

  m_models <- lapply(m_forms, function(x) lm(as.formula(x), data = df, weights = weights))
  names(m_models) <- M
  y_model <- lm(as.formula(y_form), data = df, weights = weights)

  if (interaction_DM) {
    NDE_part <- mapply(
      function(M_k, M_model_k) y_model$coef[[paste0(D, ":", M_k)]] * (M_model_k$coef[["(Intercept)"]] + M_model_k$coef[[D]] * dstar),
      M, m_models
    )
    NIE_part <- mapply(
      function(M_k, M_model_k) M_model_k$coef[[D]] * (y_model$coef[[M_k]] + y_model$coef[[paste0(D, ":", M_k)]] * d),
      M, m_models
    )
    CDE_part <- mapply(function(M_k, m_k) y_model$coef[[paste0(D, ":", M_k)]] * m_k, M, m)
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

# ============================================================
# Source estimator functions
# ============================================================
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

# ============================================================
# Single replication function
# ============================================================
run_one_rep <- function(rep_id, n = 2000, par = exp1_params) {
  # Set seed for reproducibility
  set.seed(rep_id)

  # Generate data
  df <- gen_exp1(n, par)

  # Initialize results
  results <- list(rep_id = rep_id)

  # ========== PATH-SPECIFIC EFFECTS ==========

  # 1. Linear estimator
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

  # 2. IPW estimator for PSE
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

  # 4. MedSim for PSE
  med_fit <- tryCatch({
    spec <- list(
      list(func = "glm", formula = M1 ~ D + C, args = list(family = binomial())),
      list(func = "glm", formula = M2 ~ D + M1 + C, args = list(family = binomial())),
      list(func = "lm", formula = Y ~ D + M1 + M2 + C)
    )
    medsim(data = df, num_sim = 5000, treatment = "D", intv_med = NULL,
           model_spec = spec, seed = rep_id, boot = FALSE)
  }, error = function(e) NULL)

  if (!is.null(med_fit)) {
    is_scalar <- vapply(med_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))
    scalars <- unlist(med_fit[is_scalar])
    grab <- function(pattern) {
      idx <- grep(pattern, names(scalars), fixed = FALSE)
      if (length(idx) == 0) return(NA_real_)
      scalars[idx[1]]
    }
    results$med_ATE <- grab("TE\\(")
    results$med_DY <- grab("PSE.*D\\s*->\\s*Y\\s*\\)")
    results$med_DM2Y <- grab("PSE.*M2\\s*->\\s*Y")
    results$med_DM1Y <- grab("PSE.*M1\\s*~>\\s*Y")
  } else {
    results$med_ATE <- results$med_DY <- results$med_DM2Y <- results$med_DM1Y <- NA
  }

  # ========== INTERVENTIONAL EFFECTS ==========

  # 1. IPW for Interventional Effects
  ipw_intv_fit <- tryCatch({
    ipwvent_inner(
      data = df, D = "D", M = "M2", Y = "Y", L = "M1",
      D_formula = as.formula("D ~ C"),
      L_formula = as.formula("M1 ~ D + C"),
      M_formula = as.formula("M2 ~ D + M1 + C"),
      stabilize = TRUE, censor = TRUE, minimal = TRUE
    )
  }, error = function(e) NULL)

  if (!is.null(ipw_intv_fit)) {
    results$ipw_intv_OE <- ipw_intv_fit$OE
    results$ipw_intv_IDE <- ipw_intv_fit$IDE
    results$ipw_intv_IIE <- ipw_intv_fit$IIE
  } else {
    results$ipw_intv_OE <- results$ipw_intv_IDE <- results$ipw_intv_IIE <- NA
  }

  # 2. RWR for Interventional Effects
  if (has_rwrmed) {
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

  # 3. MedSim for Interventional Effects
  medsim_intv_fit <- tryCatch({
    spec_intv <- list(
      list(func = "glm", formula = M1 ~ D + C, args = list(family = binomial())),
      list(func = "glm", formula = M2 ~ D + C, args = list(family = binomial())),
      list(func = "lm", formula = Y ~ D + M1 + M2 + C)
    )
    medsim(data = df, num_sim = 5000, treatment = "D", intv_med = "M2",
           model_spec = spec_intv, seed = rep_id, boot = FALSE)
  }, error = function(e) NULL)

  if (!is.null(medsim_intv_fit)) {
    is_scalar <- vapply(medsim_intv_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))
    scalars <- unlist(medsim_intv_fit[is_scalar])
    grab <- function(pattern) {
      idx <- grep(pattern, names(scalars), fixed = FALSE)
      if (length(idx) == 0) return(NA_real_)
      scalars[idx[1]]
    }
    results$med_intv_OE <- grab("OE\\(")
    results$med_intv_IDE <- grab("IDE\\(")
    results$med_intv_IIE <- grab("IIE\\(")
  } else {
    results$med_intv_OE <- results$med_intv_IDE <- results$med_intv_IIE <- NA
  }

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
