#!/usr/bin/env Rscript
# ============================================================
# Monte Carlo R Worker Script for Experiment 4
# Count Treatment (Poisson), Binary M1, Continuous M2, Ordinal Y
# Nonlinear DGP with interactions
# Runs R-based estimators and saves data CSV for medflow
# ============================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript mc_exp4_r_worker.R <task_id> <n_reps> <n_nodes> <n_cores> <output_dir>")
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
  library(SuperLearner)
  library(ranger)
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
# DGP: Nonlinear, Count Treatment, Mixed Variable Types
# ============================================================
rordlogit <- function(n, linear_pred, cuts) {
  u <- runif(n)
  p0 <- invlogit(cuts[1] - linear_pred)
  p1 <- invlogit(cuts[2] - linear_pred)
  ifelse(u < p0, 0L, ifelse(u < p1, 1L, 2L))
}

exp4_params <- list(
  D_alpha0 = 0.7, D_alpha1 = 0.3, D_alpha2 = 0.15,
  M1_beta0 = -1.5, M1_beta_D = 0.5, M1_beta_C = 0.4, M1_beta_DC = 0.25,
  M2_gamma_D = 0.4, M2_gamma_M1 = 0.6, M2_gamma_C = 0.3,
  M2_gamma_D2 = 0.15, M2_gamma_DM1 = -0.1, M2_sigma = 1.0,
  Y_delta_D = 0.3, Y_delta_M1 = 0.4, Y_delta_M2 = 0.5, Y_delta_C = 0.3,
  Y_delta_M1M2 = 0.2, Y_delta_sinM2 = 0.15, Y_delta_DC = 0.1,
  Y_cuts = c(-0.5, 0.8)
)

gen_exp4 <- function(n, par = exp4_params) {
  C <- rnorm(n)
  D <- rpois(n, exp(par$D_alpha0 + par$D_alpha1 * C + par$D_alpha2 * C^2))
  M1 <- rbinom(n, 1, plogis(par$M1_beta0 + par$M1_beta_D * D +
         par$M1_beta_C * C + par$M1_beta_DC * D * C))
  mu_M2 <- par$M2_gamma_D * D + par$M2_gamma_M1 * M1 + par$M2_gamma_C * C +
            par$M2_gamma_D2 * D^2 + par$M2_gamma_DM1 * D * M1
  M2 <- rnorm(n, mu_M2, par$M2_sigma)
  lp_Y <- par$Y_delta_D * D + par$Y_delta_M1 * M1 + par$Y_delta_M2 * M2 +
           par$Y_delta_C * C + par$Y_delta_M1M2 * M1 * M2 +
           par$Y_delta_sinM2 * sin(pi * M2 / 2) + par$Y_delta_DC * D * C
  Y <- rordlogit(n, lp_Y, par$Y_cuts)
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
source("ipwmed.R")
source("ipwpath.R")
source("linpath.R")
source("medsim.R")
source("pathimp.R")

# Source rwrlite if rwrmed is available
if (has_rwrmed) {
  source("rwrlite.R")
}

source("dmlpath_custom.R")

cat("All functions loaded.\n")

# ============================================================
# Single replication function
# ============================================================
d_val <- 5
dstar_val <- 0

# Full result schema so every replication returns one row, even on failure.
result_cols <- c(
  "rep_id",
  "lin_ATE", "lin_DY", "lin_DM2Y", "lin_DM1Y",
  "ipw_ATE", "ipw_DY", "ipw_DM2Y", "ipw_DM1Y",
  "reg_ATE", "reg_DY", "reg_DM2Y", "reg_DM1Y",
  "med_ATE", "med_DY", "med_DM2Y", "med_DM1Y",
  "dml_ATE", "dml_DY", "dml_DM2Y", "dml_DM1Y",
  "ipw_sl_ATE", "ipw_sl_DY", "ipw_sl_DM2Y", "ipw_sl_DM1Y",
  "reg_sl_ATE", "reg_sl_DY", "reg_sl_DM2Y", "reg_sl_DM1Y",
  "rwr_intv_OE", "rwr_intv_IDE", "rwr_intv_IIE",
  "med_intv_OE", "med_intv_IDE", "med_intv_IIE",
  "status", "error_message", "error_call", "error_log"
)

make_result_row <- function(
    rep_id,
    status = "ok",
    error_message = NA_character_,
    error_call = NA_character_,
    error_log = NA_character_
) {
  out <- as.list(setNames(rep(NA_real_, length(result_cols)), result_cols))
  out$rep_id <- as.integer(rep_id)
  out$status <- status
  out$error_message <- error_message
  out$error_call <- error_call
  out$error_log <- error_log
  as.data.frame(out, stringsAsFactors = FALSE, check.names = FALSE)
}

# Defensive cleanup for scalar outputs that may carry invalid names attributes
# (e.g., names = NA), which can break as.data.frame(list(...)).
strip_scalar_names <- function(x) {
  if (!is.list(x) && length(x) == 1L) {
    return(unname(x))
  }
  x
}

# Build deterministic starting values for MASS::polr from outcome prevalence.
build_polr_start <- function(formula, data, method = "logistic") {
  y_name <- all.vars(formula)[1]
  y <- data[[y_name]]
  if (!is.factor(y)) {
    y <- factor(y)
  }
  k <- nlevels(y)
  if (k < 3) {
    return(NULL)
  }

  mm <- model.matrix(terms(formula, data = data), data = data)
  if (ncol(mm) > 0 && "(Intercept)" %in% colnames(mm)) {
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  }
  n_beta <- ncol(mm)

  probs <- as.numeric(prop.table(table(y)))
  cum_probs <- cumsum(probs)[seq_len(k - 1)]
  eps <- 1e-4
  cum_probs <- pmin(pmax(cum_probs, eps), 1 - eps)

  zeta_start <- if (identical(method, "probit")) qnorm(cum_probs) else qlogis(cum_probs)
  c(rep(0, n_beta), zeta_start)
}

# Retry medsim once with explicit polr starts when default polr init fails.
run_medsim_with_polr_retry <- function(
    df_med,
    spec,
    treatment,
    intv_med,
    seed,
    rep_id,
    label,
    num_sim = 5000,
    cat_list = c("0", "5")
) {
  fit0 <- tryCatch(
    medsim(
      data = df_med,
      num_sim = num_sim,
      cat_list = cat_list,
      treatment = treatment,
      intv_med = intv_med,
      model_spec = spec,
      seed = seed,
      boot = FALSE
    ),
    error = function(e) e
  )
  if (!inherits(fit0, "error")) {
    return(fit0)
  }

  err_msg <- conditionMessage(fit0)
  if (!grepl("attempt to find suitable starting values failed", err_msg, fixed = TRUE)) {
    return(NULL)
  }

  y_idx <- length(spec)
  if (!identical(spec[[y_idx]]$func, "polr")) {
    return(NULL)
  }

  polr_method <- spec[[y_idx]]$args$method
  if (is.null(polr_method)) {
    polr_method <- "logistic"
  }

  start_vec <- tryCatch(
    build_polr_start(spec[[y_idx]]$formula, df_med, method = polr_method),
    error = function(e) NULL
  )
  if (is.null(start_vec)) {
    return(NULL)
  }

  retry_spec <- spec
  retry_args <- retry_spec[[y_idx]]$args
  retry_args$start <- start_vec
  retry_args$Hess <- TRUE
  retry_spec[[y_idx]]$args <- retry_args

  cat(sprintf("[Rep %d] %s retrying medsim polr with explicit start values.\n", rep_id, label))

  tryCatch(
    medsim(
      data = df_med,
      num_sim = num_sim,
      cat_list = cat_list,
      treatment = treatment,
      intv_med = intv_med,
      model_spec = retry_spec,
      seed = seed,
      boot = FALSE
    ),
    error = function(e) {
      cat(sprintf("[Rep %d] %s retry failed: %s\n", rep_id, label, conditionMessage(e)))
      NULL
    }
  )
}

# Save one-row replication checkpoint so completed work is recoverable even if
# task-level aggregation fails later.
rep_result_dir <- file.path(output_dir, "rep_results")
dir.create(rep_result_dir, showWarnings = FALSE, recursive = TRUE)

save_rep_checkpoint <- function(rep_id, rep_row) {
  chk_file <- file.path(rep_result_dir, sprintf("rep_%d_result.rds", rep_id))
  tryCatch(
    {
      saveRDS(rep_row, chk_file)
    },
    error = function(e) {
      cat(sprintf("[Rep %d] WARNING: failed to write checkpoint: %s\n",
                  rep_id, conditionMessage(e)))
    }
  )
}

run_one_rep <- function(rep_id, n = 20000, par = exp4_params) {
  set.seed(rep_id)

  df <- gen_exp4(n, par)
  results <- list(rep_id = rep_id)

  # ---- Save data CSV for medflow (Python) ----
  rep_dir <- file.path(output_dir, sprintf("rep_%d", rep_id))
  dir.create(rep_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(df, file.path(rep_dir, "data.csv"), row.names = FALSE)

  # ---- PSE Estimators ----

  # 1. Linear
  lin_fit <- tryCatch({
    linpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
            d = d_val, dstar = dstar_val, boot = FALSE)
  }, error = function(e) NULL)
  if (!is.null(lin_fit)) {
    results$lin_ATE <- as.numeric(lin_fit$ATE); results$lin_DY <- as.numeric(lin_fit$PSE[["D->Y"]])
    results$lin_DM2Y <- as.numeric(lin_fit$PSE[["D->M2->Y"]]); results$lin_DM1Y <- as.numeric(lin_fit$PSE[["D->M1~>Y"]])
  } else { results$lin_ATE <- results$lin_DY <- results$lin_DM2Y <- results$lin_DM1Y <- NA }

  # 2. IPW (Poisson auto-detected)
  ipw_fit <- tryCatch({
    ipwpath(data = df, D = "D", M = list("M1", "M2"), Y = "Y", C = "C",
            d = d_val, dstar = dstar_val, boot = FALSE)
  }, error = function(e) NULL)
  if (!is.null(ipw_fit)) {
    results$ipw_ATE <- as.numeric(ipw_fit$ATE); results$ipw_DY <- as.numeric(ipw_fit$PSE[["D->Y"]])
    results$ipw_DM2Y <- as.numeric(ipw_fit$PSE[["D->M2->Y"]]); results$ipw_DM1Y <- as.numeric(ipw_fit$PSE[["D->M1~>Y"]])
  } else { results$ipw_ATE <- results$ipw_DY <- results$ipw_DM2Y <- results$ipw_DM1Y <- NA }

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
  } else { results$reg_ATE <- results$reg_DY <- results$reg_DM2Y <- results$reg_DM1Y <- NA }

  # 4. MedSim (PSE, simple additive parametric specs)
  med_fit <- tryCatch({
    df_med <- df
    df_med$M1 <- factor(df_med$M1, ordered = FALSE)
    df_med$Y  <- factor(df_med$Y, ordered = FALSE)
    spec <- list(
      list(func = "glm",  formula = M1 ~ D + C, args = list(family = binomial())),
      list(func = "lm",   formula = M2 ~ D + M1 + C),
      list(func = "polr", formula = Y ~ D + M1 + M2 + C, args = list(method = "logistic"))
    )
    run_medsim_with_polr_retry(
      df_med = df_med,
      spec = spec,
      treatment = "D",
      intv_med = NULL,
      seed = rep_id,
      rep_id = rep_id,
      label = "MedSim-PSE"
    )
  }, error = function(e) NULL)
  if (!is.null(med_fit)) {
    scalars <- unlist(med_fit[vapply(med_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))])
    grab <- function(pattern) { idx <- grep(pattern, names(scalars)); if (length(idx) == 0) NA_real_ else scalars[idx[1]] }
    results$med_ATE <- grab("TE\\("); results$med_DY <- grab("PSE.*D\\s*->\\s*Y\\s*\\)")
    results$med_DM2Y <- grab("PSE.*M2\\s*->\\s*Y"); results$med_DM1Y <- grab("PSE.*M1\\s*~>\\s*Y")
  } else { results$med_ATE <- results$med_DY <- results$med_DM2Y <- results$med_DM1Y <- NA }

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
    results$dml_ATE <- results$dml_DY <- NA
    results$dml_DM2Y <- results$dml_DM1Y <- NA
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
    results$ipw_sl_ATE <- results$ipw_sl_DY <- NA
    results$ipw_sl_DM2Y <- results$ipw_sl_DM1Y <- NA
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
    results$reg_sl_ATE <- results$reg_sl_DY <- NA
    results$reg_sl_DM2Y <- results$reg_sl_DM1Y <- NA
  }

  # ---- Interventional Effects (no IPW interventional) ----

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

  # MedSim Interventional
  medsim_intv_fit <- tryCatch({
    df_med <- df
    df_med$M1 <- factor(df_med$M1, ordered = FALSE)
    df_med$Y  <- factor(df_med$Y, ordered = FALSE)
    spec_intv <- list(
      list(func = "glm",  formula = M1 ~ D + C, args = list(family = binomial())),
      list(func = "lm",   formula = M2 ~ D + C),
      list(func = "polr", formula = Y ~ D + M1 + M2 + C, args = list(method = "logistic"))
    )
    run_medsim_with_polr_retry(
      df_med = df_med,
      spec = spec_intv,
      treatment = "D",
      intv_med = "M2",
      seed = rep_id + 10000,
      rep_id = rep_id,
      label = "MedSim-INTV"
    )
  }, error = function(e) NULL)
  if (!is.null(medsim_intv_fit)) {
    scalars <- unlist(medsim_intv_fit[vapply(medsim_intv_fit, function(x) is.numeric(x) && length(x) == 1, logical(1))])
    grab <- function(pattern) { idx <- grep(pattern, names(scalars)); if (length(idx) == 0) NA_real_ else scalars[idx[1]] }
    results$med_intv_OE <- grab("OE\\("); results$med_intv_IDE <- grab("IDE\\("); results$med_intv_IIE <- grab("IIE\\(")
  } else { results$med_intv_OE <- results$med_intv_IDE <- results$med_intv_IIE <- NA }

  results$status <- "ok"
  results$error_message <- NA_character_
  results$error_call <- NA_character_
  results$error_log <- NA_character_
  results <- lapply(results, strip_scalar_names)

  results_df <- as.data.frame(results, stringsAsFactors = FALSE, check.names = FALSE)
  missing_cols <- setdiff(result_cols, names(results_df))
  if (length(missing_cols) > 0) {
    for (col in missing_cols) {
      results_df[[col]] <- NA
    }
  }
  results_df <- results_df[, result_cols]

  return(results_df)
}

run_one_rep_safe <- function(rep_id, n = 20000, par = exp4_params) {
  tryCatch(
    {
      rep_row <- run_one_rep(rep_id = rep_id, n = n, par = par)
      save_rep_checkpoint(rep_id = rep_id, rep_row = rep_row)
      rep_row
    },
    error = function(e) {
      err_dir <- file.path(output_dir, "error_logs")
      dir.create(err_dir, showWarnings = FALSE, recursive = TRUE)

      err_call <- if (is.null(conditionCall(e))) {
        NA_character_
      } else {
        paste(deparse(conditionCall(e), width.cutoff = 500), collapse = " ")
      }
      err_file <- file.path(err_dir, sprintf("rep_%d_error.log", rep_id))
      err_trace <- paste(
        vapply(
          sys.calls(),
          function(cl) paste(deparse(cl, width.cutoff = 500), collapse = " "),
          character(1)
        ),
        collapse = "\n"
      )
      writeLines(
        c(
          sprintf("timestamp: %s", Sys.time()),
          sprintf("rep_id: %d", rep_id),
          sprintf("message: %s", conditionMessage(e)),
          sprintf("call: %s", err_call),
          "",
          "sys.calls():",
          err_trace
        ),
        con = err_file
      )

      cat(sprintf("[Rep %d] FATAL ERROR captured: %s\n", rep_id, conditionMessage(e)))
      cat(sprintf("[Rep %d] Error log: %s\n", rep_id, err_file))

      err_row <- make_result_row(
        rep_id = rep_id,
        status = "error",
        error_message = conditionMessage(e),
        error_call = err_call,
        error_log = err_file
      )
      save_rep_checkpoint(rep_id = rep_id, rep_row = err_row)
      err_row
    }
  )
}

# ============================================================
# Run replications in parallel
# ============================================================
cat(sprintf("Starting %d replications using %d cores...\n", length(my_reps), n_cores))

results_list <- mclapply(
  my_reps,
  run_one_rep_safe,
  mc.cores = n_cores,
  # More robust: don't batch many reps on one worker process.
  mc.preschedule = FALSE
)

# Guard against any unexpected non-data.frame return.
bad_idx <- which(!vapply(results_list, is.data.frame, logical(1)))
if (length(bad_idx) > 0) {
  for (idx in bad_idx) {
    rep_id <- my_reps[idx]
    bad_class <- paste(class(results_list[[idx]]), collapse = ", ")
    bad_msg <- sprintf("Unexpected return type from run_one_rep_safe: %s", bad_class)
    cat(sprintf("[Rep %d] %s\n", rep_id, bad_msg))
    results_list[[idx]] <- make_result_row(
      rep_id = rep_id,
      status = "error",
      error_message = bad_msg
    )
  }
}

# Combine results
results_df <- bind_rows(results_list)

# Save results
output_file <- file.path(output_dir, sprintf("results_task_%d.rds", task_id))
saveRDS(results_df, output_file)

cat(sprintf("Task %d completed at %s\n", task_id, Sys.time()))
cat(sprintf("Results saved to: %s\n", output_file))
cat(sprintf("Completed %d replications\n", nrow(results_df)))
cat(sprintf("Replications with captured errors: %d\n", sum(results_df$status == "error", na.rm = TRUE)))
