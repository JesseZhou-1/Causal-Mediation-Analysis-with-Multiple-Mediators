#!/usr/bin/env Rscript
# ============================================================
# Aggregate Monte Carlo Results for Experiment 4
# Nonlinear DGP: Count Treatment, Binary M1, Continuous M2, Ordinal Y
#
# Usage:
#   Rscript aggregate_results_exp4.R <output_dir> <n_reps> [--phase=1|2]
#
# --phase=1: Aggregate R estimator results only (run after Phase 1)
# --phase=2: Merge medflow results into existing R results (run after Phase 2)
# No flag:   Aggregate everything (backward compatible)
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript aggregate_results_exp4.R <output_dir> <n_reps> [--phase=1|2]")
}

output_dir <- args[1]
n_reps <- as.integer(args[2])

# Parse optional --phase flag
phase <- 0  # 0 = aggregate everything
if (length(args) >= 3) {
  phase_arg <- args[3]
  if (grepl("^--phase=", phase_arg)) {
    phase <- as.integer(sub("--phase=", "", phase_arg))
  }
}

cat("=" , rep("=", 60), "\n", sep = "")
if (phase == 1) {
  cat("Aggregating Experiment 4 — Phase 1 (R estimators only)\n")
} else if (phase == 2) {
  cat("Aggregating Experiment 4 — Phase 2 (adding medflow results)\n")
} else {
  cat("Aggregating Experiment 4 — All Results\n")
}
cat("=" , rep("=", 60), "\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

find_medflow_result_files <- function(output_dir) {
  candidate_dirs <- c(output_dir, file.path(output_dir, "rep_results"))
  candidate_dirs <- candidate_dirs[dir.exists(candidate_dirs)]

  if (length(candidate_dirs) == 0) {
    return(character(0))
  }

  files <- unlist(lapply(candidate_dirs, function(dir_path) {
    list.files(
      dir_path,
      pattern = "^medflow_rep_[0-9]+\\.csv$",
      full.names = TRUE
    )
  }), use.names = FALSE)

  files <- unique(files[file.exists(files)])
  if (length(files) == 0) {
    return(character(0))
  }

  file_tbl <- data.frame(
    path = files,
    rep_id = as.integer(sub("^medflow_rep_([0-9]+)\\.csv$", "\\1", basename(files))),
    mtime = file.info(files)$mtime,
    in_rep_results = basename(dirname(files)) == "rep_results",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(.data$rep_id, dplyr::desc(.data$mtime), dplyr::desc(.data$in_rep_results)) %>%
    dplyr::distinct(.data$rep_id, .keep_all = TRUE) %>%
    dplyr::arrange(.data$rep_id)

  file_tbl$path
}

# ============================================================
# True values (computed via Monte Carlo with large N)
# ============================================================
invlogit <- function(x) 1 / (1 + exp(-x))

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

mu_M2_fn <- function(D_val, M1_val, C_val, par = exp4_params) {
  par$M2_gamma_D * D_val + par$M2_gamma_M1 * M1_val + par$M2_gamma_C * C_val +
  par$M2_gamma_D2 * D_val^2 + par$M2_gamma_DM1 * D_val * M1_val
}

lp_Y_fn <- function(D_val, M1_val, M2_val, C_val, par = exp4_params) {
  par$Y_delta_D * D_val + par$Y_delta_M1 * M1_val + par$Y_delta_M2 * M2_val +
  par$Y_delta_C * C_val + par$Y_delta_M1M2 * M1_val * M2_val +
  par$Y_delta_sinM2 * sin(pi * M2_val / 2) + par$Y_delta_DC * D_val * C_val
}

compute_true_pse <- function(N_truth = 2000000L, par = exp4_params, d = 5, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)

  M1_d     <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * d +
               par$M1_beta_C * C + par$M1_beta_DC * d * C))
  M1_dstar <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
               par$M1_beta_C * C + par$M1_beta_DC * dstar * C))

  M2_d_M1d      <- rnorm(N_truth, mu_M2_fn(d, M1_d, C, par), par$M2_sigma)
  M2_dstar_M1ds <- rnorm(N_truth, mu_M2_fn(dstar, M1_dstar, C, par), par$M2_sigma)
  M2_d_M1ds     <- rnorm(N_truth, mu_M2_fn(d, M1_dstar, C, par), par$M2_sigma)

  Y_d_M1d_M2d    <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_d_M1d, C, par), par$Y_cuts)
  Y_dstar        <- rordlogit(N_truth, lp_Y_fn(dstar, M1_dstar, M2_dstar_M1ds, C, par), par$Y_cuts)
  Y_d_M1ds_M2ds  <- rordlogit(N_truth, lp_Y_fn(d, M1_dstar, M2_dstar_M1ds, C, par), par$Y_cuts)
  Y_d_M1ds_M2d   <- rordlogit(N_truth, lp_Y_fn(d, M1_dstar, M2_d_M1ds, C, par), par$Y_cuts)

  c(ATE  = mean(Y_d_M1d_M2d) - mean(Y_dstar),
    DY   = mean(Y_d_M1ds_M2ds) - mean(Y_dstar),
    DM2Y = mean(Y_d_M1ds_M2d) - mean(Y_d_M1ds_M2ds),
    DM1Y = mean(Y_d_M1d_M2d) - mean(Y_d_M1ds_M2d))
}

compute_true_intv <- function(N_truth = 2000000L, par = exp4_params, d = 5, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)

  M1_d     <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * d +
               par$M1_beta_C * C + par$M1_beta_DC * d * C))
  M1_dstar <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
               par$M1_beta_C * C + par$M1_beta_DC * dstar * C))

  M2_d_natural     <- rnorm(N_truth, mu_M2_fn(d, M1_d, C, par), par$M2_sigma)
  M2_dstar_natural <- rnorm(N_truth, mu_M2_fn(dstar, M1_dstar, C, par), par$M2_sigma)

  Y_d_natural     <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_d_natural, C, par), par$Y_cuts)
  Y_dstar_natural <- rordlogit(N_truth, lp_Y_fn(dstar, M1_dstar, M2_dstar_natural, C, par), par$Y_cuts)

  # M2*(dstar): marginalize M2 over M1|D=dstar (re-draw M1 for independence)
  M1_for_M2star <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
                    par$M1_beta_C * C + par$M1_beta_DC * dstar * C))
  M2_star_dstar <- rnorm(N_truth, mu_M2_fn(dstar, M1_for_M2star, C, par), par$M2_sigma)

  Y_d_M2star <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_star_dstar, C, par), par$Y_cuts)

  c(OE  = mean(Y_d_natural) - mean(Y_dstar_natural),
    IDE = mean(Y_d_M2star) - mean(Y_dstar_natural),
    IIE = mean(Y_d_natural) - mean(Y_d_M2star))
}

cat("Computing true values (2M samples)...\n")
true_pse <- compute_true_pse()
true_intv <- compute_true_intv()

cat("\nTrue Path-Specific Effects (d=5 vs dstar=0):\n")
print(round(true_pse, 4))

cat("\nTrue Interventional Effects:\n")
print(round(true_intv, 4))

# ============================================================
# Load R results
# ============================================================
recover_from_rep_checkpoints <- function(results_df, output_dir, n_reps) {
  chk_dir <- file.path(output_dir, "rep_results")
  if (!dir.exists(chk_dir)) {
    cat("No per-rep checkpoint directory found (rep_results/).\n")
    return(results_df)
  }

  chk_files <- list.files(chk_dir, pattern = "^rep_[0-9]+_result\\.rds$", full.names = TRUE)
  cat(sprintf("Found %d per-rep checkpoint files\n", length(chk_files)))
  if (length(chk_files) == 0) {
    return(results_df)
  }

  chk_list <- lapply(chk_files, function(f) {
    tryCatch(readRDS(f), error = function(e) NULL)
  })
  n_bad <- sum(vapply(chk_list, is.null, logical(1)))
  if (n_bad > 0) {
    cat(sprintf("Warning: %d checkpoint files could not be read\n", n_bad))
  }
  chk_list <- chk_list[!vapply(chk_list, is.null, logical(1))]
  if (length(chk_list) == 0) {
    return(results_df)
  }

  chk_df <- bind_rows(chk_list)
  if (!("rep_id" %in% names(chk_df))) {
    cat("Warning: rep_id not found in checkpoint rows; skipping recovery.\n")
    return(results_df)
  }
  chk_df <- chk_df %>%
    dplyr::arrange(.data$rep_id) %>%
    dplyr::distinct(.data$rep_id, .keep_all = TRUE)

  if (!("rep_id" %in% names(results_df))) {
    cat("Warning: rep_id not found in task-level results; using checkpoints only.\n")
    results_df <- chk_df
    cat(sprintf("Total R replications after checkpoint recovery: %d (expected: %d)\n",
                nrow(results_df), n_reps))
    return(results_df)
  }

  all_cols <- union(names(results_df), names(chk_df))
  miss_in_results <- setdiff(all_cols, names(results_df))
  miss_in_chk <- setdiff(all_cols, names(chk_df))
  if (length(miss_in_results) > 0) {
    for (col in miss_in_results) results_df[[col]] <- NA
  }
  if (length(miss_in_chk) > 0) {
    for (col in miss_in_chk) chk_df[[col]] <- NA
  }
  results_df <- results_df[, all_cols]
  chk_df <- chk_df[, all_cols]

  missing_rows <- chk_df %>% dplyr::filter(!(.data$rep_id %in% results_df$rep_id))
  if (nrow(missing_rows) > 0) {
    cat(sprintf("Recovered %d replications from per-rep checkpoints\n", nrow(missing_rows)))
    results_df <- bind_rows(results_df, missing_rows)
  }

  results_df <- results_df %>%
    dplyr::arrange(.data$rep_id) %>%
    dplyr::distinct(.data$rep_id, .keep_all = TRUE)

  cat(sprintf("Total R replications after checkpoint recovery: %d (expected: %d)\n",
              nrow(results_df), n_reps))
  results_df
}

if (phase == 2) {
  # Phase 2: load previously saved R results from Phase 1
  r_results_file <- file.path(output_dir, "r_results_exp4.rds")
  if (!file.exists(r_results_file)) {
    stop("Phase 1 results not found at: ", r_results_file,
         "\nRun with --phase=1 first.")
  }
  results_df <- readRDS(r_results_file)
  cat(sprintf("\nLoaded Phase 1 R results: %d replications\n", nrow(results_df)))
} else {
  # Phase 1 or full: load from raw .rds files
  cat("\nLoading R results from:", output_dir, "\n")

  result_files <- list.files(output_dir, pattern = "^results_task_.*\\.rds$", full.names = TRUE)
  cat(sprintf("Found %d R result files\n", length(result_files)))

  if (length(result_files) == 0) {
    stop("No R result files found!")
  }

  results_list <- lapply(result_files, readRDS)
  results_df <- bind_rows(results_list)

  cat(sprintf("Total R replications: %d (expected: %d)\n", nrow(results_df), n_reps))
}

results_df <- recover_from_rep_checkpoints(results_df, output_dir, n_reps)

# Save R results so Phase 2 can load them directly
r_results_file <- file.path(output_dir, "r_results_exp4.rds")
saveRDS(results_df, r_results_file)
cat(sprintf("R results saved to: %s\n", r_results_file))

# Report captured replication-level failures (if worker includes status columns)
if ("status" %in% names(results_df)) {
  n_err <- sum(results_df$status == "error", na.rm = TRUE)
  cat(sprintf("Captured failed replications in R worker output: %d\n", n_err))
  if (n_err > 0) {
    err_cols <- c("rep_id", "status", "error_message", "error_call", "error_log")
    err_cols <- err_cols[err_cols %in% names(results_df)]
    err_df <- results_df %>%
      dplyr::filter(.data$status == "error") %>%
      dplyr::select(dplyr::all_of(err_cols))
    err_file <- file.path(output_dir, "rep_errors_exp4.csv")
    write.csv(err_df, err_file, row.names = FALSE)
    cat(sprintf("Saved replication error summary to: %s\n", err_file))
    cat(sprintf("First failed rep_ids: %s\n", paste(head(err_df$rep_id, 20), collapse = ", ")))
  }
}

# ============================================================
# Load and merge medflow results (Phase 2 or full)
# ============================================================
if (phase != 1) {
  cat("\nLoading medflow results...\n")

  medflow_files <- find_medflow_result_files(output_dir)
  cat(sprintf("Found %d medflow result files\n", length(medflow_files)))

  if (length(medflow_files) > 0) {
    medflow_list <- lapply(medflow_files, function(path) {
      tryCatch(read.csv(path), error = function(e) NULL)
    })
    medflow_list <- medflow_list[!vapply(medflow_list, is.null, logical(1))]
    medflow_df <- bind_rows(medflow_list)
    if ("rep_id" %in% names(medflow_df)) {
      medflow_df <- medflow_df %>%
        dplyr::arrange(.data$rep_id) %>%
        dplyr::distinct(.data$rep_id, .keep_all = TRUE)
      cat(sprintf("Total medflow replications: %d\n", nrow(medflow_df)))

      # Merge medflow results into main results by rep_id
      results_df <- left_join(results_df, medflow_df, by = "rep_id")
    } else {
      cat("Warning: medflow files were found but none could be read into a table with rep_id.\n")
    }
  } else {
    cat("Warning: No medflow results found. Proceeding with R estimators only.\n")
  }
}

# ============================================================
# Compute Bias and RMSE for Path-Specific Effects
# ============================================================
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("PATH-SPECIFIC EFFECTS RESULTS (d=5 vs dstar=0)\n")
cat("=" , rep("=", 60), "\n\n")

if (phase == 1) {
  pse_estimators <- c("lin", "ipw", "reg", "med", "dml", "ipw_sl", "reg_sl")
  pse_estimator_labels <- c("LIN", "IPW", "REG", "MED", "DML", "IPW-SL", "REG-SL")
} else {
  pse_estimators <- c("lin", "ipw", "reg", "med", "dml", "ipw_sl", "reg_sl", "mf")
  pse_estimator_labels <- c("LIN", "IPW", "REG", "MED", "DML", "IPW-SL", "REG-SL", "MF")
}
pse_estimands <- c("ATE", "DY", "DM2Y", "DM1Y")
pse_estimand_labels <- c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y")

pse_stats <- data.frame()

for (j in seq_along(pse_estimators)) {
  est <- pse_estimators[j]
  est_label <- pse_estimator_labels[j]

  for (i in seq_along(pse_estimands)) {
    estimand <- pse_estimands[i]
    col_name <- paste0(est, "_", estimand)
    true_val <- true_pse[estimand]

    if (col_name %in% names(results_df)) {
      estimates <- results_df[[col_name]]
      n_valid <- sum(!is.na(estimates))
      if (n_valid > 0) {
        mean_est <- mean(estimates, na.rm = TRUE)
        bias <- mean_est - true_val
        variance <- var(estimates, na.rm = TRUE)
        rmse <- sqrt(mean((estimates - true_val)^2, na.rm = TRUE))

        pse_stats <- rbind(pse_stats, data.frame(
          Estimator = est_label,
          Estimand = pse_estimand_labels[i],
          Truth = round(true_val, 4),
          Mean = round(mean_est, 4),
          Bias = round(bias, 4),
          Variance = round(variance, 4),
          RMSE = round(rmse, 4),
          N_valid = n_valid
        ))
      }
    }
  }
}

print(pse_stats, row.names = FALSE)

cat("\n--- Bias Summary (PSE) ---\n")
bias_wide_pse <- pse_stats %>%
  select(Estimator, Estimand, Bias) %>%
  pivot_wider(names_from = Estimator, values_from = Bias)
print(as.data.frame(bias_wide_pse), row.names = FALSE)

cat("\n--- RMSE Summary (PSE) ---\n")
rmse_wide_pse <- pse_stats %>%
  select(Estimator, Estimand, RMSE) %>%
  pivot_wider(names_from = Estimator, values_from = RMSE)
print(as.data.frame(rmse_wide_pse), row.names = FALSE)

# ============================================================
# Compute Bias and RMSE for Interventional Effects
# ============================================================
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("INTERVENTIONAL EFFECTS RESULTS\n")
cat("=" , rep("=", 60), "\n\n")

if (phase == 1) {
  intv_estimators <- c("rwr_intv", "med_intv")
  intv_estimator_labels <- c("RWR", "MedSim")
} else {
  intv_estimators <- c("rwr_intv", "med_intv", "mf_intv")
  intv_estimator_labels <- c("RWR", "MedSim", "MedFlow")
}
intv_estimands <- c("OE", "IDE", "IIE")

intv_stats <- data.frame()

for (j in seq_along(intv_estimators)) {
  est <- intv_estimators[j]
  est_label <- intv_estimator_labels[j]

  for (estimand in intv_estimands) {
    col_name <- paste0(est, "_", estimand)
    true_val <- true_intv[estimand]

    if (col_name %in% names(results_df)) {
      estimates <- results_df[[col_name]]
      n_valid <- sum(!is.na(estimates))
      if (n_valid > 0) {
        mean_est <- mean(estimates, na.rm = TRUE)
        bias <- mean_est - true_val
        variance <- var(estimates, na.rm = TRUE)
        rmse <- sqrt(mean((estimates - true_val)^2, na.rm = TRUE))

        intv_stats <- rbind(intv_stats, data.frame(
          Estimator = est_label,
          Estimand = estimand,
          Truth = round(true_val, 4),
          Mean = round(mean_est, 4),
          Bias = round(bias, 4),
          Variance = round(variance, 4),
          RMSE = round(rmse, 4),
          N_valid = n_valid
        ))
      }
    }
  }
}

print(intv_stats, row.names = FALSE)

cat("\n--- Bias Summary (Interventional) ---\n")
bias_wide_intv <- intv_stats %>%
  select(Estimator, Estimand, Bias) %>%
  pivot_wider(names_from = Estimator, values_from = Bias)
print(as.data.frame(bias_wide_intv), row.names = FALSE)

cat("\n--- RMSE Summary (Interventional) ---\n")
rmse_wide_intv <- intv_stats %>%
  select(Estimator, Estimand, RMSE) %>%
  pivot_wider(names_from = Estimator, values_from = RMSE)
print(as.data.frame(rmse_wide_intv), row.names = FALSE)

# ============================================================
# Save summary results
# ============================================================
if (phase == 1) {
  suffix <- "_phase1"
} else if (phase == 2) {
  suffix <- ""  # Phase 2 is the final result
} else {
  suffix <- ""
}

summary_file <- file.path(output_dir, paste0("mc_summary_exp4", suffix, ".rds"))
summary_list <- list(
  true_pse = true_pse,
  true_intv = true_intv,
  pse_stats = pse_stats,
  intv_stats = intv_stats,
  n_reps = nrow(results_df),
  results_df = results_df
)
saveRDS(summary_list, summary_file)
cat(sprintf("\nSummary saved to: %s\n", summary_file))

write.csv(pse_stats, file.path(output_dir, paste0("pse_stats_exp4", suffix, ".csv")), row.names = FALSE)
write.csv(intv_stats, file.path(output_dir, paste0("intv_stats_exp4", suffix, ".csv")), row.names = FALSE)
cat(sprintf("CSV files saved to: %s\n", output_dir))

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
if (phase == 1) {
  cat("PHASE 1 AGGREGATION COMPLETE\n")
} else if (phase == 2) {
  cat("FINAL AGGREGATION COMPLETE (R + medflow)\n")
} else {
  cat("AGGREGATION COMPLETE\n")
}
cat("=" , rep("=", 60), "\n")
