#!/usr/bin/env Rscript
# ============================================================
# Aggregate Monte Carlo Results for Experiment 2
# Ordinal Treatment, Mediators, Outcome
# Computes Bias and RMSE for each estimator and estimand
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript aggregate_results_exp2.R <output_dir> <n_reps>")
}

output_dir <- args[1]
n_reps <- as.integer(args[2])

cat("=" , rep("=", 60), "\n", sep = "")
cat("Aggregating Experiment 2 Monte Carlo Results\n")
cat("=" , rep("=", 60), "\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

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

exp2_params <- list(
  D_cuts = c(-0.5, 0.5), D_beta = 0.7,
  M1_cuts = c(-0.8, 0.8), M1_beta_D = 0.9, M1_beta_C = 0.6,
  M2_cuts = c(-0.5, 1.0), M2_beta_D = 0.7, M2_beta_M1 = 1.0, M2_beta_C = 0.5,
  Y_cuts = c(-0.3, 0.9), Y_beta_D = 0.5, Y_beta_M1 = 0.6, Y_beta_M2 = 0.7, Y_beta_C = 0.4
)

compute_true_pse <- function(N_truth = 1000000L, par = exp2_params, d = 2, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)
  M1_d <- rordlogit(N_truth, par$M1_beta_D * d + par$M1_beta_C * C, par$M1_cuts)
  M1_dstar <- rordlogit(N_truth, par$M1_beta_D * dstar + par$M1_beta_C * C, par$M1_cuts)
  M2_d_M1d <- rordlogit(N_truth, par$M2_beta_D * d + par$M2_beta_M1 * M1_d + par$M2_beta_C * C, par$M2_cuts)
  M2_dstar_M1dstar <- rordlogit(N_truth, par$M2_beta_D * dstar + par$M2_beta_M1 * M1_dstar + par$M2_beta_C * C, par$M2_cuts)
  M2_d_M1dstar <- rordlogit(N_truth, par$M2_beta_D * d + par$M2_beta_M1 * M1_dstar + par$M2_beta_C * C, par$M2_cuts)
  Y_d_M1d_M2d <- rordlogit(N_truth, par$Y_beta_D * d + par$Y_beta_M1 * M1_d + par$Y_beta_M2 * M2_d_M1d + par$Y_beta_C * C, par$Y_cuts)
  Y_dstar <- rordlogit(N_truth, par$Y_beta_D * dstar + par$Y_beta_M1 * M1_dstar + par$Y_beta_M2 * M2_dstar_M1dstar + par$Y_beta_C * C, par$Y_cuts)
  Y_d_M1dstar_M2dstar <- rordlogit(N_truth, par$Y_beta_D * d + par$Y_beta_M1 * M1_dstar + par$Y_beta_M2 * M2_dstar_M1dstar + par$Y_beta_C * C, par$Y_cuts)
  Y_d_M1dstar_M2d <- rordlogit(N_truth, par$Y_beta_D * d + par$Y_beta_M1 * M1_dstar + par$Y_beta_M2 * M2_d_M1dstar + par$Y_beta_C * C, par$Y_cuts)

  c(ATE = mean(Y_d_M1d_M2d) - mean(Y_dstar),
    DY = mean(Y_d_M1dstar_M2dstar) - mean(Y_dstar),
    DM2Y = mean(Y_d_M1dstar_M2d) - mean(Y_d_M1dstar_M2dstar),
    DM1Y = mean(Y_d_M1d_M2d) - mean(Y_d_M1dstar_M2d))
}

compute_true_intv <- function(N_truth = 1000000L, par = exp2_params, d = 2, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)
  M1_d <- rordlogit(N_truth, par$M1_beta_D * d + par$M1_beta_C * C, par$M1_cuts)
  M1_dstar <- rordlogit(N_truth, par$M1_beta_D * dstar + par$M1_beta_C * C, par$M1_cuts)
  M2_d_natural <- rordlogit(N_truth, par$M2_beta_D * d + par$M2_beta_M1 * M1_d + par$M2_beta_C * C, par$M2_cuts)
  M2_dstar_natural <- rordlogit(N_truth, par$M2_beta_D * dstar + par$M2_beta_M1 * M1_dstar + par$M2_beta_C * C, par$M2_cuts)
  Y_d_natural <- rordlogit(N_truth, par$Y_beta_D * d + par$Y_beta_M1 * M1_d + par$Y_beta_M2 * M2_d_natural + par$Y_beta_C * C, par$Y_cuts)
  Y_dstar_natural <- rordlogit(N_truth, par$Y_beta_D * dstar + par$Y_beta_M1 * M1_dstar + par$Y_beta_M2 * M2_dstar_natural + par$Y_beta_C * C, par$Y_cuts)
  M1_for_M2star <- rordlogit(N_truth, par$M1_beta_D * dstar + par$M1_beta_C * C, par$M1_cuts)
  M2_star_dstar <- rordlogit(N_truth, par$M2_beta_D * dstar + par$M2_beta_M1 * M1_for_M2star + par$M2_beta_C * C, par$M2_cuts)
  Y_d_M2star <- rordlogit(N_truth, par$Y_beta_D * d + par$Y_beta_M1 * M1_d + par$Y_beta_M2 * M2_star_dstar + par$Y_beta_C * C, par$Y_cuts)

  c(OE = mean(Y_d_natural) - mean(Y_dstar_natural),
    IDE = mean(Y_d_M2star) - mean(Y_dstar_natural),
    IIE = mean(Y_d_natural) - mean(Y_d_M2star))
}

cat("Computing true values (1M samples)...\n")
true_pse <- compute_true_pse()
true_intv <- compute_true_intv()

cat("\nTrue Path-Specific Effects (d=2 vs dstar=0):\n")
print(round(true_pse, 4))

cat("\nTrue Interventional Effects:\n")
print(round(true_intv, 4))

# ============================================================
# Load and combine results
# ============================================================
cat("\nLoading results from:", output_dir, "\n")

result_files <- list.files(output_dir, pattern = "^results_task_.*\\.rds$", full.names = TRUE)
cat(sprintf("Found %d result files\n", length(result_files)))

if (length(result_files) == 0) {
  stop("No result files found!")
}

results_list <- lapply(result_files, readRDS)
results_df <- bind_rows(results_list)

cat(sprintf("Total replications: %d (expected: %d)\n", nrow(results_df), n_reps))

# ============================================================
# Compute Bias and RMSE for Path-Specific Effects
# ============================================================
cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("PATH-SPECIFIC EFFECTS RESULTS (d=2 vs dstar=0)\n")
cat("=" , rep("=", 60), "\n\n")

pse_estimators <- c("lin", "ipw", "reg", "med")
pse_estimands <- c("ATE", "DY", "DM2Y", "DM1Y")
pse_estimand_labels <- c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y")

pse_stats <- data.frame()

for (est in pse_estimators) {
  for (i in seq_along(pse_estimands)) {
    estimand <- pse_estimands[i]
    col_name <- paste0(est, "_", estimand)
    true_val <- true_pse[estimand]

    if (col_name %in% names(results_df)) {
      estimates <- results_df[[col_name]]
      n_valid <- sum(!is.na(estimates))
      mean_est <- mean(estimates, na.rm = TRUE)
      bias <- mean_est - true_val
      variance <- var(estimates, na.rm = TRUE)
      rmse <- sqrt(mean((estimates - true_val)^2, na.rm = TRUE))

      pse_stats <- rbind(pse_stats, data.frame(
        Estimator = toupper(est),
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

intv_estimators <- c("ipw_intv", "rwr_intv", "med_intv")
intv_estimator_labels <- c("IPW", "RWR", "MedSim")
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
summary_file <- file.path(output_dir, "mc_summary_exp2.rds")
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

write.csv(pse_stats, file.path(output_dir, "pse_stats_exp2.csv"), row.names = FALSE)
write.csv(intv_stats, file.path(output_dir, "intv_stats_exp2.csv"), row.names = FALSE)
cat(sprintf("CSV files saved to: %s\n", output_dir))

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("AGGREGATION COMPLETE\n")
cat("=" , rep("=", 60), "\n")
