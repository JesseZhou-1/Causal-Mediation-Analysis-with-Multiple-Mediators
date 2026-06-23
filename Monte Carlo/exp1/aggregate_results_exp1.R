#!/usr/bin/env Rscript
# ============================================================
# Aggregate Monte Carlo Results
# Computes Bias and RMSE for each estimator and estimand
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript aggregate_results_exp1.R <output_dir> <n_reps>")
}

output_dir <- args[1]
n_reps <- as.integer(args[2])

cat("=" , rep("=", 60), "\n", sep = "")
cat("Aggregating Monte Carlo Results\n")
cat("=" , rep("=", 60), "\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

# ============================================================
# True values (computed via Monte Carlo with large N)
# ============================================================
invlogit <- function(x) 1 / (1 + exp(-x))

# Modified parameters: shift intercepts to create more extreme probabilities
# where logistic curves are more nonlinear (far from 0.5)
exp1_params <- list(
  a0 = -0.2, a1 = 0.7,                           # Treatment model (unchanged)
  b10 = -1.5, b11 = 1.2, b12 = 0.8,              # M1: low baseline, strong D effect
  b20 = 1.2, b21 = 0.8, b22 = 1.5, b23 = 0.6,   # M2: high baseline, strong M1 effect
  g0 = 0.0, g1 = 0.5, g2 = 0.6, g3 = 0.7, g4 = 0.4,  # Outcome model (unchanged)
  sigma = 1.0
)

# Compute true PSE values
compute_true_pse <- function(N_truth = 1000000L, par = exp1_params) {
  set.seed(99999)
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
    DY = mean(Y_1_00 - Y_0_00),
    DM2Y = mean(Y_1_01 - Y_1_00),
    DM1Y = mean(Y_1_11 - Y_1_01)
  )
}

# Compute true interventional effects
compute_true_intv <- function(N_truth = 1000000L, par = exp1_params) {
  set.seed(99999)
  C <- rnorm(N_truth)
  pM1 <- function(d, c) invlogit(par$b10 + par$b11 * d + par$b12 * c)
  pM2 <- function(d, m1, c) invlogit(par$b20 + par$b21 * d + par$b22 * m1 + par$b23 * c)
  EY <- function(d, m1, m2, c) par$g0 + par$g1 * d + par$g2 * m1 + par$g3 * m2 + par$g4 * c

  M1_0 <- rbinom(N_truth, 1, pM1(0, C))
  M1_1 <- rbinom(N_truth, 1, pM1(1, C))
  M2_1_natural <- rbinom(N_truth, 1, pM2(1, M1_1, C))
  M2_0_natural <- rbinom(N_truth, 1, pM2(0, M1_0, C))

  Y_1_natural <- EY(1, M1_1, M2_1_natural, C)
  Y_0_natural <- EY(0, M1_0, M2_0_natural, C)

  M1_for_M2star <- rbinom(N_truth, 1, pM1(0, C))
  M2_star_0 <- rbinom(N_truth, 1, pM2(0, M1_for_M2star, C))
  Y_1_M2star0 <- EY(1, M1_1, M2_star_0, C)

  c(
    OE = mean(Y_1_natural - Y_0_natural),
    IDE = mean(Y_1_M2star0 - Y_0_natural),
    IIE = mean(Y_1_natural - Y_1_M2star0)
  )
}

cat("Computing true values (1M samples)...\n")
true_pse <- compute_true_pse()
true_intv <- compute_true_intv()

cat("\nTrue Path-Specific Effects:\n")
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
cat("PATH-SPECIFIC EFFECTS RESULTS\n")
cat("=" , rep("=", 60), "\n\n")

# Define estimators and estimands for PSE
pse_estimators <- c("lin", "ipw", "reg", "med")
pse_estimands <- c("ATE", "DY", "DM2Y", "DM1Y")
pse_estimand_labels <- c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y")

# Compute statistics for each estimator and estimand
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

# Print PSE results
print(pse_stats, row.names = FALSE)

# Summary table by estimand
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

# Define estimators and estimands for interventional effects
intv_estimators <- c("ipw_intv", "rwr_intv", "med_intv")
intv_estimator_labels <- c("IPW", "RWR", "MedSim")
intv_estimands <- c("OE", "IDE", "IIE")

# Compute statistics for each estimator and estimand
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

# Print interventional results
print(intv_stats, row.names = FALSE)

# Summary table by estimand
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
summary_file <- file.path(output_dir, "mc_summary_exp1.rds")
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

# Also save as CSV for easy viewing
write.csv(pse_stats, file.path(output_dir, "pse_stats_exp1.csv"), row.names = FALSE)
write.csv(intv_stats, file.path(output_dir, "intv_stats_exp1.csv"), row.names = FALSE)
cat(sprintf("CSV files saved to: %s\n", output_dir))

cat("\n")
cat("=" , rep("=", 60), "\n", sep = "")
cat("AGGREGATION COMPLETE\n")
cat("=" , rep("=", 60), "\n")
