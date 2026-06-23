#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript aggregate_results_exp5_variant.R <output_dir> <n_reps> <variant_name>")
}

output_dir <- args[1]
n_reps <- as.integer(args[2])
variant_name <- args[3]

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) normalizePath(sub("^--file=", "", script_arg[1])) else normalizePath("aggregate_results_exp5_variant.R")
script_dir <- dirname(script_path)

suppressPackageStartupMessages({
  library(dplyr)
})

source(file.path(script_dir, "exp5_medflow_utils.R"))

cat("=", rep("=", 60), "\n", sep = "")
cat(sprintf("Aggregating Experiment 5 Variant: %s\n", variant_name))
cat("=", rep("=", 60), "\n\n", sep = "")

results_df <- load_medflow_results(output_dir, max_rep_id = n_reps)
observed_reps <- sort(unique(results_df$rep_id))
missing_reps <- setdiff(seq_len(n_reps), observed_reps)

cat(sprintf("Observed medflow replications: %d (expected: %d)\n", nrow(results_df), n_reps))
if (length(missing_reps) > 0) {
  cat(sprintf("Missing replications: %d\n", length(missing_reps)))
}

true_pse <- compute_true_pse_exp4()
true_intv <- compute_true_intv_exp4()

pse_stats <- compute_medflow_pse_stats(results_df, run_label = variant_name, true_pse = true_pse)
intv_stats <- compute_medflow_intv_stats(results_df, run_label = variant_name, true_intv = true_intv)

cat("\nPath-Specific Effects:\n")
print(pse_stats, row.names = FALSE)

cat("\nInterventional Effects:\n")
print(intv_stats, row.names = FALSE)

summary_list <- list(
  variant_name = variant_name,
  output_dir = output_dir,
  n_reps_expected = n_reps,
  n_reps_observed = nrow(results_df),
  missing_reps = missing_reps,
  true_pse = true_pse,
  true_intv = true_intv,
  results_df = results_df,
  pse_stats = pse_stats,
  intv_stats = intv_stats
)

summary_file <- file.path(output_dir, sprintf("mc_summary_exp5_%s.rds", variant_name))
results_file <- file.path(output_dir, sprintf("medflow_results_exp5_%s.rds", variant_name))
pse_file <- file.path(output_dir, sprintf("pse_stats_exp5_%s.csv", variant_name))
intv_file <- file.path(output_dir, sprintf("intv_stats_exp5_%s.csv", variant_name))
missing_file <- file.path(output_dir, sprintf("missing_reps_exp5_%s.csv", variant_name))

saveRDS(summary_list, summary_file)
saveRDS(results_df, results_file)
write.csv(pse_stats, pse_file, row.names = FALSE)
write.csv(intv_stats, intv_file, row.names = FALSE)
write.csv(data.frame(rep_id = missing_reps), missing_file, row.names = FALSE)

cat(sprintf("\nSummary saved to: %s\n", summary_file))
cat(sprintf("Results saved to: %s\n", results_file))
cat(sprintf("CSV files saved to: %s\n", output_dir))
