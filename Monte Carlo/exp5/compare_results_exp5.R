#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript compare_results_exp5.R <exp4_output_dir> <exp5_root_dir> <n_reps> [result_prefix]")
}

exp4_output_dir <- args[1]
exp5_root_dir <- args[2]
n_reps <- as.integer(args[3])
result_prefix <- if (length(args) >= 4) args[4] else "medflow"

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) normalizePath(sub("^--file=", "", script_arg[1])) else normalizePath("compare_results_exp5.R")
script_dir <- dirname(script_path)

suppressPackageStartupMessages({
  library(dplyr)
})

source(file.path(script_dir, "exp5_medflow_utils.R"))

variant_table <- data.frame(
  variant = c("arch_up", "arch_down", "lr_up", "lr_down", "batch_up", "batch_down"),
  group = c("architecture", "architecture", "learning_rate", "learning_rate", "training_batch_size", "training_batch_size"),
  setting = c("up", "down", "up", "down", "up", "down"),
  stringsAsFactors = FALSE
)

compute_paired_delta <- function(default_df, variant_df, variant_name, group_name, setting_name, prefix_map) {
  default_keep <- c("rep_id", unname(prefix_map))
  variant_keep <- c("rep_id", unname(prefix_map))

  if (!all(default_keep %in% names(default_df)) || !all(variant_keep %in% names(variant_df))) {
    return(data.frame())
  }

  merged_df <- default_df |>
    dplyr::select(dplyr::all_of(default_keep)) |>
    dplyr::inner_join(
      variant_df |> dplyr::select(dplyr::all_of(variant_keep)),
      by = "rep_id",
      suffix = c("_default", "_variant")
    )

  if (nrow(merged_df) == 0) {
    return(data.frame())
  }

  dplyr::bind_rows(lapply(names(prefix_map), function(estimand_label) {
    base_col <- prefix_map[[estimand_label]]
    default_col <- paste0(base_col, "_default")
    variant_col <- paste0(base_col, "_variant")
    delta <- merged_df[[variant_col]] - merged_df[[default_col]]
    delta <- delta[!is.na(delta)]
    if (length(delta) == 0) {
      return(NULL)
    }

    data.frame(
      variant = variant_name,
      group = group_name,
      setting = setting_name,
      estimand = estimand_label,
      n_common = length(delta),
      mean_delta = round(mean(delta), 4),
      sd_delta = round(stats::sd(delta), 4),
      mean_abs_delta = round(mean(abs(delta)), 4),
      rmse_delta = round(sqrt(mean(delta^2)), 4),
      stringsAsFactors = FALSE
    )
  }))
}

default_df <- load_medflow_results(
  exp4_output_dir,
  max_rep_id = n_reps,
  result_prefix = result_prefix
)
if (nrow(default_df) == 0) {
  stop(sprintf(
    "No default results with prefix '%s' found in exp4 output directory.",
    result_prefix
  ))
}

default_run_label <- if (identical(result_prefix, "medflow")) {
  "default_exp4_first240"
} else {
  paste0("default_exp4_", sub("^medflow_", "", result_prefix), "_first240")
}

run_results <- list(
  list(
    run = default_run_label,
    variant = "default",
    group = "default",
    setting = "baseline",
    output_dir = exp4_output_dir,
    results_df = default_df
  )
)

for (i in seq_len(nrow(variant_table))) {
  variant_name <- variant_table$variant[i]
  variant_dir <- file.path(exp5_root_dir, variant_name)
  if (!dir.exists(variant_dir)) {
    next
  }

  variant_df <- load_medflow_results(
    variant_dir,
    max_rep_id = n_reps,
    result_prefix = result_prefix
  )
  if (nrow(variant_df) == 0) {
    next
  }

  run_results[[length(run_results) + 1]] <- list(
    run = variant_name,
    variant = variant_name,
    group = variant_table$group[i],
    setting = variant_table$setting[i],
    output_dir = variant_dir,
    results_df = variant_df
  )
}

if (length(run_results) == 1) {
  warning("Only default exp4 results were found. No exp5 variant directories with medflow results were loaded.")
}

true_pse <- compute_true_pse_exp4()
true_intv <- compute_true_intv_exp4()

pse_stats <- dplyr::bind_rows(lapply(run_results, function(run_info) {
  compute_medflow_pse_stats(run_info$results_df, run_label = run_info$run, true_pse = true_pse) |>
    dplyr::mutate(
      variant = run_info$variant,
      group = run_info$group,
      setting = run_info$setting
    )
}))

intv_stats <- dplyr::bind_rows(lapply(run_results, function(run_info) {
  compute_medflow_intv_stats(run_info$results_df, run_label = run_info$run, true_intv = true_intv) |>
    dplyr::mutate(
      variant = run_info$variant,
      group = run_info$group,
      setting = run_info$setting
    )
}))

pse_delta <- dplyr::bind_rows(lapply(run_results[-1], function(run_info) {
  compute_paired_delta(
    default_df = default_df,
    variant_df = run_info$results_df,
    variant_name = run_info$variant,
    group_name = run_info$group,
    setting_name = run_info$setting,
    prefix_map = c(
      "ATE" = "mf_ATE",
      "D->Y" = "mf_DY",
      "D->M2->Y" = "mf_DM2Y",
      "D->M1~>Y" = "mf_DM1Y"
    )
  )
}))

intv_delta <- dplyr::bind_rows(lapply(run_results[-1], function(run_info) {
  compute_paired_delta(
    default_df = default_df,
    variant_df = run_info$results_df,
    variant_name = run_info$variant,
    group_name = run_info$group,
    setting_name = run_info$setting,
    prefix_map = c(
      "OE" = "mf_intv_OE",
      "IDE" = "mf_intv_IDE",
      "IIE" = "mf_intv_IIE"
    )
  )
}))

cat("Path-Specific Effect Comparison:\n")
print(pse_stats, row.names = FALSE)

cat("\nInterventional Effect Comparison:\n")
print(intv_stats, row.names = FALSE)

cat("\nPaired Delta vs default (PSE):\n")
print(pse_delta, row.names = FALSE)

cat("\nPaired Delta vs default (Interventional):\n")
print(intv_delta, row.names = FALSE)

summary_list <- list(
  exp4_output_dir = exp4_output_dir,
  exp5_root_dir = exp5_root_dir,
  n_reps = n_reps,
  result_prefix = result_prefix,
  true_pse = true_pse,
  true_intv = true_intv,
  run_results = run_results,
  pse_stats = pse_stats,
  intv_stats = intv_stats,
  pse_delta = pse_delta,
  intv_delta = intv_delta
)

saveRDS(summary_list, file.path(exp5_root_dir, "exp5_hyperparam_comparison.rds"))
write.csv(pse_stats, file.path(exp5_root_dir, "exp5_hyperparam_pse_stats.csv"), row.names = FALSE)
write.csv(intv_stats, file.path(exp5_root_dir, "exp5_hyperparam_intv_stats.csv"), row.names = FALSE)
write.csv(pse_delta, file.path(exp5_root_dir, "exp5_hyperparam_pse_delta_vs_default.csv"), row.names = FALSE)
write.csv(intv_delta, file.path(exp5_root_dir, "exp5_hyperparam_intv_delta_vs_default.csv"), row.names = FALSE)
