## Legacy exp4 truth helper for exp5 comparisons.
## Matches the original exp4 contrast d = 2 vs dstar = 0.
## Note: the Monte Carlo truth does not depend on whether the finite-sample run
## used n = 10000 or n = 20000; n only affects the realized experiment, not the
## population truth under the DGP.

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
  ) |>
    dplyr::arrange(.data$rep_id, dplyr::desc(.data$mtime), dplyr::desc(.data$in_rep_results)) |>
    dplyr::distinct(.data$rep_id, .keep_all = TRUE) |>
    dplyr::arrange(.data$rep_id)

  file_tbl$path
}

load_medflow_results <- function(output_dir, max_rep_id = Inf) {
  medflow_files <- find_medflow_result_files(output_dir)
  if (length(medflow_files) == 0) {
    return(data.frame(rep_id = integer()))
  }

  medflow_list <- lapply(medflow_files, function(path) {
    tryCatch(read.csv(path), error = function(e) NULL)
  })
  medflow_list <- medflow_list[!vapply(medflow_list, is.null, logical(1))]

  if (length(medflow_list) == 0) {
    return(data.frame(rep_id = integer()))
  }

  results_df <- dplyr::bind_rows(medflow_list)
  if (!("rep_id" %in% names(results_df))) {
    return(data.frame(rep_id = integer()))
  }

  results_df <- results_df |>
    dplyr::arrange(.data$rep_id) |>
    dplyr::distinct(.data$rep_id, .keep_all = TRUE)

  if (is.finite(max_rep_id)) {
    results_df <- results_df |>
      dplyr::filter(.data$rep_id <= max_rep_id)
  }

  results_df |>
    dplyr::arrange(.data$rep_id)
}

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

compute_true_pse_exp4 <- function(N_truth = 2000000L, par = exp4_params, d = 2, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)

  M1_d <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * d +
    par$M1_beta_C * C + par$M1_beta_DC * d * C))
  M1_dstar <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
    par$M1_beta_C * C + par$M1_beta_DC * dstar * C))

  M2_d_M1d <- rnorm(N_truth, mu_M2_fn(d, M1_d, C, par), par$M2_sigma)
  M2_dstar_M1ds <- rnorm(N_truth, mu_M2_fn(dstar, M1_dstar, C, par), par$M2_sigma)
  M2_d_M1ds <- rnorm(N_truth, mu_M2_fn(d, M1_dstar, C, par), par$M2_sigma)

  Y_d_M1d_M2d <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_d_M1d, C, par), par$Y_cuts)
  Y_dstar <- rordlogit(N_truth, lp_Y_fn(dstar, M1_dstar, M2_dstar_M1ds, C, par), par$Y_cuts)
  Y_d_M1ds_M2ds <- rordlogit(N_truth, lp_Y_fn(d, M1_dstar, M2_dstar_M1ds, C, par), par$Y_cuts)
  Y_d_M1ds_M2d <- rordlogit(N_truth, lp_Y_fn(d, M1_dstar, M2_d_M1ds, C, par), par$Y_cuts)

  c(
    ATE = mean(Y_d_M1d_M2d) - mean(Y_dstar),
    DY = mean(Y_d_M1ds_M2ds) - mean(Y_dstar),
    DM2Y = mean(Y_d_M1ds_M2d) - mean(Y_d_M1ds_M2ds),
    DM1Y = mean(Y_d_M1d_M2d) - mean(Y_d_M1ds_M2d)
  )
}

compute_true_intv_exp4 <- function(N_truth = 2000000L, par = exp4_params, d = 2, dstar = 0) {
  set.seed(99999)
  C <- rnorm(N_truth)

  M1_d <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * d +
    par$M1_beta_C * C + par$M1_beta_DC * d * C))
  M1_dstar <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
    par$M1_beta_C * C + par$M1_beta_DC * dstar * C))

  M2_d_natural <- rnorm(N_truth, mu_M2_fn(d, M1_d, C, par), par$M2_sigma)
  M2_dstar_natural <- rnorm(N_truth, mu_M2_fn(dstar, M1_dstar, C, par), par$M2_sigma)

  Y_d_natural <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_d_natural, C, par), par$Y_cuts)
  Y_dstar_natural <- rordlogit(N_truth, lp_Y_fn(dstar, M1_dstar, M2_dstar_natural, C, par), par$Y_cuts)

  M1_for_M2star <- rbinom(N_truth, 1, plogis(par$M1_beta0 + par$M1_beta_D * dstar +
    par$M1_beta_C * C + par$M1_beta_DC * dstar * C))
  M2_star_dstar <- rnorm(N_truth, mu_M2_fn(dstar, M1_for_M2star, C, par), par$M2_sigma)

  Y_d_M2star <- rordlogit(N_truth, lp_Y_fn(d, M1_d, M2_star_dstar, C, par), par$Y_cuts)

  c(
    OE = mean(Y_d_natural) - mean(Y_dstar_natural),
    IDE = mean(Y_d_M2star) - mean(Y_dstar_natural),
    IIE = mean(Y_d_natural) - mean(Y_d_M2star)
  )
}

compute_medflow_pse_stats <- function(results_df, run_label, true_pse = compute_true_pse_exp4()) {
  estimands <- c("ATE", "DY", "DM2Y", "DM1Y")
  estimand_labels <- c("ATE", "D->Y", "D->M2->Y", "D->M1~>Y")

  stats_list <- lapply(seq_along(estimands), function(i) {
    estimand <- estimands[i]
    col_name <- paste0("mf_", estimand)
    if (!(col_name %in% names(results_df))) {
      return(NULL)
    }

    estimates <- results_df[[col_name]]
    n_valid <- sum(!is.na(estimates))
    if (n_valid == 0) {
      return(NULL)
    }

    true_val <- true_pse[[estimand]]
    data.frame(
      Run = run_label,
      Estimand = estimand_labels[i],
      Truth = round(true_val, 4),
      Mean = round(mean(estimates, na.rm = TRUE), 4),
      Bias = round(mean(estimates, na.rm = TRUE) - true_val, 4),
      Variance = round(var(estimates, na.rm = TRUE), 4),
      RMSE = round(sqrt(mean((estimates - true_val)^2, na.rm = TRUE)), 4),
      N_valid = n_valid,
      stringsAsFactors = FALSE
    )
  })

  out <- dplyr::bind_rows(stats_list)

  if (all(c("mf_DM1Y", "mf_DM2Y") %in% names(results_df))) {
    mnie_estimates <- ifelse(
      is.na(results_df$mf_DM1Y) | is.na(results_df$mf_DM2Y),
      NA_real_,
      results_df$mf_DM1Y + results_df$mf_DM2Y
    )
    n_valid <- sum(!is.na(mnie_estimates))
    if (n_valid > 0) {
      true_mnie <- true_pse[["DM1Y"]] + true_pse[["DM2Y"]]
      out <- dplyr::bind_rows(
        out,
        data.frame(
          Run = run_label,
          Estimand = "MNIE",
          Truth = round(true_mnie, 4),
          Mean = round(mean(mnie_estimates, na.rm = TRUE), 4),
          Bias = round(mean(mnie_estimates, na.rm = TRUE) - true_mnie, 4),
          Variance = round(var(mnie_estimates, na.rm = TRUE), 4),
          RMSE = round(sqrt(mean((mnie_estimates - true_mnie)^2, na.rm = TRUE)), 4),
          N_valid = n_valid,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  out
}

compute_medflow_intv_stats <- function(results_df, run_label, true_intv = compute_true_intv_exp4()) {
  estimands <- c("OE", "IDE", "IIE")

  stats_list <- lapply(estimands, function(estimand) {
    col_name <- paste0("mf_intv_", estimand)
    if (!(col_name %in% names(results_df))) {
      return(NULL)
    }

    estimates <- results_df[[col_name]]
    n_valid <- sum(!is.na(estimates))
    if (n_valid == 0) {
      return(NULL)
    }

    true_val <- true_intv[[estimand]]
    data.frame(
      Run = run_label,
      Estimand = estimand,
      Truth = round(true_val, 4),
      Mean = round(mean(estimates, na.rm = TRUE), 4),
      Bias = round(mean(estimates, na.rm = TRUE) - true_val, 4),
      Variance = round(var(estimates, na.rm = TRUE), 4),
      RMSE = round(sqrt(mean((estimates - true_val)^2, na.rm = TRUE)), 4),
      N_valid = n_valid,
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(stats_list)
}
