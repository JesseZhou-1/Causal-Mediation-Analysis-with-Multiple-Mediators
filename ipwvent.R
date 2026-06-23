#' Helper: extract P(D=d_val|...) from polr predicted probabilities
#'
#' Defined here for standalone use; identical to the version in ipwmed.R.
#' If both files are sourced, the later definition simply overwrites the earlier.
#'
#' @param polr_model A fitted polr model object.
#' @param newdata A data frame for prediction.
#' @param d_val The treatment level whose probability to extract.
#'
#' @return A numeric vector of predicted probabilities.
#' @noRd
polr_prob <- function(polr_model, newdata, d_val) {
  probs <- predict(polr_model, newdata = newdata, type = "probs")
  col_name <- as.character(d_val)
  if (col_name %in% colnames(probs)) return(probs[, col_name])
  stop(paste("Level", d_val, "not found in polr predictions"))
}


#' Helper: extract P(variable = level_val | ...) for all observations
#'
#' Works with both glm (binary) and polr (ordinal) models.
#'
#' @param model A fitted glm or polr model object.
#' @param newdata A data frame for prediction.
#' @param level_val The level whose probability to extract.
#'
#' @return A numeric vector of predicted probabilities.
#' @noRd
get_level_prob <- function(model, newdata, level_val) {
  if (inherits(model, "glm")) {
    prob1 <- predict(model, newdata = newdata, type = "response")
    if (level_val == 1) prob1 else 1 - prob1
  } else {
    polr_prob(model, newdata, level_val)
  }
}


#' Helper: extract P(variable = obs_values[i] | ...) for each observation i
#'
#' Works with both glm (binary) and polr (ordinal) models.
#'
#' @param model A fitted glm or polr model object.
#' @param newdata A data frame for prediction.
#' @param obs_values A numeric vector of observed values (one per observation).
#'
#' @return A numeric vector of predicted probabilities for the observed values.
#' @noRd
get_obs_prob <- function(model, newdata, obs_values) {
  if (inherits(model, "glm")) {
    prob1 <- predict(model, newdata = newdata, type = "response")
    ifelse(obs_values == 1, prob1, 1 - prob1)
  } else {
    probs <- predict(model, newdata = newdata, type = "probs")
    sapply(seq_along(obs_values), function(i) {
      col <- as.character(obs_values[i])
      if (col %in% colnames(probs)) probs[i, col] else 0
    })
  }
}


#' Inverse probability weighting (IPW) estimator for interventional effects:
#' inner function
#'
#' @description
#' Estimates the overall effect (OE), interventional direct effect (IDE), and
#' interventional indirect effect (IIE) using IPW with a treatment-induced
#' confounder L. Supports both binary (D in {0,1}) and ordinal/multi-valued
#' treatments via auto-detection: GLM for binary variables, polr for ordinal.
#'
#' @param data A data frame.
#' @param D Character scalar: name of the treatment variable.
#' @param M Character scalar: name of the focal mediator.
#' @param Y Character scalar: name of the outcome variable.
#' @param L Character scalar: name of the treatment-induced confounder.
#' @param D_formula A formula for the treatment model (D ~ C).
#' @param L_formula A formula for the confounder model (L ~ D + C).
#' @param M_formula A formula for the mediator model (M ~ D + L + C).
#' @param d Numeric scalar: treatment value of interest (default 1).
#' @param dstar Numeric scalar: reference treatment value (default 0).
#' @param stabilize Logical: stabilize the IPW weights (default TRUE).
#' @param censor Logical: censor the IPW weights (default TRUE).
#' @param censor_low Numeric: lower censoring quantile (default 0.01).
#' @param censor_high Numeric: upper censoring quantile (default 0.99).
#' @param minimal Logical: return only effect estimates (default FALSE).
#'
#' @return A list with OE, IDE, IIE, and optionally weights and model objects.
#' @noRd
ipwvent_inner <- function(
    data,
    D,
    M,
    Y,
    L,
    D_formula,
    L_formula,
    M_formula,
    d = 1,
    dstar = 0,
    stabilize = TRUE,
    censor = TRUE,
    censor_low = 0.01,
    censor_high = 0.99,
    minimal = FALSE
) {
  df <- data
  n <- nrow(df)

  # store numeric values before any factor conversion
  D_numeric <- as.numeric(as.character(df[[D]]))
  L_numeric <- as.numeric(as.character(df[[L]]))
  M_numeric <- as.numeric(as.character(df[[M]]))
  Y_numeric <- as.numeric(as.character(df[[Y]]))

  # auto-detect binary vs ordinal for each variable
  D_vals <- sort(unique(D_numeric))
  L_vals <- sort(unique(L_numeric))
  M_vals <- sort(unique(M_numeric))

  is_binary_D <- length(D_vals) == 2 && all(D_vals %in% c(0, 1))
  is_binary_L <- length(L_vals) == 2 && all(L_vals %in% c(0, 1))
  is_binary_M <- length(M_vals) == 2 && all(M_vals %in% c(0, 1))

  # fit treatment model
  if (is_binary_D) {
    D_model <- glm(D_formula, data = df, family = binomial(link = "logit"))
  } else {
    df[[D]] <- factor(df[[D]], ordered = TRUE)
    D_model <- MASS::polr(D_formula, data = df, method = "logistic", Hess = FALSE)
  }

  # fit confounder model
  if (is_binary_L) {
    L_model <- glm(L_formula, data = df, family = binomial(link = "logit"))
  } else {
    df[[L]] <- factor(df[[L]], ordered = TRUE)
    L_model <- MASS::polr(L_formula, data = df, method = "logistic", Hess = FALSE)
  }

  # fit mediator model
  if (is_binary_M) {
    M_model <- glm(M_formula, data = df, family = binomial(link = "logit"))
  } else {
    df[[M]] <- factor(df[[M]], ordered = TRUE)
    M_model <- MASS::polr(M_formula, data = df, method = "logistic", Hess = FALSE)
  }

  # P(D|C) for d and dstar
  pD_d_C <- get_level_prob(D_model, df, d)
  pD_dstar_C <- get_level_prob(D_model, df, dstar)

  # marginal P(D)
  pD_d <- mean(D_numeric == d)
  pD_dstar <- mean(D_numeric == dstar)

  # marginalized P(M_obs|D=d_val,C) = sum_l P(M_obs|D=d_val,L=l,C) * P(L=l|D=d_val,C)
  L_levels <- L_vals

  compute_marg_pM <- function(d_val) {
    temp <- df
    if (is_binary_D) {
      temp[[D]] <- d_val
    } else {
      temp[[D]] <- factor(d_val, levels = levels(df[[D]]), ordered = TRUE)
    }

    marg_pM <- rep(0, n)
    for (l_val in L_levels) {
      temp_l <- temp
      if (is_binary_L) {
        temp_l[[L]] <- l_val
      } else {
        temp_l[[L]] <- factor(l_val, levels = levels(df[[L]]), ordered = TRUE)
      }

      # P(L=l|D=d_val,C)
      pL_l <- get_level_prob(L_model, temp, l_val)

      # P(M_obs|D=d_val,L=l,C)
      pM_obs <- get_obs_prob(M_model, temp_l, M_numeric)

      marg_pM <- marg_pM + pM_obs * pL_l
    }
    marg_pM
  }

  marg_pM_d <- compute_marg_pM(d)
  marg_pM_dstar <- compute_marg_pM(dstar)

  # P(M_obs|D_obs,L_obs,C) — fully observed
  pM_obs_DLC <- get_obs_prob(M_model, df, M_numeric)
  pM_obs_DLC <- pmax(pM_obs_DLC, 1e-10) # floor to avoid division by zero

  # subgroups
  group_d <- D_numeric == d
  group_dstar <- D_numeric == dstar

  # weights
  # w1: E[Y(dstar,M(dstar))] — use D=dstar units
  w1 <- ifelse(group_dstar, marg_pM_dstar / (pD_dstar_C * pM_obs_DLC), 0)
  # w2: E[Y(d,M(d))] — use D=d units
  w2 <- ifelse(group_d, marg_pM_d / (pD_d_C * pM_obs_DLC), 0)
  # w3: E[Y(d,M*(dstar))] — use D=d units, reweight M to dstar marginal
  w3 <- ifelse(group_d, marg_pM_dstar / (pD_d_C * pM_obs_DLC), 0)

  # stabilize
  if (stabilize) {
    w1 <- w1 * pD_dstar
    w2 <- w2 * pD_d
    w3 <- w3 * pD_d
  }

  # censor
  if (censor) {
    w1[group_dstar] <- trimQ(w1[group_dstar], low = censor_low, high = censor_high)
    w2[group_d] <- trimQ(w2[group_d], low = censor_low, high = censor_high)
    w3[group_d] <- trimQ(w3[group_d], low = censor_low, high = censor_high)
  }

  # estimate effects
  Ehat_YdstarMdstar <- weighted.mean(Y_numeric, w1)
  Ehat_YdMd <- weighted.mean(Y_numeric, w2)
  Ehat_YdMdstar <- weighted.mean(Y_numeric, w3)

  OE  <- Ehat_YdMd - Ehat_YdstarMdstar
  IDE <- Ehat_YdMdstar - Ehat_YdstarMdstar
  IIE <- Ehat_YdMd - Ehat_YdMdstar

  if (minimal) {
    out <- list(OE = OE, IDE = IDE, IIE = IIE)
  } else {
    out <- list(
      OE = OE, IDE = IDE, IIE = IIE,
      weights1 = w1, weights2 = w2, weights3 = w3,
      model_D = D_model, model_L = L_model, model_M = M_model
    )
  }
  return(out)
}
