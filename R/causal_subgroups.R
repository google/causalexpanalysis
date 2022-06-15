# Copyright 2022 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Subgroup Analysis
#'
#' Perform a subgroup analysis for a single metric over different slices defined
#' by categorical covariates, with the following two cases:
#' (1) Estimate mean responses for each slice, if a treatment variable is not
#'     provided.
#' (2) Estimate average treatment effect (ATE) for each slice, if a treatment
#'     variable is provided.
#'
#' Optionally takes the model-estimated predictions of individual responses or
#' individual treatment effects (ITEs) to perform the estimation of the mean
#' responses or ATEs for each slice. If a model is not provided, we provide
#' the sample mean for (1) and the sample difference between Treatment and
#' Control groups for (2). We will refer to subgroups as slices.
#'
#' @param data A data frame containing the response metric and the slicing
#'   covariates, and optionally, the treatment variable, if the ATE is
#'   estimated.
#' @param y Names of metric.
#' @param model An optional model object for predicting the metric for (1) or
#'   estimating the ATE of the metric for (2). Currently, only the
#'   \code{causal_forest} model object of the grf R package is supported.
#' @param vars A character vector containing the names of the slicing
#'   covariates.
#' @param treat Name of treatment assignment indicator variable.
#' @param estimand The causal estimand to be estimated by the causal model,
#'   either \code{"ate"} for the average treatment effect on
#'   the entire population (ATE), \code{"att"} for the average treatment effect
#'   on the treated (ATT), or \code{"atu"} or \code{"atc"} for the average
#'   treatment effect on the untreated or control (ATU or ATC), or \code{"ato"}
#'   for the average treatment effect on the overlap population. Only applicable
#'   when the \code{treat} argument is not NULL.
#' @estimator The estimator to use for averaging over the ITEs to form an ATE
#'   estimate, whether \code{"aipw"} for the augmented inverse probability
#'   weighting (AIPW) method or \code{"tmle"} for the targeted minimum
#'   loss-based estimation (TMLE) adjustment method.
#' @param sep The separator character to separate the name of the slicing
#'   variable and its slice value.
#' @param label The name of column in the output table for the slices.
#' @param include_all Whether to include the slice containing all the
#'   observations, sometimes referred to as the global slice.
#' @param conf_level Confidence level for inference.
#' @param slice_name The name of the column containing the slice names.
#' @param percent_vars_resp The variables to percentize in the output table
#'   if we are estimating the mean response.
#' @param percent_vars_ate The variables to percentize in the output table
#'   if we are estimating the ATE.
#' @param digits Number of decimals for table statistics.
#' @return A data frame containing the table of estimation statistics for each
#'   slice.

SubgroupAnalysis <- function(data, y, model = NULL, vars = NULL, treat = NULL,
                             estimand = c("ate", "att", "atc", "atu", "ato"),
                             estimator = c("aipw", "tmle"), sep = "",
                             label = "slice", include_all = TRUE,
                             conf_level = 0.95, slice_name = "slice",
                             percent_vars_resp = c(
                               "mean", "est", "lb", "ub", "se"
                             ), percent_vars_ate = c(
                               "mC", "mT", "diff", "eff", "lb", "ub", "se",
                               "smd"
                             ), digits = digits) {
  y <- y[1]
  if (inherits(model, "causal_forest")) {
    estimand <- match.arg(estimand)
    estimator <- match.arg(estimator)
    target_sample <- switch(estimand,
      ate = "all",
      att = "treated",
      atc = "control",
      atu = "control",
      ato = "overlap"
    )
    method <- toupper(estimator)
  }
  estimation_stats <- list()
  if (include_all) {
    if (!is.null(treat) && inherits(model, "causal_forest")) {
      est_stats <- average_treatment_effect(
        forest = model, target.sample = target_sample, method = method
      )
      est <- est_stats[1]
      se <- est_stats[2]
    } else {
      est <- NULL
      se <- NULL
    }
    estimation_stats$all <- EstimationStats(
      data = data, y = y, treat = treat, est = est, se = se,
      conf_level = conf_level
    )
  }
  for (v in vars) {
    unique_values <- unique(data[[v]])
    for (uv in unique_values) {
      if (is.null(treat) && inherits(model, "causal_forest")) {
        est_stats <- average_treatment_effect(
          forest = model, subset = data[[v]] == uv,
          target.sample = target_sample, method = method
        )
        est <- est_stats[1]
        se <- est_stats[2]
      } else {
        est <- NULL
        se <- NULL
      }
      estimation_stats[[paste(v, uv, sep = sep)]] <- EstimationStats(
        data = data[data[[v]] == uv, ], y = y, treat = treat, est = est,
        se = se, conf_level = conf_level
      )
    }
  }
  if (is.null(treat)) {
    percent_vars <- percent_vars_resp
  } else {
    percent_vars <- percent_vars_ate
  }
  out <- Unsplit(estimation_stats, label = slice_name) %>%
    Percentize(vars = percent_vars, digits = digits)
  return(out)
}

#' Subgroup Analysis over Metrics
#'
#' A wrapper function to run the \code{SubgroupAnalysis} function for multiple
#' metrics. See documentation for the \code{SubgroupAnalysis} function for
#' further details.
#'
#' @param data A data frame containing the response metric and the slicing
#'   covariates, and optionally, the treatment variable, if the ATE is
#'   estimated.
#' @param y A character vector of names of metric.
#' @param models Either a single model object or a list of model objects.
#' @param treat Name of treatment assignment indicator variable.
#' @param percent_vars_resp The variables to percentize in the output table
#'   if we are estimating the mean response.
#' @param percent_vars_ate The variables to percentize in the output table
#'   if we are estimating the ATE.
#' @param digits Number of decimals for table statistics.
#' @return A data frame containing table of estimation statistics for each
#'   slice if there is only one metric, and a list of such data frames if
#'   there is more than one metric.

SubgroupAnalysisOverMetrics <- function(data, y, models = NULL, treat = NULL,
                                        percent_vars_resp = c(
                                          "mean", "est", "lb", "ub", "se"
                                        ), percent_vars_ate = c(
                                          "mC", "mT", "diff", "eff", "lb", "ub",
                                          "se", "smd"
                                        ), digits = 2, ...) {
  if (is.null(treat)) {
    percent_vars <- percent_vars_resp
  } else {
    percent_vars <- percent_vars_ate
  }
  if (length(y) > 1) {
    out <- list()
    for (i in seq_along(y)) {
      out[[y[i]]] <- SubgroupAnalysis(
        data = data, y = y[i], model = models[i], treat = treat, ...
      )
    }
    out <- out %>%
      MyApply(function(df) Percentize(df, x = percent_vars, digits = digits))
  } else {
    out <- SubgroupAnalysis(data, y = y, model = models, treat = treat, ...) %>%
      Percentize(x = percent_vars, digits = digits)
  }
  return(out)
}

#' Estimation Statistics
#'
#' A helper for the \code{SubgroupAnalysis} function to compute the estimation
#' statistics for each slice. Returns the table as a single-row data frame, with
#' the following statistics, depending on whether the mean response or the
#' average treatment effect (ATE) is estimated:
#' (1) Mean Response: Observation count \code{n}, the sample mean \code{mean},
#'     the estimated mean \code{est}, the confidence interval (CI) lower and
#'     upper bounds \code{lb} and \code{ub}, and the standard error (SE) of the
#'     estimate.
#' (2) Average Treatment Effect: Observation counts for Control and Treatment
#'     groups \code{nC} and \code{nT}, the sample means for the Control and
#'     Treatment groups \code{mC} and \code{mT}, the sample mean difference
#'     between the Treatment and Control groups \code{diff}, the estimated
#'     average treatment effect \code{eff}, its confidence interval (CI) lower
#'     and upper bounds \code{lb} and \code{ub}, and the standard error (SE) of
#'     the treatment effect estimate.
#'
#' @param data A data frame containing the response metric and the slicing
#'   covariates, and optionally, the treatment variable, if the ATE is
#'   estimated.
#' @param y A character vector of names of metric.
#' @param treat Name of treatment assignment indicator variable.
#' @param est The value of an external model estimate.
#' @param se The value of the standard error (SE) of the external model
#'   estimate.
#' @param conf_level Confidence level for inference.
#' @param resp_name The name of the estimated mean response in the output table.
#' @param eff_name The name of the estimated ATE in the output table.
#' @return A single-row data frame containing the estimation statistics.

EstimationStats <- function(data, y, treat = NULL, est = NULL, se = NULL,
                            conf_level = 0.95, resp_name = "est",
                            eff_name = "eff") {
  z_half_alpha <- qnorm((1 + conf_level) / 2)
  # If no treatment variable, compute average response. Otherwise, compute
  # the average treatment effect.
  if (is.null(treat)) {
    out <- data.frame(
      n = nrow(data),
      mean = mean(data[[y]]),
      sd = sd(data[[y]])
    )
    if (is.null(est)) {
      out <- out %>%
        mutate(est = mean, se = sd / sqrt(n))
    } else {
      out <- out %>%
        mutate(est = est, se = se)
    }
  } else {
    xC <- data[as.logical(1 - data[[treat]]), y]
    xT <- data[as.logical(data[[treat]]), y]
    out <- data.frame(
      nC = length(xC),
      nT = length(xT),
      mC = mean(xC),
      mT = mean(xT),
      sC = sd(xC),
      sT = sd(xT)
    ) %>%
      mutate(
        diff = mT - mC,
        s_pooled = EffectSizeSD(sC, sT, nC, nT, method = "pooled")
      )
    if (is.null(est)) {
      out <- out %>%
        mutate(
          est = diff,
          se = s_pooled * sqrt(1 / nC + 1 / nT),
          smd = est / s_pooled
        )
    } else {
      out <- out %>%
        mutate(
          est = est,
          se = se,
          smd = est * sqrt(1 / nC + 1 / nT) / se
        )
    }
  }
  # Compute symmetric CIs
  out <- out %>%
    dplyr::mutate(
      lb = est - z_half_alpha * se,
      ub = est + z_half_alpha * se
    )
  # Select final variables
  if (is.null(treat)) {
    out <- out %>%
      dplyr::select(n, mean, est, lb, ub, se)
    names(out)[which(names(out) == "est")] <- resp_name
  } else {
    out <- out %>%
      dplyr::select(nC, nT, mC, mT, diff, est, lb, ub, se, smd)
    names(out)[which(names(out) == "est")] <- eff_name
  }
  return(out)
}
