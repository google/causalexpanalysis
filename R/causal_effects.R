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

#' Covariate Balance Table
#'
#' Create a table to show covariate balance for a list of specified covariates
#' \code{x} using the standardized mean difference (SMD) and Welch's t-test
#' statistic for mean difference and its p-value, with options for paired and
#' unpaired scenarios. Note that the SMD is Cohen's d for measuring effect size.
#'
#' @param data A data frame containing covariates for treatment and control
#'   subjects with indicator variable for treatment assignment.
#' @param x Character vector containing names of covariates.
#' @param paired Whether the data is paired, in which case the SMD and t
#'   statistic will be computed for paired data.
#' @param treat Name of treatment assignment indicator variable.
#' @param xlab Name for table column containing covariates.
#' @param percent_vars Names for variables to format as percentages. Note that
#'   this only changes the printing format, not the actual value of the
#'   variables.
#' @param format_vars Names for variables to format as a float. Note that
#'   this only changes the printing format, not the actual value of the
#'   variables.
#' @param digits Number of decimals for table statistics.
#' @return A data frame containing the covariate balance table.

BalanceTable <- function(data, x = NULL, treat = "treat",
                         xlab = "covariates", paired = FALSE,
                         percent_vars = c("smd", "pval"), format_vars = "stat",
                         digits = 1) {
  stopifnot(length(x) > 0)
  trt <- data[[treat]]
  model_formula <- CreateFormula(treat, x)
  mat <- model.matrix(model_formula, data)
  vars <- setdiff(colnames(mat), "(Intercept)")
  out <- vars %>%
    MyApply(function(x) {
      xC <- mat[as.logical(1 - trt), x]
      xT <- mat[as.logical(trt), x]
      t_test <- t.test(xT, xC,
        paired = paired, var.equal = FALSE,
      )
      df <- data.frame(
        nC = length(xC),
        nT = length(xT),
        mC = mean(xC),
        mT = mean(xT),
        sC = sd(xC),
        sT = sd(xT)
      )
      if (paired) {
        sD <- sd(xT - xC)
        df <- df %>%
          mutate(smd = EffectSize(mT, mC, sD))
      } else {
        df <- df %>%
          mutate(smd = EffectSize(mT, mC, sT, sC))
      }
      df <- df %>%
        mutate(
          smd = ifelse(is.nan(smd), 0, smd),
          stat = ifelse(is.nan(t_test$statistic), 0, t_test$statistic),
          pval = ifelse(is.nan(t_test$p.value), 1, t_test$p.value)
        )
      return(df)
    }) %>%
    Unsplit(label = xlab) %>%
    data.frame()
  out <- FormatTable(
    data = out, format_var = format_vars, percent_vars = percent_vars,
    digits = digits
  )
  return(out)
}

#' Treatment Effect Table
#'
#' Create a table to show mean differences between the treatment and control
#' groups on any set of variables, which we will refer to as metrics. The
#' metrics can be experiment outcomes or pre-treatment covariates. In the former
#' case, the mean difference is the average treatment effect for the outcome;
#' in the latter case, the mean difference measures the covariate imbalance
#' between the two arms. User has the option to measure two types of effects:
#' (1) Absolute effect, where the t test is used to construct CIs and optionally
#'     show the inference statistics for the test, including t statistic and
#'     p-value, and the standardized mean difference with corresponding
#'     observation standard deviation, as a measure of effect size.
#' (2) Relative effect, where the CIs are constructed using either the Fieller
#'     or delta method.
#' The metrics passed into the function can include categorical variables,
#' either factors or character vectors, in which case, the variable is expanded
#' into a set of dummy variables using the default contrasts available in R.
#' There is also an option to adjust the SEs for paired data.
#'
#' @param data A data frame containing metrics for treatment and control
#'   subjects with indicator variable for treatment assignment. Binary metrics
#'   can be included as well.
#' @param metrics Character vector containing names of metrics.
#' @param paired Whether the data is paired, in which case the SMD and t
#'   statistic will be computed for paired data.
#' @param var_equal Whether to assume unequal variances for the metrics across
#'   treatment and control groups. Similar to the \code{var.equal} argument in
#'   the \code{stats::t.test} function. Applies only when \code{paired = FALSE}.
#' @param treat Name of treatment assignment indicator variable.
#' @param match Name for variable for matching paired subjects.
#' @param type The type of effect(s) to show in the table, whether an absolute
#'   effect, a relative effect, or both.
#' @param reff_ci_method A string for the relative effect CI method, either
#'   \code{"fieller"} for the Fieller CI, \code{"delta"} for the delta method
#'   without log transformation of the ratio, or \code{"delta_log"} for the
#'   delta method with log transformation of the ratio. See documentation for
#'   the \code{RatioCI} function for details.
#' @param conf_level Confidence level for inference.
#' @param sig_test Whether to show the inference statitics for a significance
#'   test, i.e., the standard error of the effect estimate, the t statistic, and
#'   the p-value.
#' @param eff_size Whether to show the statitics for measuring effect size,
#'   namely, the joint standard deviation for the two samples and the
#'   standardized mean difference (SMD), i.e., the mean difference divided by
#'   the joint standard deviation.
#' @param digits Number of decimals for table statistics.
#' @param match_counts Whether to show the matched counts \code{x00},\code{x01},
#'   \code{x10}, and \code{x11} in the table.
#' @return A data frame containing the effect table for binomial proportions.

EffectTable <- function(data, metrics, paired = FALSE, var_equal = FALSE,
                        treat = "treat", match = "match",
                        type = c("both", "abs", "rel"),
                        reff_ci_method = c("fieller", "delta", "delta_log"),
                        conf_level = 0.95, sig_test = TRUE, eff_size = TRUE,
                        digits = 1) {
  stopifnot(sum(data[[treat]]) > 0 & sum(1 - data[[treat]]) > 0)
  type <- match.arg(type)
  data2 <- CovariateData(data, vars = metrics, center = FALSE)
  metrics <- names(data2)
  data2$trt <- ifelse(data[[treat]], "mT", "mC")
  keys <- "trt"
  lhs <- "metric"
  if (paired) {
    data2$match <- data[[match]]
    keys <- c("match", keys)
    lhs <- c(lhs, "match")
  }
  dcast_formula <- CreateFormula(lhs, "trt")
  out0 <- data.frame(
    metric = metrics,
    nC = sum(1 - data[[treat]]),
    nT = sum(data[[treat]])
  )
  if (paired) {
    vars <- c("mC", "mT", "eff")
    out1 <- data2 %>%
      group_by_at(keys) %>%
      summarize_at(metrics, sum, na.rm = TRUE) %>%
      data.frame() %>%
      reshape2::melt(id.vars = keys) %>%
      dplyr::rename(metric = variable) %>%
      reshape2::dcast(dcast_formula) %>%
      mutate(eff = mT - mC)
    out2 <- out1 %>%
      group_by(metric) %>%
      summarize_at(vars, sd, na.rm = TRUE) %>%
      dplyr::rename(sC = mC, sT = mT, sD = eff)
    out1 <- out1 %>%
      group_by(metric) %>%
      summarize_at(vars, mean, na.rm = TRUE)
  } else {
    out1 <- data2 %>%
      group_by_at(keys) %>%
      summarize_at(metrics, mean, na.rm = TRUE) %>%
      data.frame() %>%
      reshape2::melt(id.vars = keys) %>%
      dplyr::rename(metric = variable) %>%
      reshape2::dcast(dcast_formula) %>%
      mutate(eff = mT - mC)
    out2 <- data2 %>%
      mutate(trt = ifelse(trt == "mT", "sT", "sC")) %>%
      group_by_at(keys) %>%
      summarize_at(metrics, sd, na.rm = TRUE) %>%
      data.frame() %>%
      reshape2::melt(id.vars = keys) %>%
      dplyr::rename(metric = variable) %>%
      reshape2::dcast(dcast_formula)
  }
  out <- out0 %>%
    left_join(out1, by = "metric") %>%
    left_join(out2, by = "metric") %>%
    mutate(
      mC_se = sC / sqrt(nC),
      mT_se = sT / sqrt(nT)
    )
  if (!paired) {
    out <- out %>%
      mutate(sD = sqrt(2) * EffectSizeSD(
        sC, sT, nC, nT,
        method = ifelse(var_equal, "pooled", "average")
      ))
  }
  out <- out %>%
    mutate(sd = sD / sqrt(2))
  if (paired) {
    out <- out %>%
      mutate(
        se = sD / sqrt(nT),
        df = nT - 1
      )
  } else {
    if (var_equal) {
      out <- out %>%
        mutate(
          se = sd * sqrt((1 / nC) + (1 / nT)),
          df = nC + nT - 2
        )
    } else {
      out <- out %>%
        mutate(
          se = sqrt(mC_se^2 + mT_se^2),
          df = WelchSatter(sC, sT, 1 / nT, 1 / nT, nC - 1, nT - 1)
        )
    }
  }
  if (type %in% c("abs", "both")) {
    out <- out %>%
      mutate(
        t_half_alpha = qt((1 + conf_level) / 2, df = df),
        lb = eff - t_half_alpha * se,
        ub = eff + t_half_alpha * se
      )
    if (sig_test) {
      out <- out %>%
        mutate(
          stat = ifelse(eff == 0, 0, eff / se),
          pval = ifelse(stat == 0, 1, 2 * pt(-abs(stat), df = df))
        )
    }
    if (eff_size) {
      out <- out %>%
        mutate(smd = eff / sd)
    }
  }
  if (type %in% c("rel", "both")) {
    reff_ci_method <- match.arg(reff_ci_method)
    out$reff <- with(out, eff / mC)
    if (paired) {
      cov_CT <- with(out, mT_se^2 + mC_se^2 - se^2) / 2
    } else {
      cov_CT <- 0
    }
    # Since the numerator a is the treatment response, we need to set chg = TRUE
    # to get CI bounds for the relative lift instead of the response ratio.
    ci_df <- with(out, RatioCI(
      method = reff_ci_method, a = mT, b = mC, v11 = mT_se^2, v22 = mC_se^2,
      v12 = cov_CT, chg = TRUE
    ))
    names(ci_df) <- paste0("re_", c("lb", "ub"))
    out <- out %>%
      bind_cols(ci_df)
  }
  # Format table
  table_vars <- c("metric", "nC", "nT", "mC", "mT")
  percent_vars <- c("mC", "mT")
  format_vars <- NULL
  if (type %in% c("abs", "both")) {
    eff_vars <- c("eff", "lb", "ub")
    table_vars <- c(table_vars, eff_vars)
    percent_vars <- c(percent_vars, eff_vars)
    if (sig_test) {
      sig_test_vars <- c("se", "stat", "pval")
      sig_test_percent_vars <- c("se", "pval")
      sig_test_format_vars <- setdiff(sig_test_vars, sig_test_percent_vars)
      table_vars <- c(table_vars, sig_test_vars)
      percent_vars <- c(percent_vars, sig_test_percent_vars)
      format_vars <- c(format_vars, sig_test_format_vars)
    }
    if (eff_size) {
      eff_size_vars <- c("sd", "smd")
      table_vars <- c(table_vars, eff_size_vars)
      percent_vars <- c(percent_vars, eff_size_vars)
    }
  }
  if (type %in% c("rel", "both")) {
    reff_vars <- c("reff", "re_lb", "re_ub")
    table_vars <- c(table_vars, reff_vars)
    percent_vars <- c(percent_vars, reff_vars)
  }
  out <- out %>%
    select_at(table_vars)
  out <- FormatTable(
    data = out, format_var = format_vars, percent_vars = percent_vars,
    digits = digits
  )
  return(out)
}

#' Treatment Effect Table for Binomial Proportions
#'
#' Create a table to show the average treatment effect for binomial proportions
#' of binary metrics, showing the estimated mean difference with standard error,
#' z statistic, and p-value, with options for paired and unpaired proportions
#'
#' @param data A data frame containing binary metrics for treatment and control
#'   subjects with indicator variable for treatment assignment.
#' @param metrics Character vector containing names of metrics.
#' @param paired Whether the data is paired, in which case the SMD and t
#'   statistic will be computed for paired data.
#' @param weight Name of variable containing subject weights.
#' @param treat Name of treatment assignment indicator variable.
#' @param match Name for variable for matching paired subjects.
#' @param type The type of effect(s) to show in the table, whether an absolute
#'   effect, a relative effect, or both.
#' @param conf_level Confidence level for inference.
#' @param sig_test Whether to show the inference statitics for a significance
#'   test, i.e., the standard error of the effect estimate, the t statistic, and
#'   the p-value.
#' @param match_counts Whether to show the matched counts \code{x00},\code{x01},
#'   \code{x10}, and \code{x11} in the table.
#' @return A data frame containing the effect table for binomial proportions.

EffectTableProp <- function(data, metrics, paired = FALSE, weight = "n",
                            treat = "treat", match = "match",
                            type = c("abs", "rel", "both"), conf_level = 0.95,
                            sig_test = FALSE, digits = 1,
                            match_counts = FALSE) {
  type <- match.arg(type)
  data$n <- data[[weight]]
  data$trt <- ifelse(data[[treat]], "xT", "xC")
  keys <- "trt"
  vars <- c("xC", "xT")
  lhs <- "metric"
  if (paired) {
    data$match <- data[[match]]
    keys <- c("match", keys)
    vars <- c(vars, paste0("x", c("00", "01", "10", "11")))
    lhs <- c(lhs, "match")
  }
  dcast_formula <- CreateFormula(lhs, "trt")
  out0 <- data.frame(
    metric = metrics,
    nC = sum((1 - data[[treat]]) * data[[weight]]),
    nT = sum(data[[treat]] * data[[weight]])
  )
  out1 <- data %>%
    group_by_at(keys) %>%
    summarize_at(metrics, sum, na.rm = TRUE) %>%
    data.frame() %>%
    melt(id.vars = keys) %>%
    dplyr::rename(metric = variable) %>%
    dcast(dcast_formula)
  if (paired) {
    out1 <- out1 %>%
      mutate(
        x00 = 1 * (xC == 0 & xT == 0),
        x01 = 1 * (xC == 0 & xT == 1),
        x10 = 1 * (xC == 1 & xT == 0),
        x11 = 1 * (xC == 1 & xT == 1)
      )
  }
  out1 <- out1 %>%
    group_by(metric) %>%
    summarize_at(vars, sum)
  out <- out0 %>%
    left_join(out1, by = "metric") %>%
    mutate(rC = xC / nC, rT = xT / nT)
  percent_vars <- c("rC", "rT")
  format_vars <- NULL
  if (type %in% c("abs", "both")) {
    out$eff <- with(out, rT - rC)
    if (paired) {
      out$se <- with(out, MatchedPropDiffSE(x01, x10, nT))
    } else {
      out$se <- with(out, PropDiffSE(rT, rC, nT, nC))
    }
    z_alpha <- qnorm((1 + conf_level) / 2)
    out <- out %>%
      mutate(
        stat = ifelse(eff == 0, 0,
          formattable::formattable(eff / se, digits = digits, format = "f")
        ),
        pval = ifelse(stat == 0, 1, 2 * pnorm(-abs(stat))),
        lb = eff - z_alpha * se,
        ub = eff + z_alpha * se
      )
    percent_vars <- c(percent_vars, "eff", "lb", "ub")
    if (sig_test) {
      percent_vars <- c(percent_vars, "se", "pval")
      format_vars <- "stat"
    } else {
      out <- out %>%
        dplyr::select(-se, -stat, -pval)
    }
  }
  if (type %in% c("rel", "both")) {
    out$reff <- with(out, rT / rC - 1)
    if (paired) {
      ci_df <- with(
        out,
        # Since the numerator a is the treatment response, we need to set
        # chg = TRUE to get CI bounds for the relative lift instead of the
        # response ratio.
        FiellerCIMatchedProp(
          x00, x01, x10, x11,
          conf_level = conf_level, chg = TRUE
        )
      )
    } else {
      ci_df <- with(
        out,
        FiellerCIProp(rT, rC, nT, nC, conf_level = conf_level, chg = TRUE)
      )
    }
    names(ci_df) <- paste0("re_", c("lb", "ub"))
    out <- bind_cols(out, ci_df)
    percent_vars <- c(percent_vars, "reff", "re_lb", "re_ub")
  }
  if (paired && !match_counts) {
    out <- out %>%
      dplyr::select(-x00, -x01, -x10, -x11)
  }
  out <- FormatTable(
    data = out, format_var = format_vars, percent_vars = percent_vars,
    digits = digits
  )
  return(out)
}

#' Table of Regression-adjusted Treatment Effects for Pre-post Analysis
#'
#' Estimate regression-adjusted treatment effects for multiple metrics using
#' a linear model, and summarize the effects and inference statistics in a
#' table. The metric is regressed on a treatment indicator variable and
#' potentially the following:
#' (1) Covariates not interacted with the treatment variable, specified in the
#'     \code{vars} argument, and
#' (2) Covariates interacted with the treatment variable, specified in the
#'     \code{int_vars} argument
#' All covariates are centered, i.e., have their means subtracted from their
#' values, so that:
#' (1) The intercept coefficient becomes the average control response, and
#' (2) The coefficient for treatment variable becomes the average treatment
#'     effect, i.e., the average treatment response minus the average control
#'     response.
#' Users can include the pre-period value of the metric in the model, in which
#' case a pre-post analysis is conducted. We will refer to the pre-period value
#' as the pre-period metric, and the experiment-period value as the post-period
#' metric. Whether the pre-period metric is included in the model or not, the
#' name of the post-period metric must always follow the convention
#' \code{paste(y, "_post")}, where \code{y} is the vector metric names without
#' the time suffixes of "_post" and "_pre".
#' Users can specify two model types through the \code{model} argument:
#' (1) ANCOVA (\code{ancova}), regression of metric on treatment variable and
#'     optionally other covariates, including the pre-period metric, or
#' (2) Difference-in-difference (\code{did}), regression of difference between
#'     post-period and pre-period metrics on the treatment variable and
#'     optionally other covariates.
#' The difference-in-difference model is generally fitted without additional
#' covariates, but we provide the option for the user to include them as well.
#' The ANCOVA model without interacted covariates is known as the equal-slopes
#' ANCOVA model, while that with interacted covariates is known as the
#' unequal-slopes ANCOVA model. We have Lin's (2013) model when the following
#' conditions hold:
#' (1) The covariates interacted with the treatment are the same as those not
#'     interacted with the treatment, and
#' (2) All covariates are centered.
#' That is, Lin's (2013) model is an unequal-slopes ANCOVA model with centered
#' covariates, all of which appear twice, both as a non-interaction term and as
#' a treatment-interaction term.
#'
#' @param data A data frame containing metrics for treatment and control
#'   subjects, an indicator variable for treatment assignment and optionally
#'   additional covariates and an integer variable for matching of pairs.
#' @param y Character vector containing names of metrics. These are names of the
#'   metrics without the \code{_post} and \code{_pre} suffixes.
#' @param vars Character vector containing the names of covariates to include
#'   in the regression model without interacting with the treatment variable.
#'   These do not need to include the pre-period values of the metrics specified
#'   in \code{y} if \code{include_pre = TRUE}.
#' @param int_vars Character vector containing the names of covariates to
#'   include in the linear model for interacting with the treatment variable.
#'   These do not need to include the pre-period values of the metrics specified
#'   in \code{y} if both \code{include_pre = TRUE} and
#'   \code{interact_pre = TRUE} hold.
#' @param treat Name of treatment assignment indicator variable.
#' @param unit Name of identifier variable for each study unit, e.g.,
#'   \code{obscured_gaia_id} for users in the context of Cloud experiments.
#'   Defaults to \code{NULL}. This is needed when each unit has potentially
#'   more than one measurement or observation, e.g., when \code{panel_time} is
#'   specified, in the case of longitudinal or panel data,
#' @param estimand The causal estimand to be estimated by the regression
#'   adjustment model, either \code{"ate"} for the average treatment effect on
#'   the entire population (ATE), \code{"att"} for the average treatment effect
#'   on the treated (ATT), or \code{"atu"} or \code{"atc"} for the average
#'   treatment effect on the untreated or control (ATU or ATC).
#' @param model The type of regression model, whether \code{ancova} for a
#'   an ANCOVA model, where the response in the post-period metric, or
#'   \code{did} for a difference-in-differences model, where the response is the
#'   difference between the post- and pre-period metrics.
#' @param include_pre Whether to include pre-period metrics in the regression.
#' @param interact_pre Whether to interact pre-period metrics with the treatment
#'   variable in the regression. Applicable only when \code{include_pre = TRUE}.
#' @param weights Name of variable containing inverse-probability-of-treatment
#'   weights (IPTW) for a weighted data sample. We recommend using a
#'   heteroskedasticity-corrected standard errors in this case, e.g.,
#'   \code{hc = "hc0"} instead of \code{hc = "none"}.
#' @param use_svyglm Whether to use \code{survey::svyglm} to estimate an IPTW
#'   treatment effect effect with its standard errors. The default is FALSE,
#'   in which case \code{stats::lm} is used.
#' @param paired Whether the data is paired, in which case a mixed-effects
#'   model with a random effect for the matching variable \code{match} is
#'   fitted using \code{lme4::lmer}.
#' @param var_equal Whether to assume unequal variances for the metrics across
#'   treatment and control groups. Similar to the \code{var.equal} argument in
#'   the \code{stats::t.test} function. Applies only when \code{paired = FALSE}.
#' @param match Name for variable for matching paired subjects.
#' @param cluster Name of variable for clustering the data. Defaults to NULL,
#'   in which case, a linear model will be fit using \code{stats::lm}. If
#'   \code{cluster} is not NULL, or \code{paired | !var_equal}, a linear
#'   mixed-effects model will be fit using \code{lme4::lmer}.
#' @param panel_time Name of time variable for longidudinal or panel data.
#' @param post_suffix Suffix for the post-period metric, which will have a
#'   string name equal to \code{paste0(y, post_suffix)}. Defaults to
#'   \code{"_post"}. Pass in the empty string \code{""} to omit the suffix.
#' @param conf_level Confidence level for inference.
#' @param output_models Whether to return the regression model objects together
#'   in addition to the effect table. In this case, a list will be returned
#'   containing the effect table in the first entry and a list of model objects
#'   in the second entry.
#' @param trim_fit Whether to apply the \code{TrimLinearModelFit} function to
#'   a model object. Currently works for \code{stats::lm} and
#'   \code{survey::svyglm} objects, i.e., only for the unpaired and
#'   equal-variance cases.
#' @param rel_eff Whether to compute the relative effect and its CI using the
#'   Fieller or delta method.
#' @param reff_ci_method A string for the relative effect CI method, either
#'   \code{"fieller"} for the Fieller CI, \code{"delta"} for the delta method
#'   without log transformation of the ratio, or \code{"delta_log"} for the
#'   delta method with log transformation of the ratio. See documentation for
#'   the \code{RatioCI} function for details.
#' @param hc Option for computing heteroskedasticity-corrected SEs using
#'   \code{car::hccm}. Defaults to \code{"none"}. The other values are
#'   \code{"hc0"}, \code{"hc1"}, \code{"hc2"}, \code{"hc3"}, and \code{"hc4"}.
#'   See documentation for \code{car::hccm} for details.
#' @param adjusted_means Whether to show the model-estimated average potential
#'   outcomes for, respectively, not receiving the treatment, \code{"mCa"}, and
#'   receiving the treatment, \code{"mTa"}.
#' @param ylab Label for column containing names of metrics.
#' @param percent_vars Names for variables to format as percentages. Note that
#'   this only changes the printing format, not the actual value of the
#'   variables.
#' @param format_vars Names for variables to format as a float. Note that this
#'   only changes the printing format, not the actual value of the variables.
#' @param digits Number of decimals for table statistics.
#' @return A data frame containing the table of regression-adjusted effects.
#' @references
#' \itemize{
#' \item Lin, Winston (2013). "Agnostic Notes on Regression Adjustments to
#'   Experimental Data: Reexamining Freedman’s Critique. \emph{The Annals of
#'   Applied Statistics}, 7:295–318.
#' }

AdjustedEffectTable <- function(data, y, vars = NULL, int_vars = NULL,
                                treat = "treat", unit = NULL, estimand = c(
                                  "ate", "att", "atu", "atc"
                                ), model = c("ancova", "did"),
                                include_pre = FALSE, interact_pre = FALSE,
                                weights = NULL, use_svyglm = FALSE,
                                paired = FALSE, var_equal = TRUE,
                                match = "match", cluster = NULL,
                                panel_time = NULL, post_suffix = "_post",
                                conf_level = 0.95, rel_eff = TRUE,
                                hc = c(
                                  "none", "hc0", "hc1", "hc2", "hc3", "hc4"
                                ), output_models = FALSE, trim_fit = TRUE,
                                reff_ci_method = c(
                                  "fieller", "delta", "delta_log"
                                ), adjusted_means = TRUE, ylab = "metric",
                                percent_vars = c(
                                  "mC", "mT", "diff", "mCa", "mTa", "eff", "se",
                                  "pval", "lb", "ub", "sd", "smd"
                                ), format_vars = "stat", digits = 1) {
  # Checks
  stopifnot(length(y) > 0)
  if (is.null(unit) && !is.null(panel_time)) {
    stop("Need to specify value for 'unit' if 'panel_time' is not null.")
  }
  stopifnot(is.numeric(data[[treat]]) || is.logical(data[[treat]]))
  data[[treat]] <- as.integer(data[[treat]])
  # Create unit identifier if not specified
  if (is.null(unit)) {
    unit <- "unit"
    data[[unit]] <- seq_len(dim(data)[1])
  }
  # Create post-period names of metrics
  y_post <- paste0(y, post_suffix)
  # Determine centering weights for covariates based on causal estimand
  estimand <- match.arg(estimand)
  center_weights <- switch(estimand,
    ate = NULL,
    att = data[[treat]],
    atu = 1 - data[[treat]],
    atc = 1 - data[[treat]]
  )
  # Create model data
  # (a) Outcome data
  # Only in the case of post_suffix == "", can we allow for factor variables
  # for y, in which case, it will be expanded to k - 1 variables, if it has
  # k levels. Otherwise, all the y variables will need to be numeric.
  if (nzchar(post_suffix)) {
    outcome_data <- data[y_post]
  } else {
    outcome_data <- CovariateData(data, vars = y_post, center = FALSE)
    y <- names(outcome_data)
    y_post <- y
  }
  # (b) Uncentered data
  uncentered_vars <- c(unit, treat, weights)
  if (paired) {
    uncentered_vars <- c(uncentered_vars, match)
  }
  if (!is.null(cluster)) {
    uncentered_vars <- c(uncentered_vars, cluster)
  }
  uncentered_data <- data[uncentered_vars]
  model_data <- outcome_data %>%
    bind_cols(uncentered_data)
  # (c) Pre-period data
  if (include_pre || model == "did") {
    y_pre <- paste0(y, "_pre")
    pre_data <- CovariateData(
      data = data, vars = y_pre, center = TRUE, center_weights = center_weights
    )
    model_data <- model_data %>%
      bind_cols(pre_data)
  }
  # (d) Additive covariate data (with panel time variable)
  if (!is.null(vars) || !is.null(panel_time)) {
    vars <- c(MultipleValuedVars(data, vars = vars), panel_time)
    vars_data <- CovariateData(
      data = data, vars = vars, center = TRUE, center_weights = center_weights
    )
    vars_expanded <- names(vars_data)
    model_data <- model_data %>%
      bind_cols(vars_data[setdiff(vars_expanded, names(model_data))])
  } else {
    vars_expanded <- NULL
  }
  # (e) Interacted covariate data
  if (is.null(int_vars)) {
    int_vars_expanded <- NULL
  } else {
    int_vars <- MultipleValuedVars(data, vars = int_vars)
    int_vars_data <- CovariateData(
      data = data, vars = int_vars, center = TRUE,
      center_weights = center_weights
    )
    int_vars_expanded <- names(int_vars_data)
    model_data <- model_data %>%
      bind_cols(int_vars_data[setdiff(int_vars_expanded, names(model_data))])
  }
  # Whether to fit an lmer model with random effect level for each observation,
  # needed for unpaired data with unequal group variances
  single_obs <- !paired & !var_equal

  # Generate model formulas
  model_formula <- AdjustedEffectFormula(
    y = y, treat = treat, vars = vars_expanded, int_vars = int_vars_expanded,
    model = model, include_pre = include_pre, interact_pre = interact_pre,
    paired = paired, var_equal = var_equal, match = match, cluster = cluster,
    panel_time = panel_time, unit = unit, post_suffix = post_suffix
  )
  # Fit models and extract model stats:
  # (a) If weights argument is not NULL, use survey::svyglm if use_svyglm is
  #     TRUE; otherwise, use stats::lm.
  # (b) If weights argument is NULL, use lmerTest::lmer when data is paired or
  #     group variances are unequal; otherwise, use stats::lm.
  # (c) When lmerTest::lmer is used, allow for single observation within a level
  #     when group variances are unequal for unpaired data.
  if (is.null(weights)) {
    design <- NULL
    lmer_condition <- (
      paired | !var_equal | !is.null(cluster) | !is.null(panel_time)
    )
    method <- ifelse(lmer_condition, "lmer", "lm")
  } else {
    design <- svydesign(
      id = ~1, weights = CreateFormula(x = weights), data = model_data
    )
    method <- ifelse(use_svyglm, "svyglm", "lm")
  }
  stats <- y %>%
    MyApply(function(.) {
      row <- ExtractModelStats(
        data = model_data, formula = model_formula[[.]], model = model,
        method = method, unit = unit, weights = weights, design = design,
        single_obs = single_obs, treat = treat, hc = hc,
        conf_level = conf_level, rel_eff = rel_eff,
        reff_ci_method = reff_ci_method, output_models = output_models,
        trim_fit = trim_fit
      )
      return(row)
    })
  names(stats) <- y
  # If output_models = TRUE, the stats object actually contains two objects,
  # one containing the stats and the other containing the model objects.
  # Extract the model objects into models, set the stats object to contain
  # only the stats part, and clear memory at the end.
  if (output_models) {
    models <- stats %>%
      MyApply(function(.) .[["models"]])
    stats <- stats %>%
      MyApply(function(.) .[["stats"]])
    gc()
  }
  # Convert stats object from list to a data frame and format it.
  stats <- stats %>%
    Unsplit(label = ylab) %>%
    data.frame()
  if (rel_eff) {
    percent_vars <- c(percent_vars, "reff", "re_lb", "re_ub")
  }
  stats <- FormatTable(
    data = stats, format_vars = format_vars, percent_vars = percent_vars,
    digits = digits
  )
  if (!adjusted_means) {
    stats[c("mCa", "mTa")] <- NULL
  }
  # Return both stats and models objects if output_models = TRUE. Otherwise,
  # return only the stats object.
  if (output_models) {
    return(list(stats = stats, models = models))
  } else {
    return(stats)
  }
}

#' Weighted Effect Table
#'
#' Wrapper for \code{AdjustedEffectTable} to estimate unadjusted effects from
#' inverse-probability-of-treatment-weighted (IPTW) samples, using either
#' \code{stats::lm} with Huber-White heteroskedasticity-corrected (HC) standard
#' errors or \code{survey::svyglm}, which also implements the HC errors. From
#' examples we have seen, we find that the \code{hc = "hc1"} option gives SE
#' estimates for \code{stats::lm} closest to those obtained using
#' \code{survey::svyglm}, although the other HC options tend to give very
#' similar results as well.
#'
#' @param data A data frame containing metrics for treatment and control
#'   subjects, an indicator variable for treatment assignment and optionally
#'   additional covariates and an integer variable for matching of pairs.
#' @param y Character vector containing names of variables to estimate effects
#'   for.
#' @param weights Name of variable containing inverse-probability-of-treatment
#'   weights (IPTW) for a weighted data sample. We recommend using a
#'   heteroskedasticity-corrected standard errors in this case, e.g.,
#'   \code{hc = "hc0"} instead of \code{hc = "none"}.
#' @param use_svyglm Whether to use \code{survey::svyglm} to estimate an IPTW
#'   treatment effect effect with its standard errors. The default is FALSE,
#'   in which case \code{stats::lm} is used.
#' @param hc Option for computing heteroskedasticity-corrected SEs using
#'   \code{car::hccm}. Defaults to \code{"hc1"}. The other values are
#'   \code{"hc0"}, \code{"hc2"}, \code{"hc3"}, \code{"hc4"}, and \code{"none"}.
#'   See documentation for \code{car::hccm} for details.
#' @param ... Other arguments passed into \code{AdjustedEffectTable}, not
#'   included above and not one of the following: \code{include_pre},
#'   \code{interact_pre}, and \code{post_suffix}, which are already have set
#'   values passed into the function. See documentation for
#'   \code{AdjustedEffectTable} for description of other arguments.
#' @return A data frame containing the table of unadjusted effects for the
#'   IPTW sample with relevant statistics.

WeightedEffectTable <- function(data, y, weights = "weights",
                                use_svyglm = FALSE, hc = c(
                                  "hc1", "hc0", "hc2", "hc3", "hc4", "none"
                                ), ...) {
  stopifnot(length(weights) > 0)
  hc <- match.arg(hc)
  out <- AdjustedEffectTable(
    data = data, y = y, weights = weights, use_svyglm = use_svyglm, hc = hc, ...
  )
  return(out)
}

#' Compare Model Effects
#'
#' Compare the treatment effects for different regression adjustment models,
#' including the following:
#' (1) Soriano's (2019) pre-post model
#' (2) The analysis of covariance (ANCOVA) model with equal slopes and equal
#'     error variances
#' (3) Lin's (2013) model, i.e., ANCOVA model with unequal slopes and equal
#'     error variances
#' (4) Unequal-variance model, i.e., ANOVA model unequal slopes and unequal
#'     error variances using a random effects model.
#' For models (2) through (4), we compute both Fieller and delta method
#' confidence intervals (CIs) for the relative lifts.
#'
#' @param data A data frame containing an indicator variable for treatment
#'   assignment passed into the \code{treat} argument and post-period and
#'   pre-period metrics for treatment and control subjects specified in
#'   \code{y}.
#' @param y Character vector containing names of variables to estimate effects
#'   for. This is a name of the variable without the post-period and pre-period
#'   suffixes passed into, respectively, \code{post_suffix} and {pre_suffix}.
#'   E.g., if the post- and pre-period metrics are, respectively,
#'   \code{"gce_adoption_post"} and \code{"gce_adoption_pre"}, then the name
#'   passed into the \code{y} argument should be \code{"gce_adoption"}.
#' @param treat Name of treatment assignment indicator variable.
#' @param post_suffix Suffix for the post-period metric, which will have a
#'   string name equal to \code{paste0(y, post_suffix)}. Defaults to
#'   \code{"_post"}.
#' @param post_suffix Suffix for the pre-period metric, which will have a
#'   string name equal to \code{paste0(y, pre_suffix)}. Defaults to
#'   \code{"_pre"}.
#' @param verbose Whether to print messages.
#' @return A data frame containing the table of statistics and estimated effects
#'   for the models with Fieller and delta method CIs for the relative lifts of
#'   the frequentist models.
#' #' @references
#' \itemize{
#' \item Soriano, Jacopo (2019). “Percent Change Estimation in Large Scale
#'   Online Experiments.” https://arxiv.org/pdf/1711.00562.pdf.
#' \item Lin, Winston (2013). "Agnostic Notes on Regression Adjustments to
#'   Experimental Data: Reexamining Freedman’s Critique. \emph{The Annals of
#'   Applied Statistics}, 7:295–318.
#' }

CompareModelEffects <- function(data, y, treat, post_suffix, pre_suffix, digits,
                                verbose = TRUE) {
  y_with_suffix <- paste0(y, post_suffix)
  models <- c(
    "pre_post", "ancova_fieller", "ancova_delta", "lin_fieller", "lin_delta",
    "uneq_var_fieller", "uneq_var_delta"
  )
  vars <- c(
    "metric", "nC", "nT", "mC", "mT", "diff", "eff", "lb", "ub", "pval", "reff",
    "re_lb", "re_ub"
  )
  out <- list()
  for (m in models) {
    if (verbose) {
      cat("Fitting", m, "model...\n")
    }
    if (m == "pre_post") {
      st <- system.time({
        out[[m]] <- PrePostEffectTable(
          data = data, y = y_with_suffix, treat = treat,
          post_suffix = post_suffix, pre_suffix = pre_suffix, digits = digits
        )
      })
    } else {
      reff_ci_method <- ifelse(
        grepl("_fieller", m, fixed = TRUE), "fieller", "delta"
      )
      var_equal <- !(grepl("uneq_var", m, fixed = TRUE))
      interact_pre <- !(grepl("ancova", m, fixed = TRUE))
      st <- system.time({
        out[[m]] <- AdjustedEffectTable(
          data = data, y = y, treat = treat, include_pre = TRUE,
          post_suffix = post_suffix, digits = digits,
          interact_pre = interact_pre, reff_ci_method = reff_ci_method,
          var_equal = var_equal
        )
      })
    }
    if (verbose) {
      print(st)
    }
  }
  out <- out %>%
    map(~ dplyr::select_at(., vars)) %>% # nolint
    Unsplit(label = "model")
  return(out)
}
