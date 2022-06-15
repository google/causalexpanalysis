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

#' Unadjusted Sample
#'
#' Create an unadjusted sample for analysis, namely, by including only the
#' columns for treatment, outcomes, and covariates and removing rows where
#' the treatment variable is missing.
#'
#' @param data A data frame containing the treatment variable.
#' @param treat The name of the treatment variable.
#' @param outcomes Character vector containing names of outcomes.
#' @param covariates Character vector containing names of covariates.
#' @param ids Character vector containing names of ids.
#' @return A data frame containing the original data with rows having missing
#'   values for treatment removed.

UnadjustedSample <- function(data, treat = "treat", outcomes = NULL,
                             covariates = NULL, ids = NULL) {
  vars <- unique(c(treat, outcomes, covariates, ids))
  out <- data[!is.na(data[[treat]]), vars]
  return(out)
}

#' Weighted Sample
#'
#' Given the names of the treatment variable and covariates, return a data frame
#' with the propensity score and propensity score weights (PSWs) or balancing
#' weights attached to it. Relies on the \code{FitPropensityScoreModel} function
#' to compute the propensity score and the weights. Can also take a \code{psm}
#' object returned by the \code{FitPropensityScoreModel} function in the
#' \code{data} argument and extracts the data, score, and weights from the
#' the object. See the documentation for \code{FitPropensityScoreModel} for
#' details on the methods used to estimate the propensity score and weights.
#'
#' @param data A data frame containing the treatment variable and covariates
#'   for estimating the weights for the sample, or alternatively, a \code{psm}
#'   object returned by the \code{FitPropensityScoreModel} function.
#' @param treat The name of the treatment variable.
#' @param covariates Character vector containing names of covariates.
#' @param estimand The causal estimand for the propensity score weight or
#'   balancing weight to be applied to the sample.
#' @param score Name of the variable containing the estimated propensity score.
#' @param weights Name of weight variable.
#' @param type Deprecated name for the \code{estimand} argument. Kept here
#'   for backward compatibility.
#' @param ... Additional arguments to be passed into the
#'   \code{FitPropensityScoreModel} function.
#' @return A data frame containing the original model data with the propensity
#'   score and propensity score weights or balancing weights added to it.

WeightedSample <- function(data, treat = "treat", covariates = NULL,
                           estimand = c("att", "ate", "atu", "atc", "ato"),
                           score = "score", weights = "weights", type = NULL,
                           ...) {
  if (inherits(data, "psm")) {
    psm <- data
  } else {
    if (length(type) > 0 && type %in% c("att", "ate", "atu", "atc", "ato")) {
      estimand <- type
    } else {
      estimand <- match.arg(estimand)
    }
    psm <- FitPropensityScoreModel(
      data = data, treat = treat, covariates = covariates, estimand = estimand,
      score = score, weights = weights, ...
    )
  }
  out <- psm$data
  out[[score]] <- psm$data[[psm$score]]
  out[[weights]] <- psm$data[[psm$weights]]
  return(out)
}

#' Matched Sample
#'
#' Create a matched sample containing the treated observations and corresponding
#' matches from the control or untreated observations, with the following cases:
#' (1) Given the treatment variable and a measure to derive pairwise distances,
#'     e.g., a propensity score, and
#' (2) Given the treatment variable and a set of covariates.
#'
#' @param data A data frame containing the treatment variable and covariates
#'   for matching, or alternatively, a \code{psm} object returned by the
#'   \code{FitPropensityScoreModel} function.
#' @param treat The name of the treatment variable.
#' @param covariates Character vector containing names of covariates.
#' @param distance An optional numeric vector of distance values, such that a
#'   pairwise difference for any two observations gives the distance between
#'   them. This could be a vector of propensity scores or any other measures to
#'   derive pairwise distances.
#' @param estimand The causal estimand for the matching method to be applied
#'   to the sample. Currently, only the average treatment effect on the treated
#'   (ATT) population and the average treatment effect on the untreated or
#'   control (ATU or ATC) population are available.
#' @param match Name of the variable to indicate matching.
#' @param ... Additional arguments to be passed into the MatchIt::matchit
#'   function.
#' @return A data frame containing the treated sample and corresponding matches
#'   from the control sample.

MatchedSample <- function(data, treat = "treat", covariates = 1,
                          distance = NULL, estimand = c("att", "atu", "atc"),
                          match = "match", ...) {
  estimand <- match.arg(estimand)
  if (estimand == "atu") {
    estimand <- "atc"
  }
  if (inherits(data, "psm")) {
    distance <- as.numeric(data$data[[data$score]])
    treat <- data$treat
    data <- data$data
  }
  model_formula <- CreateFormula(treat, covariates)
  if (is.null(distance)) {
    fit <- MatchIt::matchit(
      model_formula,
      data = data, estimand = toupper(estimand), ...
    )
  } else {
    if (sum(is.na(distance)) > 0) {
      stop("The distance vector cannot have any missing values.")
    }
    fit <- MatchIt::matchit(
      model_formula,
      data = data, distance = distance,
      estimand = toupper(estimand), ...
    )
  }
  out <- MatchedData(fit, data = data, match_var = "match")
  return(out)
}

#' Matched Data Given MatchIt::matchit Object
#'
#' Create from a \code{MatchIt::matchit} object a data frame containing the
#' matched sample, with treatment subjects in the first rows followed by their
#' matched control subjects, and a \code{match_var} variable to indicate the
#' matching.
#'
#' @param matchit A \code{MatchIt::matchit} object.
#' @param data An optional data frame passed into the \code{MatchIt::matchit}
#'   function to create the \code{MatchIt::matchit} object. This is needed
#'   if the \code{MatchIt::matchit} object does not contain any modeling data.
#' @param match_var Name of the variable to indicate matching.
#' @return A data frame containing the matched sample.

MatchedData <- function(matchit, data = NULL, match_var = "match") {
  if (is.null(data)) {
    stopifnot(!is.null(matchit$model$data))
    data <- matchit$model$data
  }
  trt_data <- data[row.names(matchit$match.matrix), ]
  ctrl_data <- data[as.vector(matchit$match.matrix), ]
  out <- rbind(trt_data, ctrl_data)
  out[[match_var]] <- rep(seq(nrow(trt_data)), 2)
  return(out)
}

#' Weighted Data Sample [DEPRECATED: Don't Use.]
#'
#' Given a \code{psm} object, return a data frame with
#' inverse-probability-of-treatment weights (IPTW) added to it. The type of
#' IPTW can be specified for either the average treatment effect on the treated
#' (ATT) or the average treatment effect (ATE).
#'
#' @param psm An object of class \code{psm}, which could be created by the
#'   \code{FitPropensityScoreModel} function.
#' @param weights Name of weight variable.
#' @param type Type of the IPTW, whether ATT or ATE.
#' @return A data frame containing the original model data with the IPTW added
#'   to it.

WeightedData <- function(psm, weights = "weights", type = c("att", "ate")) {
  out <- psm$data
  out[[weights]] <- InvProbTreatWeight(
    treat = out[[psm$treat]], prob = predict(psm$model, type = "response"),
    type = type
  )
  return(out)
}

#' Fit Propensity Score Model
#'
#' Fit a propensity score model given treatment variable and covariates. Note
#' that we use the phrase "propensity score model" in a very broad sense
#' to include even models where the propensity score are not estimated directly
#' but can be derived from covariate-balancing weights under suitable
#' conditions. There are three broad categories for these models:
#' (1) Models that estimate the propensity score directly using a criterion
#'     that is agnostic to a specific causal estimand, e.g., maximizing a
#'     binomial likelihood or a penalized version of the likelihood.
#' (2) Models that estimate the propensity score directly using a criterion
#'     specific to a causal estimand, e.g., a covariate-balancing criterion
#'     optimizing for the average treatment effect on the entire population
#'     (ATE) or the average treatment effect on the treated population (ATT).
#' (3) Models that estimate the covariate-balancing weights directly without
#'     estimating a propensity score. These weights are optimized specifically
#'     for a causal estimand, e.g., the ATE or the ATT. They can be used to
#'     derive a propensity score if the causal estimand is the ATE. But if the
#'     the causal estimand is the ATT, we can derive the propensity score only
#'     for the untreated/control population. Similarly, if the causal estimand
#'     is the ATU/ATC, we can derive the propensity score only for the treated
#'     population.
#'   We use the \code{method} and \code{ps_model} arguments to categorize the
#' above three broad categories of models. The following \code{method} covers
#' models that estimate the propensity score directly using a criterion
#' not specific to a causal estimand:
#' (1) \code{ps}, covers models that estimate the propensity score directly
#'     using a criterion not specific to a causal estimand: The different models
#'     used under this method are identified using the \code{ps_model} argument.
#'   The following values for the \code{method} argument cover models that
#' estimate a propensity score using a criterion specific to a causal estimand:
#' (2) \code{gbm_es}, which fits a gradient-boosted tree model using
#'     \code{twang::ps}, with \code{version = "gbm"} and the effect size (ES)
#'     stopping method using \code{stop.method = "es.mean"}.
#' (3) \code{gbm_ks}, which fits a gradient-boosted tree model using
#'     \code{twang::ps}, with \code{version = "gbm"} and the Kolmogorov-Smirnov
#'     (KS) stopping method using \code{stop.method = "ks.max"}.
#' (4) \code{xgboost_es}, which fits an XGBoost model using \code{twang::ps},
#'     with \code{version = "xgboost"} and the ES stopping method using
#'     \code{stop.method = "es.mean"}.
#' (5) \code{xgboost_ks}, which fits an XGBoost model using \code{twang::ps},
#'     with \code{version = "xgboost"} and the KS stopping method using
#'     \code{stop.method = "ks.max"}.
#' (6) \code{cbps}, the covariate-balancing propensity score (CBPS) method of
#'     Imai et al. (2014), which maximizes covariate balance instead of the
#'     binomial likelihood of the propensity score model.
#' (7) \code{cbps_np}, the non-parametric version of the CBPS for continuous
#'     treatment but applied to binary treatments.
#'   The following values for the \code{method} argument cover models that
#' estimate the covariate-balancing weights directly without estimating a
#' propensity score:
#' (8) \code{cal}, the calibration weights of Chan et al. (2017).
#' (9) \code{ebal}, the entropy-balancing weights of Hainmuller (2012).
#'   For \code{method = "ps"}, the different propensity-score models are
#' identified using the \code{ps_model} argument, which has the following
#' values:
#' (1) \code{glm}, which fits a logistic regression model using
#'     \code{stats::glm} with \code{family = binomial}.
#' (2) \code{gam}, which fits a logistic regression model with penalized
#'     iterated reweighted least squares (PIRLS) algorithm using
#'     \code{mgcv::gam} with \code{family = binomial}.
#' (3) \code{glmnet}, which fits an elastic-net logistic regression model using
#'     \code{glmnet::cv.glmnet} with \code{family = "binomial"} and the
#'     elastic-net mixing parameter \code{alpha} specified as an argument.
#'   The propensity score weights are estimated for the following causal
#' estimands, passed into the \code{estimand} argument:
#' (1) \code{att}, the average treatment effect on the treated (ATT), available
#'     for all methods, except for \code{cbps_np}.
#' (2) \code{ate}, the average treatment effect on the entire population (ATE),
#'     available for all methods except \code{ebal}.
#' (3) \code{atu} or \code{atc}, the average treatment effect on the untreated
#'     or control (ATU or ATC), available for all methods, except for
#'     \code{cbps_np}.
#' (4) \code{ato}, the average treatment effect on the overlap population (ATO),
#'     available only for the \code{ps} method, as discussed in Li et al.
#'     (2018).
#'   The function returns the a \code{psm} object containing the following:
#' (1) \code{data}, the original data with the estimated propensity score and
#'     propensity score weights added to it.
#' (2) \code{model}, the model fit object.
#' (3) \code{method}, the method used to estimate the propensity score or
#'     propensity score weights.
#' (4) \code{ps_model}, the propensity score model used for the \code{ps}
#'     method.
#' (4) \code{estimand}, the causal estimand for estimating the score, if
#'     applicable, or for deriving the propensity score weights.
#' (5) \code{treat}, the name of the treatment variable.
#' (6) \code{score}, the name of the propensity score variable.
#' (7) \code{weights}, the name of the propensity score weight or
#'     balancing weight variable.
#'
#' @param data A data frame containing the treatment variable and covariates
#'   for fitting a model for propensity score or balancing weights.
#' @param treat The name of the treatment variable.
#' @param covariates Character vector containing names of covariates.
#' @param method The method for fitting the propensity score model, as listed
#'   earlier.
#' @param estimand The causal estimand for fitting the propensity score model,
#'   if applicable, or for fitting the balancing weight model. See values
#'   listed earlier in the above documentation.
#' @param score Name of the variable in the returned data that contains the
#'   estimated propensity score.
#' @param weights Name of the variable in the returned data that contains the
#'   estimated propensity score weights or balancing weights.
#' @param alpha The elastic-net mixing parameter between 0 and 1, with
#'   \code{alpha = 0} for ridge regression and \code{alpha = 1} for lasso,
#'   applicable for \code{method = "ps"} and \code{ps_model = "glmnet"}.
#' @param lambda The criterion to select the lambda of the elastic net, either
#'   \code{"min"} for the minimum lambda value or \code{"1se"} for lambda plus
#'   one standard error, applicable for \code{method = "ps"} and
#'   \code{ps_model = "glmnet"}.
#' @return An object of class \code{psm}, described earlier.
#' \itemize{
#' \item Chan, Kwun Chuen Gary, Sheung Chi Phillip Yam, and Zheng Zhang (2016).
#'   "Globally Efficient Non-parametric Inference of Average Treatment Effects
#'   by Empirical Balancing Calibration Weighting." \emph{Journal of the Royal
#'   Statistical Society: Series B (Statistical Methodology)}, 78(3): 673-700.
#' \item Hainmueller, Jens (2012). "Entropy Balancing for Causal Effects: A
#'   Multivariate Reweighting Method to Produce Balanced Samples in
#'   Observational Studies." \emph{Political Analysis}, 20(1): 25-46.
#' \item Imai, Kosuke and Marc Ratkovic (2014). "Covariate Balancing Propensity
#'   Score." \emph{Journal of the Royal Statistical Society B}, 76(1): 243-263.
#' \item Li, Fan, Kari Lock Morgan, and Alan M. Zaslavsky (2018). "Balancing
#'   Covariates via Propensity Score Weighting." \emph{Journal of the American
#'   Statistical Association}, 113(521): 390-400.
#' }

FitPropensityScoreModel <- function(data, treat = "treat", covariates = NULL,
                                    method = c(
                                      "ps", "gbm_es", "gbm_ks", "xgboost_es",
                                      "xgboost_ks", "cbps", "cbps_np", "cal",
                                      "ebal"
                                    ), ps_model = c("glm", "gam", "glmnet"),
                                    estimand = c(
                                      "att", "ate", "atu", "atc", "ato"
                                    ), score = "score", weights = "weights",
                                    alpha = 0, lambda = c("min", "1se"),
                                    normalized = FALSE) {
  method <- match.arg(method)
  estimand <- match.arg(estimand)
  # Check arguments
  if (method != "ps") {
    ps_model <- NA
    if (method == "cbps_np" && estimand != "ate") {
      warn_msg <- glue::glue(
        "Estimand chosen for '{method}' method must be the ATE. ",
        "Defaulting to ATE.\n"
      )
      cat(warn_msg)
      estimand <- "ate"
    } else if (method == "ebal" && estimand %in% c("ate", "ato")) {
      warn_msg <- glue::glue(
        "Estimand chosen for '{method}' method has to be ATT or ATU/ATC. ",
        "Defaulting to ATT.\n"
      )
      cat(warn_msg)
      estimand <- "att"
    } else if (estimand == "ato") {
      warn_msg <- glue::glue(
        "Estimand chosen for '{method}' method cannot be the ATO. ",
        "Defaulting to ATT.\n"
      )
      cat(warn_msg)
      estimand <- "att"
    }
  } else {
    ps_model <- match.arg(ps_model)
    if (ps_model == "glmnet") {
      lambda <- match.arg(lambda)
    }
  }
  # Specify constants
  n <- nrow(data)
  twang_methods <- c("gbm_es", "gbm_ks", "xgboost_es", "xgboost_ks")
  # Create model formula, data, and estimand variables. Convert all strings
  # into factors in model data, due to requirement for twang methods.
  model_formula <- CreateFormula(treat, covariates)
  model_data <- as.data.frame(unclass(data), stringsAsFactors = TRUE)
  model_estimand <- estimand
  # Switch treatment assignment for ATU/ATC estimands for some methods
  if (estimand %in% c("atu", "atc") &&
    method %in% c(twang_methods, "cal", "ebal")) {
    model_data[[treat]] <- 1 - model_data[[treat]]
    model_estimand <- "att"
  }
  # Get expanded design matrix for covariates for some methods
  if ((method == "ps" && ps_model == "glmnet") ||
    method %in% c("cal", "ebal")) {
    covar_mat <- CovariateData(
      data = model_data, vars = covariates, center = FALSE, as_matrix = TRUE
    )
  }
  # Fit model for propensity score or balancing weights
  if (method == "ps") {
    if (ps_model == "glm") {
      model <- glm(model_formula, data = model_data, family = binomial)
    } else if (ps_model == "gam") {
      model <- gam(model_formula, data = model_data, family = binomial)
    } else if (ps_model == "glmnet") {
      model <- cv.glmnet(
        x = covar_mat, y = model_data[[treat]], alpha = alpha,
        family = "binomial"
      )
    }
  } else if (method %in% twang_methods) {
    split_method <- unlist(strsplit(method, "_", fixed = TRUE))
    version <- split_method[1]
    stop_method <- switch(split_method[2],
      es = "es.mean",
      ks = "ks.max"
    )
    model <- twang::ps(
      formula = model_formula, data = model_data, n.trees = 5000,
      interaction.depth = 2, shrinkage = 0.01,
      estimand = toupper(model_estimand), version = version,
      stop.method = stop_method, n.minobsinnode = 10, n.keep = 1, n.grid = 25,
      ks.exact = NULL, verbose = FALSE
    )
  } else if (method %in% c("cbps", "cbps_np")) {
    if (method == "cbps") {
      # Estimate the ATE with value 0 when "ato" type is specified.
      cbps_estimand <- switch(estimand,
        ate = 0,
        att = 1,
        atu = 2,
        atc = 2,
        ato = 0
      )
      model <- CBPS::CBPS(model_formula, data = model_data, ATT = cbps_estimand)
    } else if (method == "cbps_np") {
      data[[treat]] <- factor(data[[treat]])
      model <- CBPS::npCBPS(model_formula, data = model_data)
    }
  } else if (method == "cal") {
    att_estimand <- (model_estimand == "att")
    # Note that we use dummy values for the Y argument, because we only need
    # the balancing weights, not the ATE.
    model <- ATE::ATE(
      Y = rep(1, nrow(model_data)), Ti = model_data[[treat]], X = covar_mat,
      ATT = att_estimand
    )
  } else if (method == "ebal") {
    model <- ebal::ebalance(Treatment = model_data[[treat]], X = covar_mat)
  }
  # Create psm object
  out <- list(
    data = data, model = model, method = method, ps_model = ps_model,
    treat = treat, score = score, weights = weights
  )
  class(out) <- "psm"
  # Extract propensity score estimates and corresponding weights
  # (a) Implementations that compute a score only
  if (method == "ps") {
    if (ps_model == "glmnet") {
      lambda_criterion <- model[[paste0("lambda.", lambda)]]
      est_score <- predict(
        model,
        type = "response", s = lambda_criterion, newx = covar_mat
      ) %>%
        as.numeric()
    } else {
      # We need to cast the predictions into a numeric vector, because
      # the predict.gam() function returns a 1-dimensional array instead of
      # a numeric vector.
      est_score <- as.numeric(predict(model, type = "response"))
    }
    est_weights <- Score2Weight(
      treat = data[[treat]], score = est_score, estimand = estimand
    )
    # (b) Implementations that compute both score and weight
  } else if (method %in% twang_methods) {
    est_score <- model$ps[[1]]
    # For ATU/ATC estimands, the propensity score estimated would have been
    # for the assignment to the untreated/control group. So take the unit
    # complement of the score to obtain the propensity score for the assignment
    # to the treated group.
    if (estimand %in% c("atu", "atc")) {
      est_score <- 1 - est_score
    }
    # We derive the weights ourselves, even though these weights are computed
    # by the twang::ps function.
    est_weights <- Score2Weight(
      treat = data[[treat]], score = est_score, estimand = estimand
    )
  } else if (method %in% c("cbps", "cbps_np")) {
    est_score <- model$fitted.values
    est_weights <- model$weights
    # (c) Implementations that compute a weight only
  } else if (method %in% c("cal", "ebal")) {
    # The weights and scores are computed using the treatment assignment in
    # \code{model_data[[treat]]}, which could be potentially different from
    # that in \code{data[[treat]]} if \code{estimand %in% c("atu", "atc")}.
    if (method == "cal") {
      if (att_estimand) {
        est_weights <- ifelse(model_data[[treat]], 1, model$weights.q)
      } else {
        est_weights <- model$weights.p + model$weights.q
      }
    } else if (method == "ebal") {
      est_weights <- rep(1, n)
      index <- seq_len(n)[model_data[[treat]] == 0]
      est_weights[index] <- model$w
    }
    est_score <- Weight2Score(
      treat = model_data[[treat]], weight = est_weights,
      estimand = model_estimand
    )
    # But if the original estimand is ATU/ATC, we obtain the unit complement
    # of the score instead.
    if (estimand %in% c("atu", "atc")) {
      est_score <- 1 - est_score
    }
  }
  # Normalize weights if specified
  if (normalized) {
    est_weights <- est_weights / sum(est_weights)
  }
  # Update returned data with score and weights
  out$data[[score]] <- est_score
  out$data[[weights]] <- est_weights
  # Return psm object
  return(out)
}

#' Propensity Score Weights
#'
#' Given a treatment indicator variable compute one of the following:
#' (1) Compute the propensity score weights (PSWs) given
#'     values for the probability of treatment, the propensity score.
#' (2) Compute the propensity score given the PSWs.
#' This is done for the following four kinds of causal estimands:
#' (1) The average treatment effect on the entire population (ATE),
#' (2) The average treatment effect on the treated (ATT),
#' (3) The average treatment effect on the untreated or control (ATU or ATC),
#'     and
#' (4) The average treatment effect on the overlap population (ATO).
#' Note that for the ATT, the probability of treatment can only be inferred for
#' the untreated or control population, while for the ATU/ATC, it can only be
#' inferred for the treated population.
#'
#' @param treat Numeric vector with values 1 for treatment and 0 otherwise.
#' @param score The propensity score, i.e., the probability of treatment.
#' @param weight The PSWs passed in as an argument to infer the probability of
#'   treatment.
#' @param estimand The causal estimand that determines the PSW, whether
#'   \code{att} for the average treatment effect on the treated (ATT),
#'   \code{ate} for the average treatment effect on the entire population (ATE),
#'   \code{atu} or \code{atc} for the average treatment effect on the untreated
#'   or control population (ATU or ATC), and \code{ato} for the average
#'   treatment effect on the overlap population (ATO).
#' @param normalized Whether to normalize weights. Using normalized and
#'   unnormalized weights correspond to using, respectively, the Hajek and
#'   Horvitz-Thompson estimators. Applicable only for converting probabilities
#'   to weights.
#' @return A numeric vector with PSWs or the propensity scores.
#' @references
#' \itemize{
#' \item Li, Fan, Kari Lock Morgan, and Alan M. Zaslavsky (2018). "Balancing
#'   Covariates via Propensity Score Weighting." \emph{Journal of the American
#'   Statistical Association}, 113(521): 390-400.
#' }

PropensityScoreWeight <- function(treat, score = NULL, weight = NULL,
                                  estimand = c(
                                    "att", "ate", "atu", "atc", "ato"
                                  ), normalized = FALSE) {
  stopifnot(!is.null(score) || !is.null(weight))
  estimand <- match.arg(estimand)
  # Convert probabilities to weights
  if (!is.null(score)) {
    # Tilting function
    if (estimand == "ate") {
      h <- 1
    } else if (estimand == "att") {
      h <- score
    } else if (estimand %in% c("atu", "atc")) {
      h <- 1 - score
    } else if (estimand == "ato") {
      h <- score * (1 - score)
    }
    # Normalization factors
    if (normalized) {
      factor1 <- 1 / sum(h * treat / score)
      factor0 <- 1 / sum(h * (1 - treat) / (1 - score))
    } else {
      factor1 <- 1
      factor0 <- 1
    }
    # Weights
    out <- h * (factor1 * treat / score + factor0 * (1 - treat) / (1 - score))
    # Convert weights to probabilities
  } else {
    if (estimand == "ate") {
      out <- ifelse(treat == 1, 1 / weight, 1 - 1 / weight)
    } else if (estimand == "att") {
      out <- ifelse(treat == 1, NA, weight / (1 + weight))
    } else if (estimand %in% c("atu", "atc")) {
      out <- ifelse(treat == 1, weight / (1 + weight), NA)
    } else if (estimand == "ato") {
      out <- ifelse(treat == 1, 1 - weight, weight)
    }
  }
  return(out)
}
Score2Weight <- function(treat, score,
                         estimand = c("att", "ate", "atu", "atc", "ato"),
                         normalized = FALSE) {
  out <- PropensityScoreWeight(
    treat = treat, score = score, estimand = estimand, normalized = normalized
  )
  return(out)
}
Weight2Score <- function(treat, weight,
                         estimand = c("att", "ate", "atu", "atc", "ato")) {
  out <- PropensityScoreWeight(
    treat = treat, weight = weight, estimand = estimand
  )
  return(out)
}
