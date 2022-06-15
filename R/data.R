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

#' Data Slices
#'
#' Given a data frame and a vector of names of slicing variables, return a
#' named list of data frames containing:
#' (1) The global slice of the original data frame, i.e., the original data
#'     frame itself.
#' (2) Different slices of the data frame, sliced by the unique values of the
#'     slicing variables.
#' If no slicing variables are provided, the function creates a list with one
#' element, i.e., the global slice of the data frame.
#'
#' @param data A data frame containing the variables specified in \code{vars}.
#' @param vars A character vector containing the string names of the slicing
#'     variables.
#' @param sep Separator character between variable name and its unique values,
#'   used as names for the entries of the returned vector. Defaults to the
#'   period.
#' @param global_slice_name Name for the global slice.
#' @return A named list of data frames, containing the global slice and other
#'   slices, sliced by the unique values of the slicing variables in
#'   \code{vars}.

DataSlices <- function(data, vars = NULL, sep = ".",
                       global_slice_name = "global") {
  out <- list()
  out[[global_slice_name]] <- data
  if (!is.null(vars)) {
    for (var in vars) {
      unique_values <- unique(data[[var]])
      for (value in unique_values) {
        slice_name <- paste(var, value, sep = sep)
        out[[slice_name]] <- data[data[[var]] == value, ]
      }
    }
  }
  return(out)
}

#' Unique Value Count
#'
#' Given a vector of column variables in a data frame, return a vector of the
#' number of unique values for each variable.
#'
#' @param data A data frame containing the variables specified in \code{vars}.
#' @param vars A character vector containing the string names of the variables.
#' @return A named numeric vector containing the number of unique values for
#'   variable, with the variable names in \code{vars} as the entry names.

UniqueValueCount <- function(data, vars) {
  out <- vars %>%
    sapply(function(x) {
      return(length(unique(data[[x]])))
    })
  return(out)
}

#' Multiple-valued Variables
#'
#' Given a vector of variable names in a data frame, return a vector with the
#' names of only the multiple-valued variables, i.e., variables with more than
#' one unique value.
#'
#' @param data A data frame containing the variables specified in \code{vars}.
#' @param vars A character vector containing the string names of the variables.
#' @return A character vector containing the names of only the multiple-valued
#' variables, i.e., variables with more than one unique value.

MultipleValuedVars <- function(data, vars) {
  unique_value_count <- UniqueValueCount(data, vars)
  out <- names(unique_value_count)[unique_value_count > 1]
  return(out)
}

#' Center a Variable
#'
#' Center a numeric vector \code{x} by subtracting from its entries a specified
#' weighted mean:
#' (1) If the weights \code{w} are not specified, the centering mean is just the
#'     mean of the vector, i.e., \code{mean(x)}.
#' (2) If the weights are specified, the centering mean is the following
#'     weighted mean: \code{sum(w * x) / sum(w)}.
#' One common application for centering by the above weighted mean is where
#' \code{w} is the treatment indicator in an analysis of covariance (ANCOVA)
#' model, so that:
#' (1) Centering all covariates by the treatment-weighted mean effectively
#'     centers them by their average in the treatment population, and
#' (2) Makes the intercept the counterfactual response for the treatment
#'     population.
#'
#' @param x A numeric vector containing values to center.
#' @param w A numeric vector containing weights to construct the centering
#'   mean to subtract from \code{x}. If \code{w} is not \code{NULL}, the
#'   function will stop if the sum of its values is zero, i.e.,
#'   \code{\sum(w) == 0}.
#' @return A centered version of \code{x}.

Center <- function(x, w = NULL) {
  if (is.null(w)) {
    center <- mean(x)
  } else {
    stopifnot(sum(w) != 0)
    center <- sum(w * x) / sum(w)
  }
  out <- x - center
  return(out)
}

#' Covariate Data
#'
#' Generate a data frame containing the values of the covariates, including
#' indicator variables for factor levels, using the \code{model.matrix} function
#' but with the intercept column removed.
#'
#' @param data A data frame containing predictors for regression model,
#'   including factors.
#' @param vars String names of the covariates in \code{data}. If its value is
#'   \code{NULL}, all the columns in \code{data} will be considered a covariate.
#' @param center Whether to center the covariates.
#' @param center_weights A numeric vector containing the weights to form the
#'   centering mean. See documentation for \code{Center} for details.
#' @param as_matrix Whether to return the covariate data as a matrix. If
#'   \code{FALSE}, a data frame is returned instead.
#' @return A data frame or a matrix containing the covariate values with
#'   factor contrasts expanded as indicator variables.

CovariateData <- function(data, vars = NULL, center = TRUE,
                          center_weights = NULL, as_matrix = FALSE) {
  if (is.null(vars)) {
    vars <- names(data)
  }
  model_formula <- CreateFormula(x = vars)
  model_data <- as.data.frame(model.matrix(model_formula, data))
  model_data[["(Intercept)"]] <- NULL
  if (center) {
    for (x in names(model_data)) {
      model_data[[x]] <- Center(model_data[[x]], w = center_weights)
    }
  }
  if (as_matrix) {
    model_data <- as.matrix(model_data)
  }
  return(model_data)
}

#' Bucketed Data
#'
#' Data frame of numeric data aggregated into means by buckets.
#'
#' @param data A data frame containing the treatment variable for the
#'   \code{treat} argument, the numeric variables to be aggregated for the
#'   \code{vars} argument, and optionally, the bucket variable for the
#'   \code{bucket_var} argument.
#' @param vars A character vector of the string names of the numeric variables
#'   to be aggregated.
#' @param treat String name of treatment variable.
#' @param bucket_var String name of bucket variable to bucket by, if it exists.
#'   Note that the bucket variable will always be used for bucketing if
#'   specified, regardless of the value specified for \code{buckets}.
#' @param buckets Number of buckets to randomly generate, if no value is
#'   specified for \code{bucket_var}. Number must lie in the range from 1 to
#'   \code{nrow(data)}.
#' @param seed Seed for random sampling.
#' @return A data frame containing numeric data aggregated into means by
#'   buckets, with treatment values and bucket ids as keys.

BucketedData <- function(data, vars, treat = "treat", bucket_var = NULL,
                         buckets = NULL, seed = NULL) {
  if (!is.null(bucket_var) || !is.null(buckets)) {
    n <- nrow(data)
    # Generate random buckets if no bucket variable is specified.
    if (is.null(bucket_var)) {
      stopifnot(buckets >= 1, buckets <= n)
      if (!is.null(seed)) {
        set.seed(seed)
      }
      data$bucket <- sample(buckets, size = n, replace = TRUE)
      bucket_var <- "bucket"
    }
    data <- data %>%
      group_by_at(c(treat, bucket_var)) %>%
      summarize_at(vars, mean)
  }
  return(data)
}

#' Frequency Table
#'
#' Create a table of frequency counts with, optionally, frequency proportions
#' and proportion confidence intervals.
#'
#' @param data A data frame.
#' @param group A character vector containing names of variables to group by to
#'   obtain frequencies.
#' @param unit Name of variable containing the key to the observation units.
#' @param unit_label The label for the observation units in the frequency table.
#' @param over A character vector containing names of variables to group by to
#'   obtain totals to divide frequencies by; could potentially be different
#'   what is specified in \code{group}.
#' @param prop Whether to include the frequency proportions.
#' @param ci Whether to include the confidence interval for the proportions.
#' @param conf_level Confidence level for interval.
#' @param digits Number of digits to print for percent variables.
#' @param format The \code{format} argument passed into the
#'   \code{formattable::percent} function.
#' @return A data frame containing the frequency table.

FreqTable <- function(data, group, unit = NULL, unit_label = NULL, over = NULL,
                      prop = TRUE, ci = FALSE, conf_level = 0.95, digits = 1,
                      format = "f") {
  if (is.null(unit)) {
    data$unit <- seq(nrow(data))
  } else {
    data$unit <- data[[unit]]
  }
  out <- data %>%
    group_by_at(group) %>%
    summarize(n = length(unique(unit))) %>%
    ungroup()
  if (prop) {
    percent_vars <- "prop"
    if (!is.null(over)) {
      out <- out %>%
        group_by_at(over)
    }
    out <- out %>%
      mutate(
        total = sum(n),
        prop = n / total
      )
    if (ci) {
      out <- AddPropCI(out, p = "prop", n = "total", lb = "lb", ub = "ub")
      percent_vars <- c(percent_vars, "lb", "ub")
    }
  }
  if (!is.null(unit_label)) {
    out <- out %>%
      rename_at("n", function(.) unit_label)
  }
  out <- data.frame(out)
  if (prop) {
    out <- Percentize(out, percent_vars, digits = digits, format = format)
  }
  return(out)
}

#' Rate Table
#'
#' Create a table of rates with the rate denominator and numerator and,
#' optionally, confidence interval (CI) for the rates, either a CI for a
#' binomial proportion or a Poisson incidence rate.
#'
#' @param data A data frame.
#' @param metrics A character vector containing the names of variables to obtain
#'   rates for.
#' @param group A character vector containing names of variables to group by.
#' @param label The label for the column containing the names of the metrics
#'   in the rate table.
#' @param ci Whether to include the confidence interval for the proportions.
#' @param conf_level Confidence level for interval.
#' @param type Type of rate, whether \code{prop} for a binomial proportion or
#'   \code{incid} for a Poisson incidence rate.
#' @param digits Number of digits to print for percent variables.
#' @param format The \code{format} argument passed into the
#'   \code{formattable::percent} function.
#' @return A data frame containing the rate table.

RateTable <- function(data, metrics, group = NULL, label = "metrics",
                      ci = FALSE, conf_level = 0.95, type = c("prop", "incid"),
                      digits = 1, format = "f") {
  type <- match.arg(type)
  stopifnot(length(metrics) > 0)
  out <- metrics %>%
    MyApply(function(x) {
      df <- data
      if (!is.null(group)) {
        df <- df %>%
          group_by_at(group)
      }
      df <- df %>%
        summarize_at(x, list(
          n = length,
          value = function(.) sum(., na.rm = TRUE)
        )) %>%
        mutate(rate = value / n)
      return(df)
    }) %>%
    Unsplit(label) %>%
    data.frame()
  percent_vars <- "rate"
  if (ci) {
    # TODO: Implement the CI for Poisson incidence rates.
    out <- AddPropCI(out, p = "rate", n = "n", lb = "lb", ub = "ub")
    percent_vars <- c(percent_vars, "lb", "ub")
  }
  out <- Percentize(out, percent_vars, digits = digits, format = format)
  return(out)
}
