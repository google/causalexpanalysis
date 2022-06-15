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

#' Horizontal Bar or Point Plot
#'
#' A generic wrapper function to create a horizontal bar or point plot using
#' \code{ggplot2::ggplot} with options for ordering the bars or points by a
#' specified variable. The function may be extended in the future to include
#' the option for horizontal boxplots.
#'
#' @param data A data frame containing data for the plot.
#' @param x Name of x-axis variable, which should be a numeric.
#' @param y Name of the y-axis variable, which should be categorical, i.e., a
#'   character vector or a factor.
#' @param type The \code{ggplot2} geom to use, whether \code{"bar"} for a bar
#'   plot, or \code{"point"} for a point plot.
#' @param order_by Name of numeric variable to order the y-axis variable by.
#' @param desc Whether to sort by descending order. Note that this specifies a
#'   descending order from top of plot to bottom, which is in contrast to the
#'   order that ggplot specifies, which is from bottom to top of plot.
#' @param facet Name of variable to facet by.
#' @param facet_ncol Number of columns for facetted plot.
#' @param text Whether to show the text label containing value of x-axis
#'   variable.
#' @param text_nudge Number of points to nudge text label.
#' @param errorbar_vars A character vector containing names of the min and max
#'   values of the errorbar. Vector must have at least two values if not NULL.
#' @param errorbar_height Height of errorbar.
#' @param errorbar_size Size of errorbar.
#' @param multi_fill Whether to fill bars with multiple colors or a single
#'   color.
#' @param vline Optional numeric vector to specify the x-intercepts of
#'   vertical lines.
#' @param vline Line type for vertical lines.
#' @param scale_log10 Whether to apply a log10 transformation to the x-axis.
#' @param scale_percent Whether format labels on the x-axis as percentages.
#' @param xlab Label for x-axis variable.
#' @param ylab Label for y-axis variable.
#' @return A \code{ggplot2::ggplot} horizontal bar plot.

HPlot <- function(data, x, y, type = c("bar", "point"), order_by = NULL,
                  desc = TRUE, facet = NULL, facet_ncol = NULL, text = FALSE,
                  text_nudge = 0, errorbar_vars = NULL,
                  errorbar_height = 0.2, errorbar_size = 0.2,
                  multi_fill = FALSE, vline = NULL, vline_type = 2,
                  scale_log10 = FALSE, scale_percent = FALSE, xlab = NULL,
                  ylab = NULL) {
  type <- match.arg(type)
  geom_fun <- switch(type,
    bar = geom_bar,
    point = geom_point
  )
  # Create modified data frame for plotting
  plot_data <- data.frame(x = data[[x]], y = data[[y]])
  if (!is.null(facet)) {
    plot_data$facet <- data[[facet]]
  }
  if (!is.null(errorbar_vars)) {
    stopifnot(length(errorbar_vars) >= 2)
    plot_data$lb <- data[[errorbar_vars[1]]]
    plot_data$ub <- data[[errorbar_vars[2]]]
  }
  if (!is.null(order_by)) {
    # Note that desc argument specifies a descending order from top of plot
    # to bottom, but ggplot specifies an order from bottom of plot to top.
    # Hence, decreasing = !desc.
    plot_data <- plot_data[order(data[[order_by]], decreasing = !desc), ]
    plot_data$y <- factor(plot_data$y, levels = unique(plot_data$y))
  }
  # Create base plot
  p <- plot_data %>%
    ggplot(aes(x = x, y = y))
  if (multi_fill) {
    p <- p + geom_fun(aes(fill = y), stat = "identity", show.legend = FALSE)
  } else {
    p <- p + geom_fun(aes(fill = "a"), stat = "identity", show.legend = FALSE)
  }
  # Add plot features if specified
  if (!is.null(facet)) {
    p <- p + facet_wrap(~facet, ncol = facet_ncol)
  }
  if (text) {
    p <- p + geom_text(aes(label = x), nudge_x = text_nudge)
  }
  if (!is.null(errorbar_vars)) {
    p <- p + geom_errorbarh(
      aes(xmin = lb, xmax = ub),
      height = errorbar_height,
      size = errorbar_size
    )
  }
  if (!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = vline_type)
  }
  if (!is.null(xlab)) {
    xlab <- paste0("\n", xlab)
    p <- p + xlab(xlab)
  }
  if (!is.null(ylab)) {
    ylab <- paste0(ylab, "\n")
    p <- p + ylab(ylab)
  }
  if (scale_log10) {
    scale_x_fun <- scale_x_log10
  } else {
    scale_x_fun <- scale_x_continuous
  }
  if (scale_percent) {
    p <- p + scale_x_fun(labels = scales::percent)
  } else if (scale_log10) {
    p <- p + scale_x_log10()
  }
  return(p)
}
HBarPlot <- function(data, x, y, ...) {
  out <- HPlot(data, x = x, y = y, type = "bar", ...)
  return(out)
}
HPointPlot <- function(data, x, y, ...) {
  out <- HPlot(data, x = x, y = y, type = "point", ...)
  return(out)
}

#' Covariate Balance Table Plot
#'
#' Creates a horizontal bar plot of the standardized mean difference (SMD) for
#' for the covariates in a covariate balance table.
#'
#' @param data A data frame containing data for covariate balance table.
#' @param type Type of plot.
#' @param order Order for the covariates.
#' @param smd_threshold Threshold for the SMD.
#' @param conf_level Confidence level for inference.
#' @param vline_type Line type for vertical lines.
#' @param scale_log10 Whether to apply a log10 transformation to the x-axis.
#' @param ylab Label for y-axis variable.
#' @return A \code{ggplot2::ggplot} horizontal bar plot of the SMD of the
#'   covariates, or the statistic or p-value for the t-test of the mean
#'   differences of the covariates.

BalanceTablePlot <- function(data, type = c("smd", "stat", "pval"),
                             order = c("none", "asc", "desc"),
                             smd_threshold = 0.1, conf_level = 0.95,
                             vline_type = 2, scale_log10 = FALSE,
                             ylab = "Covariates") {
  type <- match.arg(type)
  alpha <- 1 - conf_level
  switch(type,
    smd = {
      x <- "smd"
      xlab <- "Standardized Mean Difference (SMD)"
      vline <- smd_threshold * c(-1, 1)
    },
    stat = {
      x <- "stat"
      xlab <- "Test Statistic"
      vline <- qnorm(1 - alpha / 2) * c(-1, 1)
    },
    pval = {
      x <- "pval"
      xlab <- "Test p-Value"
      vline <- alpha
      scale_log10 <- TRUE
    }
  )
  order <- match.arg(order)
  p <- HBarPlot(data,
    x = x, y = "covariates",
    order_by = switch(order,
      none = NULL,
      x
    ), desc = (order == "desc"),
    vline = vline, vline_type = vline_type, scale_log10 = scale_log10,
    xlab = xlab, ylab = ylab
  )
  return(p)
}

#' Mean Estimate and Treatment Effect Table Plots
#'
#' Creates a horizontal bar or point plot of one of the following measures:
#' (1) A mean estimate of a response, specified using \code{type = "est"}.
#' (2) The absolute effect, specified using \code{type = "eff"}.
#' (3) The relative ffect, specified using  \code{type = "reff"}.
#' (4) The z-statistic of the absolute effect, specified using
#'     \code{type = "stat"}.
#' (5) The p-value of the absolute effect, specified using \code{type = "pval"}.
#' (6) The standardized mean difference (SMD) as a measure of effect size,
#'     specified using \code{type = "smd"}.
#'
#' Note that this function was originally named \code{EffectTablePlot} but has
#' been renamed to \code{EstimateTablePlot}, because it has been generalized
#' to include non-effect measures, e.g., the mean estimate of a response.
#' However, for backward compatibility's sake, the \code{EffectTablePlot}
#' function name works for the same function as well.
#'
#' @param data A data frame containing data for covariate balance table.
#' @param type Type of plot, depending on the variable to plot.
#' @param type The \code{ggplot2} geom to use, whether \code{"bar"} for a bar
#'   plot, or \code{"point"} for a point plot.
#' @param order Order for the covariates.
#' @param conf_level Confidence level for inference.
#' @param ref_slice The value for the \code{slice} variable in the data for an
#'   \code{est} plot, for use in determining the x-axis value of the vertical
#'   line.
#' @param vline_type Line type for vertical lines.
#' @param text Whether to show values as text labels.
#' @param scale_log10 Whether to apply a log10 transformation to the x-axis.
#' @param y Name of y variable. Defaults to \code{"metric"}.
#' @param ylab Label for y-axis variable.
#' @param smd_threshold Threshold for minimum effect size expressed by the
#'   standardized mean difference (SMD). If the string \code{"mdes"} is passed
#'   in instead, the minimum detectable effect size (MDES) will be used as the
#'   SMD threshold.
#' @param sig_level Significance level for computing the MDES.
#' @param power Power for computing the MDES.
#' @param n_vars A character vector of length 2 containing the names of the
#'   variables in \code{data} containing, respectively, the sample sizes for
#'   the treatment and control groups.
#' @param ... Additional arguments for the \code{HBarPlot} function.
#' @return A \code{ggplot2::ggplot} horizontal bar plot of statistic or p-value
#'   for the t-test of the mean differences of the covariates.

EstimateTablePlot <- function(data, type = c(
                                "eff", "reff", "stat", "pval", "smd", "est"
                              ), plot_type = c("bar", "point"),
                              order = c("none", "asc", "desc"),
                              conf_level = 0.95, ref_slice = "all",
                              vline_type = 2, text = FALSE, scale_log10 = FALSE,
                              scale_percent = FALSE, y = "metric",
                              ylab = "Metrics", smd_ci = FALSE,
                              smd_threshold = 0.02, sig_level = 0.05,
                              power = 0.80, n_vars = c("nT", "nC"), ...) {
  type <- match.arg(type)
  alpha <- 1 - conf_level
  if (type %in% c("eff", "reff")) {
    if (type == "eff") {
      x <- "eff"
      xlab <- "Mean Difference"
      errorbar_vars <- c("lb", "ub")
    } else {
      x <- "reff"
      xlab <- "Relative Mean Difference"
      errorbar_vars <- c("re_lb", "re_ub")
    }
    vline <- 0
  } else if (type == "stat") {
    x <- "stat"
    xlab <- "Test Statistic"
    vline <- qnorm(1 - alpha / 2) * c(-1, 1)
    errorbar_vars <- NULL
  } else if (type == "pval") {
    x <- "pval"
    xlab <- "Test p-Value"
    vline <- alpha
    errorbar_vars <- NULL
    scale_log10 <- TRUE
  } else if (type == "smd") {
    if (smd_ci) {
      smd_ci <- EffectSizeCI(
        d = data$smd, n1 = data[[n_vars[1]]], n2 = data[[n_vars[2]]],
        conf_level = conf_level, ci_vars = c("smd_lb", "smd_ub")
      )
      data <- data %>%
        bind_cols(smd_ci)
      errorbar_vars <- c("smd_lb", "smd_ub")
    } else {
      errorbar_vars <- NULL
    }
    if (smd_threshold == "mdes") {
      smd_threshold <- MinDetectableEffectSize(
        n1 = data[[n_vars[1]]], n2 = data[[n_vars[2]]], sig_level = sig_level,
        power = power
      )
    }
    x <- "smd"
    xlab <- "Standardized Mean Difference (SMD)"
    vline <- smd_threshold * c(-1, 1)
  } else if (type == "est") {
    x <- "est"
    xlab <- "Mean Estimate"
    vline <- data[data$slice == ref_slice, x]
    errorbar_vars <- c("lb", "ub")
  }
  order <- match.arg(order)
  if (scale_log10) {
    xlab <- paste(xlab, "(log10 Scale)")
  }
  p <- HPlot(data,
    x = x, y = y, type = plot_type, errorbar_vars = errorbar_vars, text = text,
    order_by = switch(order,
      none = NULL,
      x
    ), desc = (order == "desc"),
    vline = vline, vline_type = vline_type, scale_log10 = scale_log10,
    scale_percent = scale_percent, xlab = xlab, ylab = ylab, ...
  )
  return(p)
}
EffectTablePlot <- function(...) {
  out <- EstimateTablePlot(...)
  return(out)
}

#' Overlap Plots
#'
#' Create an "overlap" plot to evaluate the overlap of the propensity scores for
#' different study groups. The propensity score distributions are compared using
#' histograms, density plots, or both. This is essentially a wrapper around the
#' \code{personalized::check.overlap} function, which creates a
#' \code{ggplot2::ggplot} object. But the propensity score is passed directly
#' into the \code{personalized::check.overlap} function instead of having it
#' estimated by the function. Note the following for our implementation:
#' (1) The \code{type}, \code{bins}, and \code{alpha} arguments and their
#'     default values are the same as those in the
#'     \code{personalized::check.overlap} function.
#' (2) Since \code{personalized::check.overlap} function does not handle missing
#'     values, we impute missing values in the propensity score with the median
#'     value.
#'
#' @param data A data frame containing the propensity score and the treatment
#'   assignment variable, specified in, respectively, the \code{score} and
#'   \code{treat} arguments.
#' @param score String name of the propensity score vector in \code{data}.
#' @param treat String name of the treatment assignment vector in \code{data}.
#' @param type Type of plot, either \code{"histogram"}, \code{"density"}, or
#'   \code{"both"}.
#' @param bins Number of bins for the histogram. Defaults to 50.
#' @param alpha The \code{alpha} aesthetic for the \code{ggplot} object, to
#'   determine opacity of the geom.
#' @param title Optional argument for the plot title.
#' @param xlab Optional argument for the x-axis label.
#' @param ylab Optional argument for the y-axis label.
#' @return A \code{ggplot2::ggplot} object showing the overlap plot.

OverlapPlot <- function(data, score = "score", treat = "treat",
                        type = c("histogram", "density", "both"),
                        bins = 50L, alpha = ifelse(type == "both", 0.35, 0.5),
                        title = NULL, xlab = NULL, ylab = NULL) {
  # Check arguments
  stopifnot(c(score, treat) %in% names(data))
  type <- match.arg(type)
  # Impute missing score values with median
  score_median <- median(data[[score]], na.rm = TRUE)
  score_imputed <- ifelse(is.na(data[[score]]), score_median, data[[score]])
  # Create plot
  p <- personalized::check.overlap(
    x = score_imputed, trt = data[[treat]],
    propensity.func = function(x, trt) x, type = type, bins = bins,
    alpha = alpha
  )
  # Add plot title and axis labels if specified
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (!is.null(xlab)) {
    p <- p + xlab(xlab)
  }
  if (!is.null(ylab)) {
    p <- p + ylab(ylab)
  }
  return(p)
}

#' Group Time Series Plots
#'
#' Creates a list of group-aggregated time series plots for specified groups,
#' e.g., experiment arms, with two possible plot types:
#' (1) Line plots of group-aggregated statistics over time, using the mean
#'     statistic (default) or any user-specified function and the
#'     \code{ggplot2::geom_line} function.
#' (2) Smoothed-mean plots of variables over time, using the
#'     \code{method = "gam"} argument in \code{ggplot2::geom_smooth}.
#'
#' @param data A data frame containing data for the time series plots.
#' @param y A character vector of names for the time series metrics.
#' @param x Name of the variable in the x-axis, typically a date variable.
#' @param group Name of the group variable.
#' @param type Whether the plot is a line plot or a smooth plot.
#' @param FUN The function for aggregating \code{y} variables by group for line
#'   plots only. Defaults to \code{mean}. Note that for smooth plots, the only
#'   aggregation function available is the mean.
#' @param facet Name of facet variable.
#' @param facet_all Facet value for the entire dataset, if this argument is
#'   not null. The aggregated data will be appended with a dataset aggregated
#'   over all observations in the data with the facet variable value equal to
#'   \code{facet_all}.
#' @param facet_nrow The value passed into the \code{nrow} argument of the
#'   \code{ggplot2::facet_wrap} function.
#' @param facet_scales The value passed into the \code{scales} argument of the
#'   \code{ggplot2::facet_wrap} function.
#' @return A list of \code{ggplot2::ggplot} objects containing the line or
#'   smooth plots.

GroupTimeSeriesPlots <- function(data, y, x, group = "group",
                                 type = c("line", "smooth"), FUN = mean,
                                 facet = NULL, facet_all = NULL, facet_nrow = 1,
                                 facet_scales = "fixed",
                                 scale_percent = FALSE) {
  type <- match.arg(type)
  # Internal functions
  # (a) Data aggregation for line plots
  .AggregateData <- . %>%
    group_by_at(c(x, group, facet)) %>%
    summarize_at(y, FUN)
  # (b) Function to create a single plot
  .PlotFun <- function(.) {
    p <- ggplot(plot_data, aes_string(y = ., x = x, color = group))
    if (type == "line") {
      p <- p + geom_line()
    } else {
      p <- p + geom_smooth(formula = y ~ s(x), method = "gam")
    }
    if (!is.null(facet)) {
      facet_formula <- CreateFormula(x = facet)
      p <- p +
        facet_wrap(facet_formula, nrow = facet_nrow, scales = facet_scales)
    }
    if (scale_percent) {
      p <- p + scale_y_continuous(labels = scales::percent)
    }
    return(p)
  }
  # Aggregate data for line plots but not for smooth plots
  if (type == "line") {
    plot_data <- .AggregateData(data)
  } else {
    plot_data <- data
  }
  # If facet_all argument is specified, append with data having facet variable
  # set to facet_all value (aggregated for line plots but not for smooth plots).
  if (!is.null(facet_all)) {
    data[[facet]] <- facet_all
    if (type == "line") {
      plot_data <- bind_rows(plot_data, .AggregateData(data))
    } else {
      plot_data <- bind_rows(plot_data, data)
    }
  }
  # Create plot for each y variable
  out <- MyApply(y, .PlotFun)
  return(out)
}
