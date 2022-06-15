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

#' Effect Size from Mean Difference
#'
#' Compute the effect size as a measure of practical significance by
#' standardizing a mean difference. Implements two conventions:
#' (1) Cohen's d, which divides the mean difference by the pooled standard
#'     deviation of the two samples, as shown in, e.g., Cohen (1988).
#' (2) Standardized Difference, which divides the mean difference by the
#'     square root of the average of two sample variances, as proposed by
#'     Rosenbaum and Rubin (1985).
#'
#' @param m1 Mean of first sample.
#' @param m2 Mean of second sample.
#' @param s1 Standard deviation of first sample.
#' @param s2 Standard deviation of second sample.
#' @param n1 Size of first sample.
#' @param n2 Size of second sample.
#' @param method One of two choices, whether \code{cohens} for Cohen's d or
#'   \code{sdiff} for the standardized difference of
#' @return Effect size estimate from mean difference.
#' \itemize{
#' \item Cohen, Jacob (1988). \emph{Statistical Power Analysis for the
#'    Behavioral Sciences,} 2nd ed. New Jersey: Lawrence Erlbaum Associates.
#' \item Rosenbaum, Paul R. and Donald B. Rubin DB (1985). "Constructing a
#'    Control Group Using Multivariate Matched Sampling Methods That
#'    Incorporate the Propensity Score." \emph{The American Statistician.}
#'    39: 33-38.
#' }

EffectSize <- function(m1, m2 = 0, s1, s2 = NULL, n1 = NULL, n2 = NULL,
                       method = c("cohens", "sdiff")) {
  method <- match.arg(method)
  if (is.null(s2)) {
    sd <- s1
  } else {
    sd <- EffectSizeSD(
      s1, s2, n1, n2,
      method = switch(method,
        cohens = "pooled",
        sdiff = "average"
      )
    )
  }
  out <- (m1 - m2) / sd
  return(out)
}

#' Standard Deviation for Computing Effect Size of Mean Difference
#'
#' Computes the standard deviation to divide a mean difference for obtaining
#' an effect size, using either the pooled standard deviation, which gives us
#' the Cohen's d effect size, or the average standard deviation, which gives us
#' the standardized difference effect size. See documentation for
#' \code{EffectSize} for details.
#'
#' @param s1 Standard deviation of first sample.
#' @param s2 Standard deviation of second sample.
#' @param n1 Size of first sample.
#' @param n2 Size of second sample.
#' @param method One of two choices, whether \code{pooled} for the pooled
#'   standard deviation or \code{average} for the square root of the average
#'   sample variance.
#' @return Standard deviation for effect size calculation.

EffectSizeSD <- function(s1, s2, n1 = NULL, n2 = NULL,
                         method = c("pooled", "average")) {
  method <- match.arg(method)
  sd <- switch(method,
    pooled = sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2)),
    average = sqrt((s1^2 + s2^2) / 2)
  )
  return(sd)
}

#' Minimum Detectable Effect Size
#'
#' Computes the minimum detectable effect size for the power analysis of a
#' two-sample z test.
#'
#' @param n1 Size of the first sample, which is required.
#' @param n2 Size for the second sample. If \code{NULL}, then it is
#'   to be equal to the size \code{n1} of the first sample.
#' @param sig_level Significance level of the test.
#' @param power Power of the test.
#' @return The minimum detectable effect size for the two-sample z test.

MinDetectableEffectSize <- function(n1, n2 = NULL, sig_level = 0.05,
                                    power = 0.80) {
  if (is.null(n2)) {
    n2 <- n1
  }
  out <- (qnorm(1 - sig_level / 2) + qnorm(power)) *
    sqrt(1 / n1 + 1 / n2)
  return(out)
}

#' Effect Size Confidence Interval
#'
#' Computes the confidence interval (CI) for the Cohen's d standardized mean
#' difference (SMD) effect size, i.e., the estimated effect divided by the
#' estimated pooled standard deviation, based on its asymptotic normal
#' distribution, as shown on p. 86 of Hedges and Olkin (1985). To be accurate,
#' the CI formula derived by Hedges and Olkin (1985) is for Hedges' d, which
#' multiplies Cohen's d by a bias correction factor. However, the bias
#' correction factor is close to 1 for large sample sizes, so that both Hedges'
#' d and Cohen's d and their CIs are asymptotically equivalent. Hence, the CI
#' can be used for Cohen's d for large sample sizes.

#' @param d A numeric vector containing values of the Hedges' d standardized
#'   mean difference (SMD) effect size, or Cohen's d SMD effect size for large
#'   sample sizes.
#' @param n1 Size of the first sample, which is required.
#' @param n2 Size for the second sample, which is optional. If \code{NULL}, then
#'   it is set equal to the size \code{n1} of the first sample.
#' @param conf_level Confidence level for interval.
#' @param ci_vars Names for the lower and upper bounds of the CI in the output
#'   data frame.
#' @return A data frame containing the asymptotic CIs for Hedges' d and Cohen's
#'   d SMD effect sizes.
#' @references
#' \itemize{
#' \item Hedges, Larry and Ingram Olkin (1985). \emph{The Statistical Methods
#'   for Meta-analysis}. Orlando: Academic Press.
#' }

EffectSizeCI <- function(d, n1, n2 = NULL, conf_level = 0.95,
                         ci_vars = c("lb", "ub")) {
  if (is.null(n2)) {
    n2 <- n1
  }
  z_half_alpha <- qnorm((1 + conf_level) / 2)
  se <- sqrt((1 / n1 + 1 / n2) + d^2 / (2 * (n1 + n2)))
  out <- data.frame(d - z_half_alpha * se, d + z_half_alpha * se)
  names(out) <- ci_vars
  return(out)
}

#' Welch-Satterthwaitte Effective Degrees of Freedom
#'
#' Computes the Welch-Satterthwaitte effective degrees of freedom for a linear
#' combination of two independent sample variances, i.e., k1 * s1^2 + k2 * s2^2,
#' where s1^2, s2^2 are the sample variances and k1, k2 the coefficients.
#'
#' @param s1 Standard deviation of first sample.
#' @param s2 Standard deviation of second sample.
#' @param k1 Coefficient for the first sample in the linear combination.
#' @param k2 Coefficient for the first sample in the linear combination.
#' @param df1 Degrees of freedom for the first sample variance.
#' @param df2 Degrees of freedom for the second sample variance.
#' @return Effective degrees of freedom for the linear combination of sample
#'   variances.

WelchSatter <- function(s1, s2, k1, k2, df1, df2) {
  df <- (k1 * s1^2 + k2 * s2^2)^2 / (k1^2 * s1^4 / df1 + k2^2 * s2^4 / df2)
  return(df)
}

#' Standard Error of a Sample Proportion
#'
#' Computes the standard error of the sample binomial proportion given an
#  estimate of the population proportion and the sample size.
#'
#' @param p Estimate of the population proportion.
#' @param n The sample size.
#' @return Standard error estimate for the sample binomial proportion.

PropSE <- function(p, n) {
  out <- sqrt(p * (1 - p) / n)
  return(out)
}

#' Standard Error of the Difference between Two Independent-Sample Proportions
#'
#' Computes the standard error of the difference between the sample proportions
#' from two independent samples.
#'
#' @param p1 Estimate of the population proportion for the first sample.
#' @param p2 Estimate of the population proportion for the second sample.
#' @param n1 Size of the first sample.
#' @param n2 Size of the second sample.
#' @return Standard error estimate for the difference between the two
#'   independent-sample proportion.

PropDiffSE <- function(p1, p2, n1, n2) {
  out <- sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
  return(out)
}

#' Standard Error of the Difference between Two Matched-Sample Proportions
#'
#' Computes the standard error of the difference between the sample proportions
#' from two matched or paired samples. See Agresti and Min (2004).
#'
#' @param b Count of observations with 1 for first sample and 0 for second.
#' @param c Count of observations with 0 for first sample and 1 for second.
#' @param n The sample size.
#' @return Standard error estimate for the difference between the two
#'   matched-sample proportions.

MatchedPropDiffSE <- function(b, c, n) {
  out <- sqrt(((b + c) + (b - c)^2 / n) / n^2)
  return(out)
}

#' Add a Binomial Proportion Confidence Interval to a Data Frame
#'
#' Adds the lower and upper bounds of a binomial proportion confidence interval
#' to a data frame.
#'
#' @param data A data frame containing the sample proportions and sample sizes
#'   for computing the proportion confidence interval.
#' @param p Name of column containing the sample proportions.
#' @param n Name of column containing the sample sizes.
#' @param lb Name of column to contain lower bound of confidence interval.
#' @param ub Name of column to contain upper bound of confidence interval.
#' @param conf_level Confidence level for interval.
#' @return The data frame passed into the function with the binomial proportion
#'   confidence interval in the \code{lb} and \code{ub} columns.

AddPropCI <- function(data, p = "p", n = "n", lb = "lb", ub = "ub",
                      conf_level = 0.95) {
  z_alpha <- qnorm((1 + conf_level) / 2)
  margin <- z_alpha * PropSE(data[[p]], data[[n]])
  data[[lb]] <- data[[p]] - margin
  data[[ub]] <- data[[p]] + margin
  return(data)
}

#' Confidence Interval for the Ratio of Means
#'
#' Computes the confidence interval (CI) for a ratio of two means, given two
#' sample mean estimates and their covariance matrix, with the following
#' options:
#' (1) The Fieller CI, proposed in Fieller (1954).
#' (2) The delta method CI without log transformation of the ratio.
#' (3) The delta method CI with log transformation of the ratio.
#' See Hirschberg and Lye (2010) for a comparison of the first two approaches.
#'   Note that for the delta method without log transformation of the ratio,
#' the numerator and the denominator could be negative, as in the case for the
#' Fieller CI as well. But for the delta method with log transformation of the
#' ratio, both numerator and denominator need to be positive. The delta method
#' with log transformation of the ratio would not work for cases where the
#' numerator is zero, so it is less robust than the delta method without
#' log transformation for the ratio.
#'   The function is implemented in a "tidy" fashion by returning the results as
#' a data frame. The arguments \code{a, b, v11, v22, v12, s} can be numeric
#' vectors with multiple entries.
#'
#' @param method A string for the method name, either \code{"fieller"} for the
#'   Fieller CI or \code{"delta"} for the delta method without log
#'   transformation of the ratio.
#' @param a Sample estimate for mean in the numerator of the ratio.
#' @param b Sample estimate for mean in the denominator of the ratio.
#' @param v11 Scaled variance of a.
#' @param v22 Scaled variance of b.
#' @param v12 Scaled covariance of a and b.
#' @param s The scaling factor for the covariance matrix for [a, b]. Defaults
#'   to 1.
#' @param conf_level Confidence level for interval.
#' @param df Degree of freedom for t statistic for confidence interval.
#' @param chg Whether to convert the CI bounds from a ratio to a relative
#'   change, i.e., ratio minus 1.
#' @return A data frame containing the lower and upper bounds of the CI for the
#'   ratio of means.
#' @references
#' \itemize{
#' \item Fieller, E. C. (1954). "Some Problems in Interval Estimation."
#'    \emph{Journal of the Royal Statistical Society Series B,} 16: 175-185.
#' \item Hirschberg, Joe and Jenny Lye (2010). "A Geometric Comparison of the
#'    Delta and Fieller Confidence Intervals." \emph{The American Statistician,}
#'    64(3): 234-241.
#' }

RatioCI <- function(method = c("fieller", "delta", "delta_log"), a, b, v11, v22,
                    v12 = 0, s = 1, conf_level = 0.95, df = Inf, chg = FALSE) {
  method <- match.arg(method)
  if (method == "delta_log") {
    stopifnot(a > 0, b > 0)
  }
  theta <- a / b
  t_alpha <- qt(1 - (1 - conf_level) / 2, df = df)
  if (method == "fieller") {
    g <- t_alpha * s^2 * v22 / b^2
    C <- sqrt(
      v11 - (2 * theta * v12) + (theta^2 * v22) - g * (v11 - v12^2 / v22)
    )
    A0 <- 1 / (1 - g)
    A1 <- theta - (g * v12 / v22)
    A2 <- (t_alpha * s / b) * C
    lb <- A0 * (A1 - A2)
    ub <- A0 * (A1 + A2)
  } else {
    se_log_ratio <- s * sqrt((v11 / a^2) + (v22 / b^2) - (2 * v12 / (a * b)))
    if (method == "delta") {
      est <- theta
      se <- abs(theta) * se_log_ratio
    } else if (method == "delta_log") {
      est <- log(theta)
      se <- se_log_ratio
    }
    lb <- est - t_alpha * se
    ub <- est + t_alpha * se
    if (method == "delta_log") {
      lb <- exp(lb)
      ub <- exp(ub)
    }
  }
  out <- data.frame(lb = lb, ub = ub)
  if (chg) {
    out <- out - 1
  }
  return(out)
}

#' Delta Method Confidence Interval for the Ratio of Means
#'
#' Wrapper of the \code{RatioCI} function to compute the delta method
#' confidence interval (CI) for a ratio of two means, given two sample mean
#' estimates and their covariance matrix. See documentation for the
#' the \code{RatioCI} function.
#'
#' @param ... Parameters passed into the \code{RatioCI} function. For details,
#'   see parameter description for \code{RatioCI} function.
#' @return A data frame containing the lower and upper bounds of the delta
#    method CI for the ratio of two means.

DeltaRatioCI <- function(...) {
  out <- RatioCI(method = "delta", ...)
  return(out)
}
DeltaLogRatioCI <- function(...) {
  out <- RatioCI(method = "delta_log", ...)
  return(out)
}

#' Fieller Confidence Interval for the Ratio of Means
#'
#' Wrapper of the \code{RatioCI} function to compute the Fieller confidence
#' interval (CI) for a ratio of two means, given two sample mean estimates and
#' their covariance matrix. Implemented in a "tidy" fashion by returning the
#' results as a data frame. The arguments \code{a, b, v11, v22, v12, s} can be
#' numeric vectors with multiple entries.
#'
#' @param ... Parameters passed into the \code{RatioCI} function.
#' @return A data frame containing the lower and upper bounds of the
#'   Fieller CI.

FiellerCI <- function(...) {
  out <- RatioCI(method = "fieller", ...)
  return(out)
}

#' Fieller Confidence Interval for the Ratio of Proportions
#'
#' Wrapper of the \code{FiellerCI} function to computes the Fieller confidence
#' interval for a ratio of two proportions, given two sample proportions, their
#' sample sizes, and their mutual correlation. The arguments
#' \code{p1, p2, n1, n2, r} can be numeric vectors with multiple entries.
#'
#' @param p1 Estimate for the numerator of the ratio.
#' @param p2 Estimate for the denominator of the ratio.
#' @param n1 Size for the first sample.
#' @param n2 Size for the second sample.
#' @param r Estimate of correlation between p1 and p2.
#' @param conf_level Confidence level for interval.
#' @param ... Additional arguments to pass into \code{FiellerCI}
#' @return A data frame containing the lower and upper bounds of the
#'   Fieller confidence interval.

FiellerCIProp <- function(p1, p2, n1, n2, r = 0, conf_level = 0.95, ...) {
  v11 <- p1 * (1 - p1) / n1
  v22 <- p2 * (1 - p2) / n2
  out <- FiellerCI(
    a = p1, b = p2, v11 = v11, v22 = v22, v12 = 0, s = 1,
    conf_level = conf_level, df = Inf, ...
  )
  return(out)
}

#' Fieller Confidence Interval for the Ratio of Matched Proportions
#'
#' Wrapper of the \code{FiellerCI} function to compute the Fieller confidence
#' interval for a ratio of two matched or paired proportions, given the counts
#' of the 2x2 contingency table for their paired values. The arguments
#' \code{a, b, c, d} can be numeric vectors with multiple entries. Covariance
#' between two matched proportions are computed from the relation:
#' Cov(p1, p2) = (Var(p1) + Var(p2) - Var(p1 - p2)) / 2
#'
#' @param a A numeric vector of pair counts where both variables are 0
#' @param b A numeric vector of pair counts where first variable is 0 and the
#'   second 1.
#' @param c A numeric vector of pair counts where the first variable is 1 and
#'   the second 0.
#' @param d A numeric vector of pair counts where both variables are 1.
#' @param conf_level Confidence level for interval.
#' @param ... Additional arguments to pass into \code{FiellerCI}
#' @return A data frame containing the lower and upper bounds of the
#'   Fieller confidence interval.

FiellerCIMatchedProp <- function(a, b, c, d, conf_level = 0.95, ...) {
  n <- a + b + c + d
  p1 <- (b + d) / n
  p2 <- (c + d) / n
  v11 <- p1 * (1 - p1) / n
  v22 <- p2 * (1 - p2) / n
  v12 <- (v11 + v22 - MatchedPropDiffSE(b, c, n)^2) / 2
  V <- matrix(c(v11, rep(v12, 2), v22), ncol = 2)
  out <- FiellerCI(
    a = p1, b = p2, v11 = v11, v22 = v22, v12 = v12, s = 1,
    conf_level = conf_level, df = Inf, ...
  )
  return(out)
}

#' Effective Sample Size
#'
#' Compute the effective sample size based on the relative efficiency of the
#' estimation of a weighted mean. See http://go/effective-sample-size for
#' details.
#'
#' @param w A numeric vector of weights for the weighted mean.
#' @param all_stats Whether to return all the relevant statistics. If
#'   \code{TRUE}, a data frame will be returned with the effective sample size,
#'   the nominal sample size, and the relative efficiency. Otherwise, the scalar
#'   value of the effective sample size will be returned.
#' @return Either a data frame containing the effective and nominal sample
#'   sizes with their ratio, or a scalar value of the effective sample size,
#'   depending on the \code{all_stats} argument.

EffectiveSampleSize <- function(w, all_stats = FALSE) {
  ess <- sum(w)^2 / sum(w^2)
  if (all_stats) {
    out <- data.frame(n_eff = ess, n = length(w)) %>%
      mutate(efficiency = n_eff / n)
  } else {
    out <- ess
  }
  return(out)
}
