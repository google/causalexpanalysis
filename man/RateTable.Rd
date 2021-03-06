% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{RateTable}
\alias{RateTable}
\title{Rate Table}
\usage{
RateTable(
  data,
  metrics,
  group = NULL,
  label = "metrics",
  ci = FALSE,
  conf_level = 0.95,
  type = c("prop", "incid"),
  digits = 1,
  format = "f"
)
}
\arguments{
\item{data}{A data frame.}

\item{metrics}{A character vector containing the names of variables to obtain
rates for.}

\item{group}{A character vector containing names of variables to group by.}

\item{label}{The label for the column containing the names of the metrics
in the rate table.}

\item{ci}{Whether to include the confidence interval for the proportions.}

\item{conf_level}{Confidence level for interval.}

\item{type}{Type of rate, whether \code{prop} for a binomial proportion or
\code{incid} for a Poisson incidence rate.}

\item{digits}{Number of digits to print for percent variables.}

\item{format}{The \code{format} argument passed into the
\code{formattable::percent} function.}
}
\value{
A data frame containing the rate table.
}
\description{
Create a table of rates with the rate denominator and numerator and,
optionally, confidence interval (CI) for the rates, either a CI for a
binomial proportion or a Poisson incidence rate.
}
