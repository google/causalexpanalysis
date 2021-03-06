% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{WelchSatter}
\alias{WelchSatter}
\title{Welch-Satterthwaitte Effective Degrees of Freedom}
\usage{
WelchSatter(s1, s2, k1, k2, df1, df2)
}
\arguments{
\item{s1}{Standard deviation of first sample.}

\item{s2}{Standard deviation of second sample.}

\item{k1}{Coefficient for the first sample in the linear combination.}

\item{k2}{Coefficient for the first sample in the linear combination.}

\item{df1}{Degrees of freedom for the first sample variance.}

\item{df2}{Degrees of freedom for the second sample variance.}
}
\value{
Effective degrees of freedom for the linear combination of sample
  variances.
}
\description{
Computes the Welch-Satterthwaitte effective degrees of freedom for a linear
combination of two independent sample variances, i.e., k1 * s1^2 + k2 * s2^2,
where s1^2, s2^2 are the sample variances and k1, k2 the coefficients.
}
