% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{Formatize}
\alias{Formatize}
\title{Apply Printing Formats to Columns in a Data Frame}
\usage{
Formatize(data, x = NULL, ...)
}
\arguments{
\item{data}{A data frame.}

\item{x}{A character vector containing names of columns to format.}

\item{...}{Additional arguments to pass to \code{formattable::formattable}.}
}
\value{
The data frame passed into the function with the
 \code{formattable::formattable} printing formats applied to the selected
 columns.
}
\description{
Apply printing formats from \code{formattable::formattable} to a selection of
columns in a data frame.
}
