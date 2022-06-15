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

#' Split Data Frame
#'
#' Wrapper function around \code{base::split} to split data frames into a list
#' of sub data frames using one or more grouping variables.
#'
#' @param data A data frame.
#' @param split A character vector containing the name(s) of the grouping
#'   variable(s) for splitting.
#' @param drop_split Whether the grouping variable(s) should be dropped from
#'   resulting list of data frames.
#' @param ... Additional arguments to pass to \code{base::split}.
#' @return A list of data frames split by the grouping variable(s).

Split <- function(data, split = NULL, drop_split = TRUE, ...) {
  stopifnot(length(split) > 0)
  split_list <- data[split]
  if (drop_split) {
    data[split] <- NULL
  }
  split_data <- split(data, split_list, ...)
  return(split_data)
}

#' Unsplit a List of Data Frames
#'
#' Combines a named list of similar data frames into a single data frame where
#' the first column contains the list names. Function is so named, because it is
#' the reverse operation of the \code{Split} function. It functions similarly to
#' the \code{base::unsplit} function but does not use it. Data frames are
#' converted to \code{data.table} objects for performance.
#'
#' @param data A data frame.
#' @param label The name of the first column added to data frame unsplit from
#'   the list.
#' @param drop_split Whether the grouping variable(s) should be dropped from
#'   resulting list of data frames.
#' @param ... Additional arguments to pass to \code{base::split}.
#' @return A list of data frames split by the grouping variable(s).

Unsplit <- function(data, label = "label") {
  # Add labels
  for (d in names(data)) {
    if (!(inherits(data[[d]], "data.table"))) {
      data.table::setDT(data[[d]])
    }
    existing_names <- names(data[[d]])
    data[[d]][[label]] <- d
    setcolorder(data[[d]], c(label, existing_names))
  }
  # Combine
  out <- data.frame(rbindlist(data))
  return(out)
}

#' Replace NAs in Data Frame
#'
#' Wrapper for the \code{tidyr::replace_na} function to make it easier for a
#' common use case, i.e., replace the NAs in most variables with a single value.
#' The highlight utility of this function is that we could apply the single
#' replacement value to all vectors of the same class as the replacement value,
#' e.g., all numeric vectors if the replacement value is \code{0}, all character
#' vectors if the replacement value is \code{"Missing"}, etc.
#'
#' @param data A data frame.
#' @param vars A character vector of names of variables to replace NAs for.
#' @param except A character vector of names of variables not to replace NAs
#'   for.
#' @param replace A single character string for the replacement value, or
#'   a named list of replacement strings. The names of the named list indicate
#'   the variables to apply the replacement values to. For a named list
#'   argument, \code{tidyr::replace_na} will be applied to the entire data frame
#'   with this argument without any changes.
#' @return A list of data frames split by the grouping variable(s).

ReplaceNA <- function(data, vars = NULL, except = NULL, replace = 0) {
  # If the replace argument is a list, apply the tidyr::replace_na function
  # to the entire data frame. Otherwise, apply it to selected vectors
  # in the data frame.
  if (inherits(replace, "list")) {
    out <- tidyr::replace_na(data, replace)
  } else {
    replace <- replace[1]
    # If the vars argument is empty, find the set of vectors of the same class
    # as the replace argument.
    if (is.null(vars)) {
      if (is.numeric(replace)) {
        FUN <- is.numeric
      } else if (is.character(replace)) {
        FUN <- is.character
      } else if (is.logical(replace)) {
        FUN <- is.logical
      }
      vars <- setdiff(names(data)[sapply(data, FUN)], except)
    }
    out <- data %>%
      mutate_at(vars, tidyr::replace_na, replace)
  }
  return(out)
}

#' Combine Strings
#'
#' Create all possible combinations of concatenations from two vectors of
#' strings, with the following variants:
#' (1) \code{Combine}, the main function for this operation.
#' (2) \code{Combine0}, the special case of \code{Combine} with \code{sep = ""}.
#' (3) \code{CombinePaths}, the special case of \code{Combine} with
#'   \code{sep = "/"}, for combining file paths.
#'
#' @param x The "left" vector of strings.
#' @param y The "right" vector of strings.
#' @param sep The separator string for concatenation.
#' @param lex Whether the combination obeys lexical order or not.
#' @return A vector of all possible concatenations of values from the two
#'   vectors.

Combine <- function(x, y, sep = "_", lex = TRUE) {
  if (lex) {
    out <- paste(
      rep(x, each = length(y)),
      rep(y, times = length(x)),
      sep = sep
    )
  } else {
    out <- paste(
      rep(x, times = length(y)),
      rep(y, each = length(x)),
      sep = sep
    )
  }
  return(out)
}
Combine0 <- function(x, y, lex = TRUE) {
  return(Combine(x, y, sep = "", lex = lex))
}
CombinePaths <- function(x, y, lex = TRUE) {
  return(Combine(x, y, sep = "/", lex = lex))
}

#' Create Formula Object from String Names
#'
#' Creates a formula object from string names of variables on left-hand side
#' (LHS) and right-hand side (RHS) of the formula.
#'
#' @param y A character vector containing names of the variables on the LHS of
#'   the formula.
#' @param x A character vector containing names of the variables on the RHS of
#'   the formula.
#' @param collapse A character for separating variables on either side of the
#'   the formula. Defaults to \code{+}.
#' @return A formula object with the LHS and RHS variables specified.

CreateFormula <- function(y = NULL, x = NULL, collapse = "+") {
  if (is.null(y) && is.null(x)) {
    stop("At least one of y and x needs to be non-null.")
  }
  if (is.null(y)) {
    y <- ""
  }
  y_formula <- paste(y, collapse = collapse)
  if (is.null(x)) {
    out <- as.formula(y_formula)
  } else {
    x_formula <- paste(x, collapse = collapse)
    out <- as.formula(paste(y_formula, x_formula, sep = "~"))
  }
  return(out)
}

#' Customized Apply Function
#'
#' A wrapper function around \code{base::sapply} but with the \code{simplify}
#' argument set to FALSE instead of TRUE. Useful for converted a character
#' vector into a list with the vector as names.
#'
#' @param x A vector.
#' @param fun The function to apply.
#' @param ... Additional arguments to pass to \code{base::sapply}.
#' @param use_names Value to pass into the \code{USE.NAMES} argument for
#'   \code{base::sapply}.
#' @return The result of applying the \code{base::sapply} function.

MyApply <- function(x, fun, ..., use_names = TRUE) {
  return(sapply(X = x, FUN = fun, ..., simplify = FALSE, USE.NAMES = use_names))
}

#' Binarize Numeric Vectors
#'
#' Convert a numeric vector into a binary variable taking values 1 and 0,
#' depending on whether the vector is greater than a threshold or otherwise.
#'
#' @param x A numeric vector.
#' @param threshold Numeric value for thresholding.
#' @param strict Whether inequality condition is strict, i.e., greater than or
#'   greater than or equal to.
#' @return A vector containing the binary indicator variable.

Binarize <- function(x, threshold = 0, strict = TRUE) {
  if (strict) {
    out <- (x > threshold)
  } else {
    out <- (x >= threshold)
  }
  return(1 * out)
}

#' Apply Printing Formats to Columns in a Data Frame
#'
#' Apply printing formats from \code{formattable::formattable} to a selection of
#' columns in a data frame.
#'
#' @param data A data frame.
#' @param x A character vector containing names of columns to format.
#' @param ... Additional arguments to pass to \code{formattable::formattable}.
#' @return The data frame passed into the function with the
#'  \code{formattable::formattable} printing formats applied to the selected
#'  columns.

Formatize <- function(data, x = NULL, ...) {
  if (!is.null(x)) {
    data[x] <- data[x] %>%
      lapply(formattable::formattable, ...)
  }
  return(data)
}

#' Apply Percent Printing Format to Columns in a Data Frame
#'
#' Apply percent printing format from \code{formattable::percent} to a selection
#' of columns in a data frame.
#'
#' @param data A data frame.
#' @param x A character vector containing names of columns to format.
#' @param ... Additional arguments to pass to \code{formattable::percent}.
#' @return The data frame passed into the function with the
#'  \code{formattable::percent} printing format applied to the selected columns.

Percentize <- function(data, x = NULL, ...) {
  if (!is.null(x)) {
    data[x] <- data[x] %>%
      lapply(formattable::percent, ...)
  }
  return(data)
}

#' Apply Numeric and Percent Printing Formats to Columns in a Data Frame
#'
#' Apply numeric printing format from \code{formattable::formattable} and
#' percent printing format from \code{formattable::percent} to a selection
#' of columns in a data frame.
#'
#' @param data A data frame.
#' @param format_vars A character vector containing names of columns to apply
#'   numeric printing format to.
#' @param percent_vars A character vector containing names of columns to apply
#'   percent printing format to.
#' @param digits The number of digits for the printing numbers.
#' @return The data frame passed into the function with the
#'   numeric and percent printing formats applied to the selected columns.

FormatTable <- function(data, format_vars = NULL, percent_vars = NULL,
                        digits = 1) {
  if (!is.null(format_vars)) {
    data <- Formatize(data, format_vars, digits = digits, format = "f")
  }
  if (!is.null(percent_vars)) {
    data <- Percentize(data, percent_vars, digits = digits, format = "f")
  }
  return(data)
}

#' Distinct Value Frequency Table
#'
#' Create a table of frequencies of units having k number of distinct values of
#' a specific variable \code{x}.
#'
#' @param data A data frame.
#' @param x Variable that we would like to obtain frequencies for its distinct
#'   values.
#' @param unit Variable specifying the IDs for the frequency units.
#' @return A data frame containing frequencies for each number of distinct
#'   values of the variable \code{x}.

DistinctValueFreqTable <- function(data, x, unit) {
  out <- data %>%
    group_by_at(unit) %>%
    summarize_at(x, n_distinct) %>%
    ungroup() %>%
    group_by_at(x) %>%
    summarize_at(unit, n_distinct) %>%
    data.frame()
  return(out)
}
