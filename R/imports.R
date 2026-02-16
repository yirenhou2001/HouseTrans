#' Package imports and internal aliases
#'
#' Internal namespace declarations and lightweight helpers.
#'
#' @name imports
#' @keywords internal
#' @import data.table
#' @import rstan
#' @import ggplot2
#' @importFrom stats rnorm rbinom runif rgamma dgamma pgamma qgamma approx median setNames sd quantile
#' @importFrom data.table setnames as.data.table rbindlist
#' @importFrom utils globalVariables modifyList packageName
#' @importFrom dplyr %>% mutate group_by summarise summarize filter count left_join arrange
#' @importFrom dplyr slice ungroup row_number n select distinct rename case_when lag
#' @importFrom dplyr bind_rows inner_join pull n_distinct everything all_of any_of
#' @importFrom tidyr pivot_longer
#' @importFrom tibble tibble
#' @importFrom rlang sym .data
#' @importFrom utils head
#' @importFrom scales percent pretty_breaks breaks_pretty
#' @noRd
NULL

#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @importFrom dplyr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
#' @noRd
NULL
