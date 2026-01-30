#' Package imports and internal aliases
#'
#' Internal namespace declarations and lightweight helpers.
#'
#' @name imports
#' @keywords internal
#' @import data.table
#' @import rstan
#' @importFrom stats rnorm rbinom runif rgamma setNames sd
#' @importFrom data.table setnames
#' @importFrom utils globalVariables modifyList
#' @importFrom dplyr %>% mutate group_by summarise summarize filter count left_join arrange slice ungroup row_number n
#' @importFrom ggplot2 ggplot aes geom_col labs theme_classic facet_wrap
#' @importFrom ggplot2 scale_color_manual theme_bw geom_point geom_boxplot position_jitter
#' @importFrom tidyr pivot_longer
#' @importFrom tibble tibble
#' @importFrom utils head
#' @noRd
NULL


#' Null-coalescing operator
#'
#' Returns \code{a} unless it is \code{NULL}, otherwise returns \code{b}.
#'
#' @param a,b Objects to choose from.
#' @return \code{a} if non-\code{NULL}; otherwise \code{b}.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Rescaled Gaussian infectivity profile
#'
#' Computes a Gaussian-shaped infectivity profile rescaled to have maximum 1.
#'
#' @param t Numeric vector of time points.
#' @param peak_day Numeric; day of peak infectivity.
#' @param width Numeric; width parameter (standard deviation) of the Gaussian.
#'
#' @return Numeric vector of rescaled infectivity values in \[0, 1\].
#' @keywords internal
#' @noRd
g_rescaled <- function(t, peak_day, width) {
  if (length(t) == 0) return(numeric(0))
  vals <- exp(-0.5 * ((t - peak_day) / width)^2)
  mx <- max(vals, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) return(vals)
  vals / mx
}


#' Simulate viral load trajectory
#'
#' Computes log10 viral load over time using a double-exponential model.
#'
#' @param t Numeric vector of time points (days since infection).
#' @param v_p Numeric; log10 peak viral load.
#' @param t_p Numeric; time to peak (days).
#' @param lambda_g Numeric; growth rate parameter.
#' @param lambda_d Numeric; decay rate parameter.
#'
#' @return Numeric vector of log10 viral loads.
#' @keywords internal
#' @noRd
simulate_viral_load_trajectory <- function(t, v_p, t_p, lambda_g, lambda_d) {
  vt <- 2 * 10^v_p / (exp(-lambda_g * (t - t_p)) + exp(lambda_d * (t - t_p)))
  return(log10(vt))
}


#' Simulate Ct (cycle threshold) trajectory
#'
#' Computes Ct values over time using a piecewise linear model.
#'
#' @param t Numeric vector of time points (days since infection).
#' @param Cpeak Numeric; peak (minimum) Ct value.
#' @param r Numeric; rise rate (Ct increase per day before peak).
#' @param d Numeric; decay rate (Ct increase per day after peak).
#' @param t_peak Numeric; time to peak (days).
#'
#' @return Numeric vector of Ct values.
#' @keywords internal
#' @noRd
simulate_Ct_trajectory <- function(t, Cpeak, r, d, t_peak) {
  ct <- ifelse(t <= t_peak, Cpeak + r * (t_peak - t), Cpeak + d * (t - t_peak))
  return(ct)
}


#' Draw random viral load parameters for a role
#'
#' Retrieves viral load trajectory parameters for a given role from a parameter list.
#'
#' @param role Character; one of "adult", "child", "toddler", "elderly".
#' @param params Named list of role-specific VL parameters (default: \code{default_VL_params}).
#'
#' @return Named list with elements \code{v_p}, \code{t_p}, \code{lambda_g}, \code{lambda_d}.
#' @keywords internal
#' @noRd
draw_random_VL_params <- function(role, params = default_VL_params) {
  p <- params[[role]]
  if (is.null(p)) stop("Unknown role in VL_params: ", role)
  list(v_p = p$v_p, t_p = p$t_p, lambda_g = p$lambda_g, lambda_d = p$lambda_d)
}


#' Draw random Ct parameters for a role
#'
#' Retrieves Ct trajectory parameters for a given role from a parameter list.
#'
#' @param role Character; one of "adult", "child", "toddler", "elderly".
#' @param params Named list of role-specific Ct parameters (default: \code{default_Ct_params}).
#'
#' @return Named list with elements \code{Cpeak}, \code{r}, \code{d}, \code{t_peak}.
#' @keywords internal
#' @noRd
draw_random_Ct_params <- function(role, params = default_Ct_params) {
  p <- params[[role]]
  if (is.null(p)) stop("Unknown role in Ct_params: ", role)
  list(Cpeak = p$Cpeak, r = p$r, d = p$d, t_peak = p$t_peak)
}


#' Generate household roles based on country profile
#'
#' Samples household composition (adults, children, toddlers, elderly) based
#' on probability parameters in a country/region profile.
#'
#' @param country_profile Named list with \code{prob_single_parent} (numeric),
#'   \code{prob_siblings} (length-3 probability vector for 0/1/2 additional children),
#'   \code{prob_elderly} (length-3 probability vector for 0/1/2 elderly).
#'   Defaults are used for missing elements.
#'
#' @return Character vector of roles (e.g., \code{c("adult", "adult", "child", "toddler")}).
#' @keywords internal
#' @noRd
generate_household_roles <- function(country_profile = list()) {
  base_profile <- list(
    prob_single_parent = 0,
    prob_siblings = c(0.40, 0.50, 0.10),
    prob_elderly = c(0.5, 0.25, 0.25)
  )
  profile <- modifyList(base_profile, country_profile)

  n_adults <- sample(1:2, 1, prob = c(profile$prob_single_parent, 1 - profile$prob_single_parent))
  roles <- rep("adult", n_adults)
  roles <- c(roles, "child")

  n_siblings_plus <- sample(0:2, 1, prob = profile$prob_siblings)
  if (n_siblings_plus > 0) roles <- c(roles, rep("toddler", n_siblings_plus))

  n_grandparents <- sample(0:2, 1, prob = profile$prob_elderly)
  if (n_grandparents > 0) roles <- c(roles, rep("elderly", n_grandparents))

  return(sample(roles))
}


#' Default viral load trajectory parameters by role
#'
#' A named list with role-specific parameters for simulating viral load
#' trajectories. Each role contains: \code{v_p} (peak viral load log10),
#' \code{t_p} (time to peak), \code{lambda_g} (growth rate),
#' \code{lambda_d} (decay rate).
#'
#' @format Named list with elements \code{adult}, \code{child}, \code{toddler}, \code{elderly}.
#' @keywords internal
#' @noRd
default_VL_params <- list(
  adult   = list(v_p = 4.14, t_p = 5.09, lambda_g = 2.31, lambda_d = 2.71),
  child   = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  toddler = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  elderly = list(v_p = 2.95, t_p = 5.1,  lambda_g = 3.15, lambda_d = 0.87)
)


#' Default Ct (cycle threshold) trajectory parameters by role
#'
#' A named list with role-specific parameters for simulating Ct value
#' trajectories. Each role contains: \code{Cpeak} (peak Ct value),
#' \code{r} (rise rate), \code{d} (decay rate), \code{t_peak} (time to peak).
#'
#' @format Named list with elements \code{adult}, \code{child}, \code{toddler}, \code{elderly}.
#' @keywords internal
#' @noRd
default_Ct_params <- list(
  child   = list(Cpeak = 33.3, r = 2.11, d = 1.38, t_peak = 5.06),
  toddler = list(Cpeak = 34,   r = 1.26, d = 1.27, t_peak = 4.75),
  adult   = list(Cpeak = 33,   r = 1.49, d = 1.22, t_peak = 5.14),
  elderly = list(Cpeak = 33,   r = 1.49, d = 1.22, t_peak = 5.14)
)


#' Default household composition profile
#'
#' A named list specifying probabilities for household composition:
#' \code{prob_single_parent} (probability of single-parent household),
#' \code{prob_siblings} (probabilities for 0, 1, 2 additional siblings/toddlers),
#' \code{prob_elderly} (probabilities for 0, 1, 2 elderly members).
#'
#' @format Named list with probability vectors.
#' @keywords internal
#' @noRd
default_household_profile <- list(
  prob_single_parent = 0,
  prob_siblings = c(0.40, 0.50, 0.10),
  prob_elderly = c(0.9, 0.08, 0.02)
)


#' Normalize role labels to canonical set
#'
#' Maps common variants to canonical roles: adult, child, toddler, elderly.
#'
#' @param x Character vector of role labels.
#' @return Character vector of normalized labels.
#' @keywords internal
#' @noRd
.norm_role <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x == "infant"]  <- "toddler"
  x[x == "baby"]    <- "toddler"
  x[x == "sibling"] <- "child"
  x[x == "kid"]     <- "child"
  x[x == "parent"]  <- "adult"
  x[x == "mother"]  <- "adult"
  x[x == "father"]  <- "adult"
  x[x == "elder"]   <- "elderly"
  x[x == "grandparent"] <- "elderly"
  x
}
