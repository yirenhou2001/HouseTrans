#' @title Helper Functions and Defaults
#' @description Utility functions, default parameters, and internal helpers
#' @name helpers
#' @keywords internal
NULL


# ==============================================================================
# NULL-COALESCING OPERATOR
# ==============================================================================

#' Null-coalescing operator
#'
#' Returns \code{a} unless it is \code{NULL}, otherwise returns \code{b}.
#'
#' @param a,b Objects to choose from.
#' @return \code{a} if non-\code{NULL}; otherwise \code{b}.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


# ==============================================================================
# VIRAL LOAD / CT TRAJECTORY FUNCTIONS
# ==============================================================================

#' Simulate viral load trajectory (log10)
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
  return(log10(pmax(1e-9, vt)))
}


#' Simulate Ct (cycle threshold) trajectory
#'
#' Computes Ct values over time using a piecewise linear model.
#' Note: Lower Ct = higher viral load.
#'
#' @param t Numeric vector of time points (days since infection).
#' @param Cpeak Numeric; peak (minimum) Ct value.
#' @param r Numeric; rise rate (Ct decrease per day before peak).
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


#' Draw viral load parameters for a role
#'
#' Retrieves viral load trajectory parameters for a given role from a parameter list.
#'
#' @param role Character; one of "adult", "infant", "toddler", "elderly".
#' @param params Named list of role-specific VL parameters.
#'
#' @return Named list with elements \code{v_p}, \code{t_p}, \code{lambda_g}, \code{lambda_d}.
#' @keywords internal
#' @noRd
draw_random_VL_params <- function(role, params = default_VL_params) {
  p <- params[[role]]
  if (is.null(p)) stop("Unknown role in VL_params: ", role)
  list(v_p = p$v_p, t_p = p$t_p, lambda_g = p$lambda_g, lambda_d = p$lambda_d)
}


#' Draw Ct parameters for a role
#'
#' Retrieves Ct trajectory parameters for a given role from a parameter list.
#'
#' @param role Character; one of "adult", "infant", "toddler", "elderly".
#' @param params Named list of role-specific Ct parameters.
#'
#' @return Named list with elements \code{Cpeak}, \code{r}, \code{d}, \code{t_peak}.
#' @keywords internal
#' @noRd
draw_random_Ct_params <- function(role, params = default_Ct_params) {
  p <- params[[role]]
  if (is.null(p)) stop("Unknown role in Ct_params: ", role)
  list(Cpeak = p$Cpeak, r = p$r, d = p$d, t_peak = p$t_peak)
}


#' Compute theoretical Ct trajectory
#'
#' @param t Days since infection.
#' @param p List with Cpeak, r, d, t_peak.
#' @return Ct value.
#' @keywords internal
#' @noRd
.get_theoretical_ct <- function(t, p) {
  ifelse(t <= p$t_peak,
         p$Cpeak + p$r * (p$t_peak - t),
         p$Cpeak + p$d * (t - p$t_peak))
}


#' Compute theoretical log10 viral load trajectory
#'
#' @param t Days since infection.
#' @param p List with v_p, t_p, lambda_g, lambda_d.
#' @return Log10 viral load.
#' @keywords internal
#' @noRd
.get_theoretical_vl <- function(t, p) {
  numerator <- 2 * 10^(p$v_p)
  denominator <- exp(-p$lambda_g * (t - p$t_p)) + exp(p$lambda_d * (t - p$t_p))
  vt <- numerator / denominator
  return(log10(pmax(1e-9, vt)))
}


#' Calculate curve value for imputation
#'
#' @param t Days since infection.
#' @param p Parameters list.
#' @param is_ct Logical; TRUE for Ct, FALSE for VL.
#' @return Curve value.
#' @keywords internal
#' @noRd
.calc_curve_val <- function(t, p, is_ct) {
  if (is_ct) {
    return(ifelse(t <= p$t_peak,
                  p$Cpeak + p$r * (p$t_peak - t),
                  p$Cpeak + p$d * (t - p$t_peak)))
  } else {
    val <- 2 * 10^p$v_p / (exp(-p$lambda_g * (t - p$t_p)) + exp(p$lambda_d * (t - p$t_p)))
    return(log10(pmax(1e-9, val)))
  }
}


# ==============================================================================
# HOUSEHOLD GENERATION
# ==============================================================================

#' Generate household roles based on country profile
#'
#' Samples household composition (adults, infants, toddlers, elderly) based
#' on probability parameters in a country/region profile.
#'
#' @param country_profile Named list with:
#'   \describe{
#'     \item{prob_adults}{Length-3 probability vector for 0/1/2 adults (default: c(0, 0, 1) = always 2 adults)}
#'     \item{prob_infant}{Probability of having an infant, 0-1 (default: 1.0 = always 1 infant)}
#'     \item{prob_siblings}{Length-3 probability vector for 0/1/2 toddlers/siblings (default: c(0.1, 0.7, 0.2))}
#'     \item{prob_elderly}{Length-3 probability vector for 0/1/2 elderly (default: c(0.7, 0.15, 0.15))}
#'   }
#'   Defaults are used for missing elements.
#'
#' @return Character vector of roles (e.g., \code{c("adult", "adult", "infant", "toddler")}).
#' @keywords internal
#' @noRd
generate_household_roles <- function(country_profile = list()) {
  # Default Profile
  # prob_adults: Default is 0% chance of 0, 0% chance of 1, 100% chance of 2
  # prob_infant: Default is 100% chance of 1
  base_profile <- list(
    prob_adults   = c(0.0, 0, 1),
    prob_infant   = 1.0,
    prob_siblings = c(0.1, 0.7, 0.2),
    prob_elderly  = c(0.7, 0.15, 0.15)
  )

  if (is.null(country_profile)) country_profile <- list()

  # Merge lists safely (prefer user input over base)
  profile <- utils::modifyList(base_profile, country_profile)

  roles <- c()

  # 1. Add Adults (0, 1, or 2)
  n_adults <- sample(0:2, 1, prob = profile$prob_adults)
  if (n_adults > 0) roles <- c(roles, rep("adult", n_adults))

  # 2. Add Infant (0 or 1) - Controlled by prob_infant
  has_infant <- stats::rbinom(1, 1, profile$prob_infant)
  if (has_infant == 1) roles <- c(roles, "infant")

  # 3. Add Toddlers (Siblings)
  n_toddlers <- sample(0:2, 1, prob = profile$prob_siblings)
  if (n_toddlers > 0) roles <- c(roles, rep("toddler", n_toddlers))

  # 4. Add Elderly
  n_elderly <- sample(0:2, 1, prob = profile$prob_elderly)
  if (n_elderly > 0) roles <- c(roles, rep("elderly", n_elderly))

  return(roles)
}


# ==============================================================================
# ROLE NORMALIZATION
# ==============================================================================

#' Normalize role labels to canonical set
#'
#' Maps common variants to canonical roles: adult, infant, toddler, elderly.
#'
#' @param x Character vector of role labels.
#' @return Character vector of normalized labels.
#' @keywords internal
#' @noRd
.norm_role <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x == "child"]   <- "infant"
  x[x == "baby"]    <- "infant"
  x[x == "sibling"] <- "toddler"
  x[x == "kid"]     <- "toddler"
  x[x == "parent"]  <- "adult"
  x[x == "mother"]  <- "adult"
  x[x == "father"]  <- "adult"
  x[x == "elder"]   <- "elderly"
  x[x == "grandparent"] <- "elderly"
  x
}


# ==============================================================================
# PRIOR PARSING
# ==============================================================================

#' Parse prior specification to Stan format
#'
#' Converts user-friendly prior specification to Stan-compatible format.
#'
#' @param p_list List with \code{dist} and \code{params} elements.
#' @param def_type Default prior type (1=Normal, 2=Uniform, 3=LogNormal).
#' @param def_params Default parameter values.
#'
#' @return List with \code{type} (integer) and \code{params} (numeric vector).
#' @keywords internal
#' @noRd
.parse_prior <- function(p_list, def_type, def_params) {
  if (is.null(p_list)) return(list(type = def_type, params = def_params))

  type_int <- 1L  # Default: Normal
  if (!is.null(p_list$dist)) {
    dist <- tolower(p_list$dist)
    if (dist == "normal") type_int <- 1L
    else if (dist == "uniform") type_int <- 2L
    else if (dist == "lognormal") type_int <- 3L
  }

  params <- if (!is.null(p_list$params)) p_list$params else def_params
  list(type = type_int, params = params)
}


# ==============================================================================
# DEFAULT PARAMETERS
# ==============================================================================

#' Default viral load trajectory parameters by role
#'
#' A named list with role-specific parameters for simulating viral load
#' trajectories. Each role contains: \code{v_p} (peak viral load log10),
#' \code{t_p} (time to peak), \code{lambda_g} (growth rate),
#' \code{lambda_d} (decay rate).
#'
#' @format Named list with elements \code{adult}, \code{infant}, \code{toddler}, \code{elderly}.
#' @keywords internal
#' @noRd
default_VL_params <- list(
  adult   = list(v_p = 4.14, t_p = 5.09, lambda_g = 2.31, lambda_d = 2.71),
  infant  = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  toddler = list(v_p = 5.84, t_p = 4.09, lambda_g = 2.82, lambda_d = 1.01),
  elderly = list(v_p = 2.95, t_p = 5.1,  lambda_g = 3.15, lambda_d = 0.87)
)


#' Default Ct (cycle threshold) trajectory parameters by role
#'
#' A named list with role-specific parameters for simulating Ct value
#' trajectories. Each role contains: \code{Cpeak} (peak Ct value),
#' \code{r} (rise rate), \code{d} (decay rate), \code{t_peak} (time to peak).
#'
#' @format Named list with elements \code{adult}, \code{infant}, \code{toddler}, \code{elderly}.
#' @keywords internal
#' @noRd
default_Ct_params <- list(
  infant  = list(Cpeak = 33.3, r = 2.11, d = 1.38, t_peak = 5.06),
  toddler = list(Cpeak = 34,   r = 1.26, d = 1.27, t_peak = 4.75),
  adult   = list(Cpeak = 33,   r = 1.49, d = 1.22, t_peak = 5.14),
  elderly = list(Cpeak = 33,   r = 1.49, d = 1.22, t_peak = 5.14)
)


#' Default household composition profile
#'
#' A named list specifying probabilities for household composition:
#' \describe{
#'   \item{prob_adults}{Length-3 probability vector for 0/1/2 adults}
#'   \item{prob_infant}{Probability (0-1) of having an infant}
#'   \item{prob_siblings}{Length-3 probability vector for 0/1/2 toddlers}
#'   \item{prob_elderly}{Length-3 probability vector for 0/1/2 elderly}
#' }
#'
#' @format Named list with probability vectors.
#' @keywords internal
#' @noRd
default_household_profile <- list(
  prob_adults   = c(0.0, 0.0, 1.0),
  prob_infant   = 1.0,
  prob_siblings = c(0.10, 0.70, 0.20),
  prob_elderly  = c(0.70, 0.15, 0.15)
)


#' Default priors for Stan model
#'
#' @keywords internal
#' @noRd
default_priors <- list(
  beta1 = list(dist = "normal", params = c(-5, 1)),
  beta2 = list(dist = "normal", params = c(-5, 1)),
  alpha = list(dist = "normal", params = c(-6, 2)),
  covariates = list(dist = "normal", params = c(0, 1)),
  gen_shape = list(dist = "lognormal", params = c(1.5, 0.5)),
  gen_rate = list(dist = "lognormal", params = c(0.0, 0.5)),
  ct50 = list(dist = "normal", params = c(35.0, 3.0)),
  slope = list(dist = "lognormal", params = c(0.4, 0.5))
)


# ==============================================================================
# VIRAL DATA IMPUTATION
# ==============================================================================

#' Fill Missing Viral Data (Ct or Log10) based on Episode Start
#'
#' This function imputes missing viral data during confirmed episodes using
#' theoretical trajectories defined by the parameters.
#'
#' @param df Dataframe containing 'person_id', 'episode_id', 'date', 'role_name', and the viral column.
#' @param viral_col_name String. Name of the column containing viral data (e.g. "ct_value").
#' @param type String. Either "ct_value" or "log10".
#' @param params_list List of parameters for the curves (role-specific).
#' @param detection_limit Numeric. Value to assign for non-infected days.
#'
#' @return The original dataframe with NAs in the viral column filled.
#' @keywords internal
#' @noRd
fill_missing_viral_data <- function(df, viral_col_name, type = c("ct_value", "log10"),
                                    params_list, detection_limit) {

  type <- match.arg(type)

  if (!viral_col_name %in% names(df)) stop(paste("Column", viral_col_name, "not found in dataframe."))

  viral_sym <- rlang::sym(viral_col_name)

  df_imputed <- df %>%
    dplyr::group_by(person_id, episode_id) %>%
    dplyr::mutate(
      episode_start = min(date, na.rm = TRUE),
      days_since_inf = as.numeric(date - episode_start),
      imputed_val = dplyr::case_when(
        !is.na(!!viral_sym) ~ as.numeric(!!viral_sym),
        episode_id == 0     ~ detection_limit,
        TRUE ~ {
          role <- unique(role_name)[1]
          if (is.na(role)) role <- "adult"
          p <- params_list[[role]]
          if (is.null(p)) p <- params_list[["adult"]]

          if (type == "ct_value") {
            theo <- .get_theoretical_ct(days_since_inf, p)
            pmin(detection_limit, pmax(10, theo))
          } else {
            theo <- .get_theoretical_vl(days_since_inf, p)
            pmax(detection_limit, theo)
          }
        }
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(!!viral_col_name := imputed_val) %>%
    dplyr::select(-episode_start, -days_since_inf, -imputed_val)

  return(df_imputed)
}
