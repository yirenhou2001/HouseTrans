#' Simulate households (RSV/VL engine)
#'
#' Generates synthetic household data for the RSV/VL–Stan pipeline.
#' Returns a list with stacked data frame, per-household list, and
#' long testing records (diagnostic_df).
#'
#' @param n_households Integer; number of households.
#' @param start_date,end_date \code{Date} study window.
#' @param max_days Integer; simulation horizon (days).
#' @param seasonal_forcing_list Optional named list (\code{adult, child, elderly, toddler})
#'   for seasonal forcing.
#' @param alpha_comm_by_role Numeric; baseline community acquisition rate.
#' @param beta1,beta2 Numeric; transmission coefficients.
#' @param delta Numeric; household size scaling exponent.
#' @param phi_by_role Named numeric vector; susceptibility multipliers by role.
#' @param kappa_by_role Named numeric vector; infectivity multipliers by role.
#' @param latent_shape,latent_scale Numeric; gamma parameters for latent period.
#' @param infectious_shape,infectious_scale Numeric; gamma parameters for infectious period.
#' @param resolve_shape,resolve_scale Numeric; gamma parameters for resolution period.
#' @param peak_day,width Numeric; infectivity profile parameters.
#' @param ptrans_threshold Numeric; transmission potential threshold.
#' @param detect_threshold_log10 Numeric; viral load detection threshold in log10.
#' @param detect_threshold_Ct Numeric; Ct value detection threshold.
#' @param surveillance_interval Integer; days between scheduled tests.
#' @param test_daily Logical; switch to daily testing after first detection.
#' @param viral_testing Character; one of \code{"viral load"} or \code{"Ct"}.
#' @param V_ref,V_rho Numeric; viral load reference and power for transmission.
#' @param Ct_50,Ct_delta Numeric; Ct-based infectivity parameters.
#' @param VL_params_list Named list; role-specific VL trajectory parameters.
#' @param Ct_params_list Named list; role-specific Ct trajectory parameters.
#' @param household_profile_list Named list; household composition probabilities.
#' @param verbose Logical; print progress information.
#'
#' @return A list with \code{hh_df} (stacked data frame), \code{households}
#'   (per-household list), and \code{diagnostic_df} (long testing records).
#' @export
simulate_households <- function(
    n_households,
    start_date = as.Date("2024-01-01"),
    end_date   = as.Date("2024-12-31"),
    max_days = 365,
    seasonal_forcing_list = NULL,
    alpha_comm_by_role = 5e-3,
    beta1 = 0.3,
    beta2 = 0.05,
    delta = 0,
    phi_by_role = c(adult = 1.0, child = 7.0, toddler = 7.0, elderly = 4.0),
    kappa_by_role = c(adult = 1.0, child = 1.5, toddler = 1.5, elderly = 1.0),
    latent_shape = 3,
    latent_scale = 0.5,
    infectious_shape = 3,
    infectious_scale = 1,
    resolve_shape = 1.5,
    resolve_scale = 0.5,
    peak_day = 1,
    width = 4,
    ptrans_threshold = 0.5,
    detect_threshold_log10 = 1,
    detect_threshold_Ct = 40,
    surveillance_interval = 7,
    test_daily = FALSE,
    viral_testing = "viral load",
    V_ref = 3,
    V_rho = 2.5,
    Ct_50 = 40,
    Ct_delta = 2,
    VL_params_list = default_VL_params,
    Ct_params_list = default_Ct_params,
    household_profile_list = default_household_profile,
    verbose = FALSE
) {
  simulate_multiple_households_comm(
    n_households = n_households,
    alpha_comm_by_role = alpha_comm_by_role,
    beta1 = beta1,
    beta2 = beta2,
    delta = delta,
    phi_by_role = phi_by_role,
    kappa_by_role = kappa_by_role,
    latent_shape = latent_shape,
    latent_scale = latent_scale,
    infectious_shape = infectious_shape,
    infectious_scale = infectious_scale,
    resolve_shape = resolve_shape,
    resolve_scale = resolve_scale,
    peak_day = peak_day,
    width = width,
    max_days = max_days,
    verbose = verbose,
    seasonal_forcing_list = seasonal_forcing_list,
    ptrans_threshold = ptrans_threshold,
    detect_threshold_log10 = detect_threshold_log10,
    detect_threshold_Ct = detect_threshold_Ct,
    surveillance_interval = surveillance_interval,
    test_daily = test_daily,
    viral_testing = viral_testing,
    V_ref = V_ref,
    V_rho = V_rho,
    Ct_50 = Ct_50,
    Ct_delta = Ct_delta,
    start_date = as.character(start_date),
    VL_params_list = VL_params_list,
    Ct_params_list = Ct_params_list,
    household_profile_list = household_profile_list
  )
}


#' Simulate a single household with community and within-household transmission
#'
#' Simulates infection dynamics for one household over a specified time horizon,
#' incorporating community acquisition, within-household transmission based on
#' viral load or Ct values, and testing/detection logic.
#'
#' @param hh_id Character or integer; household identifier.
#' @param roles Character vector; roles for each household member.
#' @param alpha_comm_by_role Numeric; baseline community acquisition rate.
#' @param beta1,beta2 Numeric; transmission coefficients.
#' @param delta Numeric; household size scaling exponent.
#' @param phi_by_role Named numeric vector; susceptibility multipliers by role.
#' @param kappa_by_role Named numeric vector; infectivity multipliers by role.
#' @param latent_shape,latent_scale Numeric; gamma parameters for latent period.
#' @param infectious_shape,infectious_scale Numeric; gamma parameters for infectious period.
#' @param resolve_shape,resolve_scale Numeric; gamma parameters for resolution period.
#' @param peak_day,width Numeric; infectivity profile parameters.
#' @param max_days Integer; simulation horizon (days).
#' @param perfect_detection Logical; assume perfect detection when VL/Ct exceeds threshold.
#' @param contact_mat Optional matrix; contact structure within household.
#' @param verbose Logical; print debug information.
#' @param seasonal_forcing_list Named list; seasonal forcing vectors by role.
#' @param ptrans_threshold Numeric; transmission potential threshold.
#' @param detect_threshold_log10 Numeric; VL detection threshold (log10).
#' @param detect_threshold_Ct Numeric; Ct detection threshold.
#' @param surveillance_interval Integer; days between scheduled tests.
#' @param test_daily Logical; switch to daily testing after first detection.
#' @param viral_testing Character; \code{"viral load"} or \code{"Ct"}.
#' @param V_ref,V_rho Numeric; viral load reference and power for transmission.
#' @param Ct_50,Ct_delta Numeric; Ct-based infectivity parameters.
#' @param VL_params_input Named list; role-specific VL parameters.
#' @param Ct_params_input Named list; role-specific Ct parameters.
#'
#' @return A data frame with one row per person containing infection timing,
#'   detection, viral load trajectories, and testing results.
#' @keywords internal
simulate_one_household_comm <- function(
    hh_id,
    roles,
    alpha_comm_by_role,
    beta1 = 0.3,
    beta2 = 0.05,
    delta = 0,
    phi_by_role = c(adult = 1.0, child = 7.0, toddler = 7.0, elderly = 4.0),
    kappa_by_role = c(adult = 1.0, child = 1.5, toddler = 1.5, elderly = 1.0),
    latent_shape = 3,
    latent_scale = 0.5,
    infectious_shape = 3,
    infectious_scale = 1,
    resolve_shape = 1.5,
    resolve_scale = 0.5,
    peak_day = 1,
    width = 4,
    max_days = 365,
    perfect_detection = TRUE,
    contact_mat = NULL,
    verbose = FALSE,
    seasonal_forcing_list = NULL,
    ptrans_threshold = 0.5,
    detect_threshold_log10 = 1,
    detect_threshold_Ct = 40,
    surveillance_interval = 7,
    test_daily = FALSE,
    viral_testing = "viral load",
    V_ref = 3,
    V_rho = 2.5,
    Ct_50 = 40,
    Ct_delta = 2,
    VL_params_input = default_VL_params,
    Ct_params_input = default_Ct_params
) {

  n <- length(roles)
  infection_time <- rep(NA_integer_, n)
  infectious_start <- rep(NA_integer_, n)
  infectious_end <- rep(NA_integer_, n)
  detection_time <- rep(NA_integer_, n)
  infection_resolved <- rep(NA_integer_, n)
  viral_loads_log10 <- vector("list", n)
  vl_trajs <- vector("list", n)

  latent_period <- pmax(1, ceiling(rgamma(n, shape = latent_shape, scale = latent_scale)))
  infectious_period <- pmax(1, ceiling(rgamma(n, shape = infectious_shape, scale = infectious_scale)))
  resolve_period <- pmax(1, ceiling(rgamma(n, shape = resolve_shape, scale = resolve_scale)))

  if (is.null(contact_mat)) {
    contact_mat <- matrix(1, n, n)
    diag(contact_mat) <- 0
  }

  if (is.null(seasonal_forcing_list)) {
    seasonal_forcing_list <- list(
      adult = rep(0.1, max_days),
      child = rep(0.1, max_days),
      elderly = rep(0.1, max_days),
      toddler = rep(0.1, max_days)
    )
  }

  household_detected <- FALSE
  index_case_id <- NA_integer_

  g_profile <- exp(-0.5 * ((1:max_days - peak_day) / width)^2)
  g_profile <- g_profile / max(g_profile)

  for (t in 1:max_days) {
    is_scheduled_day <- ((t - 1) %% surveillance_interval == 0)

    if (!household_detected) {
      test_today <- is_scheduled_day
    } else {
      test_today <- if (test_daily) TRUE else is_scheduled_day
    }

    p_jt_vec <- rep(0, n)

    for (j in seq_len(n)) {
      if (!is.na(infection_time[j]) && infection_time[j] < t) next

      lambda_house_j <- 0
      for (i in seq_len(n)) {
        if (i == j) next
        if (is.na(infection_time[i])) next
        if (is.na(infectious_start[i]) || is.na(infectious_end[i])) next
        if (t < infectious_start[i] || t > infectious_end[i]) next

        tau <- t - infectious_start[i] + 1
        gval <- if (tau >= 1 && tau <= length(g_profile)) 1 else 0
        rel_day_idx <- t - infection_time[i] + 1

        if (viral_testing == "viral load") {
          log10V <- if (!is.na(infection_time[i]) && rel_day_idx >= 1 &&
                        !is.null(vl_trajs[[i]]) && rel_day_idx <= length(vl_trajs[[i]])) {
            vl_trajs[[i]][rel_day_idx]
          } else { 0 }

          V_it <- max(0, ifelse(is.na(log10V), 0, log10V))
          scaling_n <- (1 / max(n, 1))^delta
          term1 <- beta1 * gval
          term2 <- beta2 * (V_it / V_ref)^V_rho
          h_ij_t <- scaling_n * kappa_by_role[roles[i]] * (term1 + term2)
          lambda_house_j <- lambda_house_j + h_ij_t * contact_mat[j, i]

        } else if (viral_testing == "Ct") {
          log10V <- if (!is.na(infection_time[i]) && rel_day_idx >= 1 &&
                        !is.null(vl_trajs[[i]]) && rel_day_idx <= length(vl_trajs[[i]])) {
            vl_trajs[[i]][rel_day_idx]
          } else { 45 }

          V_it <- ifelse(is.na(log10V), 45, log10V)
          scaling_n <- (1 / max(n, 1))^delta
          term1 <- beta1 * gval
          term2 <- beta2 * 1 / (1 + exp((V_it - Ct_50) / Ct_delta))
          h_ij_t <- scaling_n * kappa_by_role[roles[i]] * (term1 + term2)
          lambda_house_j <- lambda_house_j + h_ij_t * contact_mat[j, i]
        }
      }

      alpha_base_j <- alpha_comm_by_role
      season_mult_j <- seasonal_forcing_list[[roles[j]]]
      alpha_comm_j_baseline <- alpha_base_j * season_mult_j[min(t, length(season_mult_j))]
      lambda_jt <- phi_by_role[roles[j]] * (alpha_comm_j_baseline + lambda_house_j)
      lambda_jt <- pmin(pmax(lambda_jt, 0), 1e6)
      p_jt <- 1 - exp(-lambda_jt)
      p_jt <- pmin(pmax(p_jt, 1e-20), 1 - 1e-20)

      if (is.na(infection_time[j])) p_jt_vec[j] <- p_jt
    }

    new_infections <- runif(n) < p_jt_vec
    for (j in which(new_infections)) {
      if (is.na(index_case_id)) index_case_id <- j
      infection_time[j] <- t
      t_seq <- seq(1, 20)
      t_rel <- t_seq - 1

      if (viral_testing == "viral load") {
        p <- draw_random_VL_params(roles[j], VL_params_input)
        vl_trajs[[j]] <- simulate_viral_load_trajectory(t_rel, p$v_p, p$t_p, p$lambda_g, p$lambda_d)
        if (!is.null(vl_trajs[[j]]) && length(vl_trajs[[j]]) > 0) {
          ptrans_vec <- (vl_trajs[[j]] / V_ref)^V_rho
          inf_from_vl_idx <- which(ptrans_vec > ptrans_threshold)
        }
      } else {
        p <- draw_random_Ct_params(roles[j], Ct_params_input)
        vl_trajs[[j]] <- simulate_Ct_trajectory(t_rel, p$Cpeak, p$r, p$d, p$t_peak)
        if (!is.null(vl_trajs[[j]]) && length(vl_trajs[[j]]) > 0) {
          infectivity_vec <- 1 / (1 + exp((vl_trajs[[j]] - Ct_50) / Ct_delta))
          inf_from_vl_idx <- which(infectivity_vec > ptrans_threshold)
        }
      }

      if (exists("inf_from_vl_idx") && length(inf_from_vl_idx) > 0) {
        inf_start_from_vl <- infection_time[j] + min(inf_from_vl_idx) - 1
        inf_end_from_vl <- infection_time[j] + max(inf_from_vl_idx) - 1
        infectious_start[j] <- pmax(1, pmin(inf_start_from_vl, max_days))
        infectious_end[j] <- pmax(1, pmin(inf_end_from_vl, max_days))
      } else {
        infectious_start[j] <- NA_integer_
        infectious_end[j] <- NA_integer_
      }

      if (!is.na(infectious_end[j])) {
        infection_resolved[j] <- pmin(infectious_end[j] + resolve_period[j], max_days)
      } else {
        infection_resolved[j] <- pmin(infection_time[j] + resolve_period[j], max_days)
      }
    }

    if (test_today) {
      for (i in seq_len(n)) {
        if (is.na(infection_time[i])) next
        vl_rel_idx <- t - infection_time[i] + 1
        val_today <- NA_real_
        if (!is.null(vl_trajs[[i]]) && vl_rel_idx >= 1 && vl_rel_idx <= length(vl_trajs[[i]])) {
          val_today <- vl_trajs[[i]][vl_rel_idx]
        }

        if (is.na(detection_time[i]) && perfect_detection && !is.na(val_today)) {
          is_positive <- FALSE
          if (viral_testing == "viral load") {
            if (val_today >= detect_threshold_log10) is_positive <- TRUE
          } else {
            if (val_today <= detect_threshold_Ct) is_positive <- TRUE
          }
          if (is_positive) {
            detection_time[i] <- t
            household_detected <- TRUE
          }
        }
      }
    }
  }

  # Generate test days for output
  first_det <- if (all(is.na(detection_time))) NA_integer_ else min(detection_time, na.rm = TRUE)
  test_days <- integer(0)

  for (d in seq_len(max_days)) {
    is_scheduled <- ((d - 1) %% surveillance_interval == 0)
    is_post_detection <- !is.na(first_det) && d >= first_det

    if (is_post_detection) {
      if (test_daily) {
        test_days <- c(test_days, d)
      } else {
        if (is_scheduled) test_days <- c(test_days, d)
      }
    } else {
      if (is_scheduled) test_days <- c(test_days, d)
    }
  }

  for (i in seq_len(n)) {
    if (is.null(vl_trajs[[i]]) || all(is.na(vl_trajs[[i]]))) {
      viral_loads_log10[[i]] <- if (viral_testing == "viral load") 0 else 45
    } else {
      vals <- sapply(test_days, function(dd) {
        if (is.na(infection_time[i]) || dd < infection_time[i] || dd > infection_resolved[i]) {
          return(if (viral_testing == "viral load") 0 else 45)
        }
        rel_day <- dd - infection_time[i] + 1
        if (rel_day >= 1 && rel_day <= length(vl_trajs[[i]])) return(vl_trajs[[i]][rel_day])
        else return(if (viral_testing == "viral load") 0 else 45)
      })
      names(vals) <- as.character(test_days)
      viral_loads_log10[[i]] <- vals
    }
  }

  hh_df <- data.frame(
    hh_id = hh_id,
    person_id = seq_len(n),
    role = roles,
    infection_time = ifelse(is.infinite(infection_time), NA_integer_, infection_time),
    infectious_start = ifelse(is.infinite(infectious_start), NA_integer_, infectious_start),
    infectious_end = ifelse(is.infinite(infectious_end), NA_integer_, infectious_end),
    detection_time = ifelse(is.infinite(detection_time), NA_integer_, detection_time),
    infection_resolved = ifelse(is.infinite(infection_resolved), NA_integer_, infection_resolved),
    stringsAsFactors = FALSE
  )

  hh_df$vl_full_trajectory <- vl_trajs
  hh_df$viral_loads_test_days <- viral_loads_log10
  attr(hh_df, "test_days") <- test_days
  attr(hh_df, "params") <- list(beta1 = beta1, beta2 = beta2, delta = delta, V_rho = V_rho, V_ref = V_ref)
  hh_df
}


#' Simulate multiple households with community transmission
#'
#' Wrapper that calls \code{\link{simulate_one_household_comm}} for each household
#' and assembles results into a combined data frame and diagnostic testing records.
#'
#' @inheritParams simulate_one_household_comm
#' @param n_households Integer; number of households to simulate.
#' @param start_date Character or Date; start date for actual calendar dates.
#' @param VL_params_list Named list; role-specific VL parameters.
#' @param Ct_params_list Named list; role-specific Ct parameters.
#' @param household_profile_list Named list; household composition probabilities.
#'
#' @return A list with \code{hh_df}, \code{households}, and \code{diagnostic_df}.
#' @keywords internal
simulate_multiple_households_comm <- function(
    n_households,
    alpha_comm_by_role = 5e-3,
    beta1 = 0.3,
    beta2 = 0.05,
    delta = 0,
    phi_by_role = c(adult = 1.0, child = 7.0, toddler = 7.0, elderly = 4.0),
    kappa_by_role = c(adult = 1.0, child = 1.5, toddler = 1.5, elderly = 1.0),
    latent_shape = 3,
    latent_scale = 0.5,
    infectious_shape = 3,
    infectious_scale = 1,
    resolve_shape = 1.5,
    resolve_scale = 0.5,
    peak_day = 1,
    width = 4,
    max_days = 365,
    verbose = FALSE,
    seasonal_forcing_list = NULL,
    ptrans_threshold = 0.5,
    detect_threshold_log10 = 1,
    detect_threshold_Ct = 40,
    surveillance_interval = 7,
    test_daily = FALSE,
    viral_testing = "viral load",
    V_ref = 3,
    V_rho = 2.5,
    Ct_50 = 40,
    Ct_delta = 2,
    start_date = "2024-01-01",
    VL_params_list = default_VL_params,
    Ct_params_list = default_Ct_params,
    household_profile_list = default_household_profile
) {

  start_date_obj <- as.Date(start_date)
  households <- vector("list", n_households)

  for (h in seq_len(n_households)) {
    roles <- generate_household_roles(household_profile_list)

    hh <- simulate_one_household_comm(
      hh_id = paste0("HH", h),
      roles = roles,
      alpha_comm_by_role = alpha_comm_by_role,
      beta1 = beta1, beta2 = beta2, delta = delta,
      phi_by_role = phi_by_role, kappa_by_role = kappa_by_role,
      latent_shape = latent_shape, latent_scale = latent_scale,
      infectious_shape = infectious_shape, infectious_scale = infectious_scale,
      resolve_shape = resolve_shape, resolve_scale = resolve_scale,
      peak_day = peak_day, width = width, max_days = max_days,
      verbose = verbose, seasonal_forcing_list = seasonal_forcing_list,
      ptrans_threshold = ptrans_threshold,
      detect_threshold_log10 = detect_threshold_log10,
      detect_threshold_Ct = detect_threshold_Ct,
      surveillance_interval = surveillance_interval,
      test_daily = test_daily, viral_testing = viral_testing,
      V_ref = V_ref, V_rho = V_rho, Ct_50 = Ct_50, Ct_delta = Ct_delta,
      VL_params_input = VL_params_list, Ct_params_input = Ct_params_list
    )
    households[[h]] <- hh
  }

  hh_df <- do.call(rbind, households)
  rownames(hh_df) <- NULL

  # Generate diagnostic DataFrame
  all_records <- list()
  counter <- 1

  for (h_idx in seq_along(households)) {
    hh <- households[[h_idx]]
    hh_id_val <- hh$hh_id[1]
    test_days <- attr(hh, "test_days")

    if (length(test_days) == 0) next

    n_people <- nrow(hh)
    for (i in 1:n_people) {
      p_id <- hh$person_id[i]
      role <- hh$role[i]
      raw_vals <- hh$viral_loads_test_days[[i]]

      sample_vals <- if (length(raw_vals) == length(test_days)) raw_vals else rep(raw_vals[1], length(test_days))

      test_result <- if (viral_testing == "viral load") {
        as.integer(sample_vals >= detect_threshold_log10)
      } else {
        as.integer(sample_vals <= detect_threshold_Ct)
      }

      actual_dates <- start_date_obj + (test_days - 1)

      person_df <- data.frame(
        hh_id = hh_id_val,
        person_id = p_id,
        role = role,
        sample_date = actual_dates,
        day_index = test_days,
        pcr_sample = as.numeric(sample_vals),
        test_result = test_result,
        stringsAsFactors = FALSE
      )

      all_records[[counter]] <- person_df
      counter <- counter + 1
    }
  }

  diagnostic_df <- if (length(all_records) > 0) {
    df <- do.call(rbind, all_records)
    rownames(df) <- NULL
    df
  } else {
    data.frame()
  }

  list(hh_df = hh_df, households = households, diagnostic_df = diagnostic_df)
}


#' Normalize various household inputs into a per-household list
#'
#' Accepts (i) an object with \code{$households}, (ii) a list of data frames,
#' or (iii) a long data frame with a recognizable household id column, and
#' returns a list of per-household data frames.
#'
#' @param obj Mixed household representation.
#' @return List of per-household data frames.
#' @keywords internal
.normalize_households_input <- function(obj) {
  hh_attr <- attr(obj, "households", exact = TRUE)
  if (!is.null(hh_attr)) {
    if (!length(hh_attr) || !all(vapply(hh_attr, is.data.frame, logical(1)))) {
      stop("Attribute 'households' exists but is not a list of data.frames.", call. = FALSE)
    }
    return(hh_attr)
  }

  if (is.list(obj) && !is.data.frame(obj) && !is.null(obj$households)) {
    hh <- obj$households
    if (!length(hh) || !all(vapply(hh, is.data.frame, logical(1)))) {
      stop("`$households` exists but is not a list of data.frames.", call. = FALSE)
    }
    return(hh)
  }

  if (is.list(obj) && length(obj) && all(vapply(obj, is.data.frame, logical(1)))) {
    return(obj)
  }

  if (inherits(obj, c("data.frame", "tbl_df", "tbl", "data.table"))) {
    df <- as.data.frame(obj)
    candidates <- c("hh_id", "HH", "household", "household_id", "householdID", "hh", "hid")
    key <- candidates[candidates %in% names(df)]
    if (length(key) == 0L) {
      stop("Cannot find a household id column. Please include one of: hh_id, HH, household, household_id.", call. = FALSE)
    }
    key <- key[1]
    df$hh_id <- as.character(df[[key]])
    households <- split(df, df$hh_id, drop = TRUE)
    names(households) <- unique(df$hh_id)
    return(households)
  }

  stop("Expected (i) list(hh_df=..., households=...), or (ii) a list of per-household data.frames, or (iii) a long data.frame with a household id column.", call. = FALSE)
}
