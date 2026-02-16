#' @title Household Transmission Simulator
#' @description Core simulation engine for household infection dynamics
#' @name simulator
#' @keywords internal
NULL


# ==============================================================================
# 1. CORE SIMULATION ENGINE (Single Household)
# ==============================================================================

#' Simulate a Single Household
#'
#' @param hh_id Household identifier.
#' @param roles Character vector of roles for each person.
#' @param alpha_comm_by_role Community acquisition rate.
#' @param beta1,beta2 Transmission coefficients.
#' @param delta Household size scaling exponent.
#' @param phi_by_role,kappa_by_role Susceptibility/infectivity by role.
#' @param infectious_shape,infectious_scale Gamma parameters for infectious period.
#' @param waning_shape,waning_scale Gamma parameters for immunity waning.
#' @param peak_day,width Infectivity profile parameters.
#' @param max_days Maximum simulation days.
#' @param test_weekly_before_detection Logical.
#' @param perfect_detection Logical.
#' @param contact_mat Optional contact matrix.
#' @param verbose Logical.
#' @param seasonal_forcing_list Seasonal forcing by role.
#' @param detect_threshold_log10 Detection threshold (log10 VL).
#' @param detect_threshold_Ct Detection threshold (Ct).
#' @param surveillance_interval Days between tests.
#' @param test_daily Test daily after detection.
#' @param viral_testing "viral load" or "Ct".
#' @param V_ref,V_rho VL parameters.
#' @param Ct_50,Ct_delta Ct parameters.
#' @param VL_params_input,Ct_params_input Trajectory parameters.
#' @param susc_modifiers_vec Susceptibility modifiers.
#' @param inf_modifiers_vec Infectivity modifiers.
#' @param covariate_data Covariate data frame.
#' @param max_infections Maximum infections per person.
#'
#' @return List with hh_df and diagnostic_df.
#' @keywords internal
simulate_one_household_comm <- function(hh_id,
                                        roles,
                                        alpha_comm_by_role,
                                        beta1, beta2, delta,
                                        phi_by_role, kappa_by_role,
                                        infectious_shape, infectious_scale,
                                        waning_shape, waning_scale,
                                        peak_day, width,
                                        max_days,
                                        test_weekly_before_detection,
                                        perfect_detection,
                                        contact_mat = NULL,
                                        verbose,
                                        seasonal_forcing_list,
                                        detect_threshold_log10,
                                        detect_threshold_Ct,
                                        surveillance_interval,
                                        test_daily,
                                        viral_testing,
                                        V_ref, V_rho,
                                        Ct_50, Ct_delta,
                                        VL_params_input,
                                        Ct_params_input,
                                        susc_modifiers_vec = NULL,
                                        inf_modifiers_vec = NULL,
                                        covariate_data = NULL,
                                        max_infections = Inf) {

  n <- length(roles)

  infection_history       <- vector("list", n)
  infectious_end_history  <- vector("list", n)
  immunity_end_history    <- vector("list", n)

  infection_counts  <- integer(n)
  current_status    <- integer(n) # 0=S, 1=I, 2=R
  time_next_state   <- rep(NA_integer_, n)

  current_vl_traj   <- vector("list", n)
  current_inf_start <- rep(NA_integer_, n)
  detection_time    <- rep(NA_integer_, n)

  # Defaults
  if (is.null(susc_modifiers_vec)) susc_modifiers_vec <- rep(1.0, n)
  if (is.null(inf_modifiers_vec))  inf_modifiers_vec  <- rep(1.0, n)

  scaling_n <- (1.0 / max(n, 1))^delta
  phi_vec   <- phi_by_role[roles]
  kappa_vec <- kappa_by_role[roles]

  if (is.null(contact_mat)) {
    contact_mat <- matrix(1, n, n)
    diag(contact_mat) <- 0
  }

  household_detected <- FALSE

  # ==========================
  # MAIN TIME LOOP
  # ==========================
  for (t in 1:max_days) {

    # --- A. Update States ---
    for (i in 1:n) {
      if (current_status[i] != 0 && !is.na(time_next_state[i]) && t >= time_next_state[i]) {
        st <- current_status[i]
        if (st == 1) { # I -> R (Recovery)
          current_status[i] <- 2
          infectious_end_history[[i]] <- c(infectious_end_history[[i]], t)
          dur <- pmax(1, ceiling(stats::rgamma(1, shape = waning_shape, scale = waning_scale)))
          time_next_state[i] <- t + dur
          immunity_end_history[[i]] <- c(immunity_end_history[[i]], t + dur)
          current_vl_traj[i] <- list(NULL)
          current_inf_start[i] <- NA
        } else if (st == 2) { # R -> S
          current_status[i] <- 0
          time_next_state[i] <- NA
        }
      }
    }

    is_scheduled_day <- ((t - 1) %% surveillance_interval == 0)
    test_today <- if (household_detected && test_daily) TRUE else is_scheduled_day

    # --- B. Calculate Force ---
    total_hh_force <- rep(0.0, n)
    infectors <- which(current_status == 1)

    if (length(infectors) > 0) {
      infectivity_values <- numeric(length(infectors))
      for (k in seq_along(infectors)) {
        idx <- infectors[k]
        rel_day <- t - current_inf_start[idx] + 1
        val <- NA_real_
        traj <- current_vl_traj[[idx]]
        if (!is.null(traj) && rel_day >= 1 && rel_day <= length(traj)) val <- traj[rel_day]

        term1 <- beta1 * 1.0
        term2 <- 0.0
        if (!is.na(val)) {
          if (viral_testing == "viral load") {
            val_clean <- max(0, val)
            term2 <- beta2 * (val_clean / V_ref)^V_rho
          } else {
            val_clean <- ifelse(val > 45, 45, val)
            term2 <- beta2 * 1.0 / (1.0 + exp((val_clean - Ct_50) / Ct_delta))
          }
        }
        # Apply Infectivity Modifier
        infectivity_values[k] <- scaling_n * kappa_vec[idx] * inf_modifiers_vec[idx] * (term1 + term2)
      }
      for (k in seq_along(infectors)) {
        src <- infectors[k]
        force <- infectivity_values[k]
        contacts <- which(contact_mat[, src] > 0)
        total_hh_force[contacts] <- total_hh_force[contacts] + (force * contact_mat[contacts, src])
      }
    }

    # --- C. New Infections ---
    targets <- which(current_status == 0 & infection_counts < max_infections)

    if (length(targets) > 0) {
      alpha_comm_val <- numeric(length(targets))
      for (k in seq_along(targets)) {
        tgt <- targets[k]
        role_name <- roles[tgt]
        season_val <- seasonal_forcing_list[[role_name]][t]
        alpha_comm_val[k] <- alpha_comm_by_role * season_val
      }

      # Apply Susceptibility Modifier
      lambda_vec <- (phi_vec[targets] * susc_modifiers_vec[targets]) * (alpha_comm_val + total_hh_force[targets])
      lambda_vec <- pmin(lambda_vec, 1e6)
      prob_inf <- 1.0 - exp(-lambda_vec)

      is_infected <- stats::runif(length(targets)) < prob_inf
      infected_indices <- targets[is_infected]

      for (j in infected_indices) {
        infection_history[[j]] <- c(infection_history[[j]], t)
        infection_counts[j] <- infection_counts[j] + 1
        current_status[j] <- 1 # I

        dur_I <- pmax(1, ceiling(stats::rgamma(1, shape = infectious_shape, scale = infectious_scale)))
        time_next_state[j] <- t + dur_I
        current_inf_start[j] <- t

        # Curve length is exactly infection duration
        t_seq <- 0:dur_I

        if (viral_testing == "viral load") {
          p <- draw_random_VL_params(roles[j], VL_params_input)
          traj <- simulate_viral_load_trajectory(t_seq, p$v_p, p$t_p, p$lambda_g, p$lambda_d)
        } else {
          p <- draw_random_Ct_params(roles[j], Ct_params_input)
          traj <- simulate_Ct_trajectory(t_seq, p$Cpeak, p$r, p$d, p$t_peak)
        }
        current_vl_traj[[j]] <- traj
      }
    }

    # --- D. Testing ---
    if (test_today) {
      shedders <- which(current_status == 1)
      for (i in shedders) {
        rel_day <- t - current_inf_start[i] + 1
        traj <- current_vl_traj[[i]]
        if (!is.null(traj) && rel_day >= 1 && rel_day <= length(traj)) {
          val <- traj[rel_day]
          is_pos <- FALSE
          if (viral_testing == "viral load") {
            if (val >= detect_threshold_log10) is_pos <- TRUE
          } else {
            if (val <= detect_threshold_Ct) is_pos <- TRUE
          }
          if (is_pos && is.na(detection_time[i])) {
            detection_time[i] <- t
            household_detected <- TRUE
          }
        }
      }
    }
  }

  # ==============================================================================
  # EXPORT RESULTS
  # ==============================================================================

  hh_rows <- list()
  for (i in 1:n) {
    hist <- infection_history[[i]]

    # Base Row
    base_row <- data.frame(
      hh_id = hh_id, person_id = i, role = roles[i],
      stringsAsFactors = FALSE
    )

    # Merge covariates (remove duplicate person_id)
    if (!is.null(covariate_data)) {
      covs_clean <- covariate_data[i, !names(covariate_data) %in% "person_id", drop = FALSE]
      base_row <- cbind(base_row, covs_clean)
    }

    if (is.null(hist)) {
      base_row$infection_time <- NA
      base_row$infectious_end <- NA
      base_row$resolved_time  <- NA
      hh_rows[[length(hh_rows) + 1]] <- base_row
    } else {
      for (k in seq_along(hist)) {
        row <- base_row
        row$infection_time <- hist[k]
        row$infectious_end <- if (k <= length(infectious_end_history[[i]])) infectious_end_history[[i]][k] else NA
        row$resolved_time  <- if (k <= length(immunity_end_history[[i]])) immunity_end_history[[i]][k] else NA
        hh_rows[[length(hh_rows) + 1]] <- row
      }
    }
  }
  hh_df <- do.call(rbind, hh_rows)

  # --- Diagnostic DF ---
  test_days <- seq(1, max_days, by = surveillance_interval)
  if (any(!is.na(detection_time)) && test_daily) {
    first_det <- min(detection_time, na.rm = TRUE)
    if (!is.infinite(first_det)) test_days <- unique(sort(c(test_days, seq(first_det, max_days))))
  }
  n_tests <- length(test_days)
  res_matrix <- matrix(if (viral_testing == "viral load") 0 else 45, nrow = n, ncol = n_tests)
  ep_matrix <- matrix(0, nrow = n, ncol = n_tests)

  for (i in 1:n) {
    hist <- infection_history[[i]]
    if (!is.null(hist)) {
      for (k in seq_along(hist)) {
        inf_t   <- hist[k]
        inf_end <- if (k <= length(infectious_end_history[[i]])) infectious_end_history[[i]][k] else (inf_t + 10)
        p <- if (viral_testing == "viral load") draw_random_VL_params(roles[i], VL_params_input) else draw_random_Ct_params(roles[i], Ct_params_input)

        # Stop window exactly at infection end
        relevant_indices <- which(test_days >= inf_t & test_days <= inf_end)

        if (length(relevant_indices) > 0) {
          d_vals <- test_days[relevant_indices]
          rel_vals <- d_vals - inf_t + 1
          if (viral_testing == "viral load") {
            v_vals <- simulate_viral_load_trajectory(rel_vals, p$v_p, p$t_p, p$lambda_g, p$lambda_d)
            res_matrix[i, relevant_indices] <- pmax(res_matrix[i, relevant_indices], v_vals)
          } else {
            ct_vals <- simulate_Ct_trajectory(rel_vals, p$Cpeak, p$r, p$d, p$t_peak)
            res_matrix[i, relevant_indices] <- pmin(res_matrix[i, relevant_indices], ct_vals)
          }
        }

        if (length(relevant_indices) > 0) ep_matrix[i, relevant_indices] <- k
      }
    }
  }

  diag_df_list <- vector("list", n)
  for (i in 1:n) {
    vals <- res_matrix[i, ]
    eps  <- ep_matrix[i, ]
    results <- if (viral_testing == "viral load") as.integer(vals >= detect_threshold_log10) else as.integer(vals <= detect_threshold_Ct)
    diag_df_list[[i]] <- data.frame(
      hh_id = hh_id, person_id = i, role = roles[i],
      day_index = test_days, pcr_sample = vals, test_result = results,
      episode_id = eps,
      stringsAsFactors = FALSE
    )
  }
  diagnostic_df <- do.call(rbind, diag_df_list)

  list(hh_df = hh_df, diagnostic_df = diagnostic_df)
}


# ==============================================================================
# 2. EXPORTED WRAPPER (Multiple Households with Generic Covariates)
# ==============================================================================

#' Simulate Multiple Households with Community Transmission
#'
#' Simulates household infection dynamics across multiple households with
#' support for generic covariates, reinfections, and viral load trajectories.
#'
#' @param n_households Integer; number of households to simulate.
#' @param surveillance_df Optional data frame with 'date' and 'cases' columns.
#' @param start_date,end_date Character or Date; study period bounds.
#' @param alpha_comm_by_role Numeric; community acquisition rate.
#' @param beta1,beta2 Numeric; transmission coefficients.
#' @param delta Numeric; household size scaling exponent.
#' @param phi_by_role,kappa_by_role Named numeric vectors.
#' @param infectious_shape,infectious_scale Numeric; Gamma parameters.
#' @param waning_shape,waning_scale Numeric; Gamma parameters.
#' @param peak_day,width Numeric; infectivity profile.
#' @param verbose Logical; print progress.
#' @param seasonal_forcing_list Named list of forcing vectors.
#' @param detect_threshold_log10 Numeric; detection threshold (log10 VL).
#' @param detect_threshold_Ct Numeric; detection threshold (Ct).
#' @param surveillance_interval Integer; days between tests.
#' @param test_daily Logical; test daily after detection.
#' @param viral_testing Character; "viral load" or "Ct".
#' @param V_ref,V_rho Numeric; VL parameters.
#' @param Ct_50,Ct_delta Numeric; Ct parameters.
#' @param VL_params_list,Ct_params_list Named lists of trajectory parameters.
#' @param household_profile_list Named list for household composition.
#' @param perfect_detection Logical.
#' @param contact_mat Optional contact matrix.
#' @param covariates_config List of covariate configurations.
#' @param seed Integer; random seed.
#' @param max_infections Integer; max infections per person.
#'
#' @return List with hh_df and diagnostic_df.
#' @examples
#' \dontrun{
#' household_profile <- list(
#'   prob_adults   = c(0, 0, 1),
#'   prob_infant   = 1.0,
#'   prob_siblings = c(0, .8, .2),
#'   prob_elderly  = c(0.7, 0.1, 0.2)
#' )
#'
#' sim_res <- simulate_multiple_households_comm(
#'   n_households = 100,
#'   viral_testing = "viral load",
#'   infectious_shape = 10,
#'   infectious_scale = 1,
#'   waning_shape = 6,
#'   waning_scale = 10,
#'   surveillance_interval = 4,
#'   start_date = “2024-07-01”,
#'   end_date = “2025-06-30”,
#'   surveillance_df = surveillance_data,
#'   seed = 123,
#'   household_profile_list = household_profile
#' )
#' }
#' @export
simulate_multiple_households_comm <- function(n_households = 50,
                                              surveillance_df = NULL,
                                              start_date = "2024-07-01",
                                              end_date = "2025-06-30",
                                              alpha_comm_by_role = 5e-4,
                                              beta1 = 8e-3, beta2 = 8e-3, delta = 0,
                                              phi_by_role = c(adult = 1, infant = 4, toddler = 5, elderly = 1),
                                              kappa_by_role = c(adult = 1, infant = 1, toddler = 1.2, elderly = 1),
                                              infectious_shape = 3, infectious_scale = 1,
                                              waning_shape = 16, waning_scale = 10,
                                              peak_day = 1, width = 4,
                                              verbose = FALSE,
                                              seasonal_forcing_list = NULL,
                                              detect_threshold_log10 = 1e-6,
                                              detect_threshold_Ct = 99,
                                              surveillance_interval = 1,
                                              test_daily = FALSE,
                                              viral_testing = "viral load",
                                              V_ref = 3.0, V_rho = 2.5,
                                              Ct_50 = 40, Ct_delta = 2,
                                              VL_params_list = NULL,
                                              Ct_params_list = NULL,
                                              household_profile_list = NULL,
                                              perfect_detection = TRUE,
                                              contact_mat = NULL,
                                              covariates_config = NULL,
                                              seed = NULL,
                                              max_infections = Inf) {

  if (!is.null(seed)) set.seed(seed)

  # Date logic
  d_start <- as.Date(start_date)
  d_end <- as.Date(end_date)
  max_days <- as.integer(d_end - d_start) + 1
  if (max_days <= 0) stop("end_date must be after start_date")

  # Seasonality
  final_forcing_list <- NULL
  if (!is.null(surveillance_df)) {
    target_dates <- seq(from = d_start, to = d_end, by = "day")
    interp_res <- stats::approx(x = as.Date(surveillance_df$date), y = surveillance_df$cases, xout = target_dates, rule = 2)
    daily_vec <- interp_res$y
    if (max(daily_vec, na.rm = TRUE) > 0) daily_vec <- daily_vec / max(daily_vec, na.rm = TRUE)
    daily_vec[is.na(daily_vec)] <- 0
    final_forcing_list <- list(adult = daily_vec, infant = daily_vec, toddler = daily_vec, elderly = daily_vec)
  } else if (!is.null(seasonal_forcing_list)) {
    final_forcing_list <- seasonal_forcing_list
  }
  if (is.null(final_forcing_list)) {
    final_forcing_list <- list(
      adult = rep(1, max_days), infant = rep(1, max_days),
      toddler = rep(1, max_days), elderly = rep(1, max_days)
    )
  }

  if (is.null(household_profile_list)) {
    household_profile_list <- list(
      prob_single_parent = 0,
      prob_siblings = c(0.10, 0.50, 0.40),
      prob_elderly = c(0.9, 0.08, 0.02)
    )
  }
  if (is.null(VL_params_list)) VL_params_list <- default_VL_params
  if (is.null(Ct_params_list)) Ct_params_list <- default_Ct_params

  households <- vector("list", n_households)

  for (h in seq_len(n_households)) {
    roles <- generate_household_roles(household_profile_list)
    n_hh <- length(roles)

    # --- GENERATE COVARIATES ---
    hh_covariates <- data.frame(person_id = 1:n_hh)

    # Init Modifiers
    hh_susc_modifiers <- rep(1.0, n_hh)
    hh_inf_modifiers  <- rep(1.0, n_hh)

    if (!is.null(covariates_config)) {
      for (cov in covariates_config) {
        # Determine Status
        col_vec <- numeric(n_hh)
        for (i in 1:n_hh) {
          prob <- cov$coverage[[roles[i]]]
          if (is.null(prob)) prob <- 0
          col_vec[i] <- stats::rbinom(1, 1, prob)
        }
        hh_covariates[[cov$name]] <- col_vec

        # Apply Efficacy
        eff_type <- tolower(cov$effect_on)
        multiplier <- (1.0 - (col_vec * cov$efficacy))

        if (eff_type %in% c("susceptibility", "both")) hh_susc_modifiers <- hh_susc_modifiers * multiplier
        if (eff_type %in% c("infectivity", "both")) hh_inf_modifiers <- hh_inf_modifiers * multiplier
      }
    }

    hh <- simulate_one_household_comm(
      hh_id = paste0("HH", h),
      roles = roles,
      alpha_comm_by_role = alpha_comm_by_role,
      beta1 = beta1, beta2 = beta2, delta = delta,
      phi_by_role = phi_by_role, kappa_by_role = kappa_by_role,
      infectious_shape = infectious_shape, infectious_scale = infectious_scale,
      waning_shape = waning_shape, waning_scale = waning_scale,
      peak_day = peak_day, width = width,
      max_days = max_days,
      verbose = verbose,
      seasonal_forcing_list = final_forcing_list,
      detect_threshold_log10 = detect_threshold_log10,
      detect_threshold_Ct = detect_threshold_Ct,
      surveillance_interval = surveillance_interval,
      test_daily = test_daily,
      viral_testing = viral_testing,
      V_ref = V_ref, V_rho = V_rho,
      Ct_50 = Ct_50, Ct_delta = Ct_delta,
      VL_params_input = VL_params_list,
      Ct_params_input = Ct_params_list,
      perfect_detection = perfect_detection,
      contact_mat = contact_mat,
      test_weekly_before_detection = TRUE,
      susc_modifiers_vec = hh_susc_modifiers,
      inf_modifiers_vec  = hh_inf_modifiers,
      covariate_data = hh_covariates,
      max_infections = max_infections
    )
    households[[h]] <- hh
  }

  hh_list <- lapply(households, function(x) x$hh_df)
  hh_df <- dplyr::bind_rows(hh_list)

  diag_list <- lapply(households, function(x) x$diagnostic_df)
  diag_list <- diag_list[!sapply(diag_list, is.null)]
  diagnostic_df <- if (length(diag_list) > 0) do.call(rbind, diag_list) else data.frame()

  list(hh_df = hh_df, diagnostic_df = diagnostic_df)
}
