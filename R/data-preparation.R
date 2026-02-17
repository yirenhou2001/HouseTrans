#' @title Data Preparation for Stan Model
#' @description Functions to prepare household data for Stan inference
#' @name data_preparation
#' @keywords internal
NULL


#' Prepare Data for Stan Model
#'
#' Prepares household observation data for Stan model fitting, including
#' imputation of episode timing, covariate matrices, and viral load data.
#'
#' @param df_clean Dataframe with observation data (must contain 'episode_id').
#' @param surveillance_df Dataframe with columns 'date' and 'cases'.
#' @param role_levels Character vector; canonical role levels.
#' @param study_start_date,study_end_date Date objects; study period.
#' @param seasonal_forcing_list Named list of forcing vectors.
#' @param use_vl_data Logical; whether to use viral load data.
#' @param covariates_susceptibility Character vector; susceptibility covariate names.
#' @param covariates_infectivity Character vector; infectivity covariate names.
#' @param recovery_params List of Gamma parameters (shape, scale) by role.
#' @param imputation_params List of viral curve parameters by role.
#' @param priors List of flexible priors (dist, params).
#' @param seed Integer; random seed.
#'
#' @return A list of data formatted for Stan.
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
#' person_covariates <- sim_res$hh_df %>%
#'   dplyr::select(hh_id, person_id, vacc_status) %>%   dplyr::distinct()
#'
#' df_for_stan <- sim_res$diagnostic_df %>%
#'   dplyr::left_join(person_covariates, by = c("hh_id", "person_id"))
#'
#' stan_input <- prepare_stan_data(
#'   df_clean = df_for_stan,
#'   use_vl_data      = 1,
#'   study_start_date = as.Date(study_start),
#'   study_end_date   = as.Date(study_end),
#'   surveillance_df  = surveillance_data,
#'   study_start_date = as.Date(study_start),
#'   study_end_date   = as.Date(study_end),
#'   imputation_params = VL_params_list)
#' }
#' @export
prepare_stan_data <- function(df_clean,
                              surveillance_df = NULL,
                              role_levels = c("adult", "infant", "toddler", "elderly"),
                              study_start_date = as.Date("2024-07-01"),
                              study_end_date = as.Date("2025-07-01"),
                              seasonal_forcing_list = NULL,
                              use_vl_data = TRUE,
                              covariates_susceptibility = NULL,
                              covariates_infectivity = NULL,
                              recovery_params = NULL,
                              imputation_params = NULL,
                              priors = list(),
                              seed = 123) {

  if (!is.null(seed)) set.seed(seed)
  T_max <- as.integer(study_end_date - study_start_date) + 1
  if (T_max <= 0) stop("study_end_date must be after study_start_date")

  # =========================================================
  # 1. PARSE FLEXIBLE PRIORS
  # =========================================================

  # --- A. Transmission Priors ---
  p_beta1 <- .parse_prior(priors$beta1, 1, c(-5, 1))
  p_beta2 <- .parse_prior(priors$beta2, 1, c(-5, 1))
  p_alpha <- .parse_prior(priors$alpha, 1, c(-6, 2))
  p_cov   <- .parse_prior(priors$covariates, 1, c(0, 1))

  # --- B. Viral Dynamics Priors ---
  p_shape <- .parse_prior(priors$gen_shape, 1, c(5.0, 2.0))
  p_rate  <- .parse_prior(priors$gen_rate, 1, c(1.0, 0.5))
  p_ct50  <- .parse_prior(priors$ct50, 1, c(35.0, 3.0))
  p_slope <- .parse_prior(priors$slope, 1, c(1.5, 1.0))

  # =========================================================
  # 2. STANDARDIZE COLUMN NAMES
  # =========================================================
  if ("hh_id" %in% names(df_clean) && !"familyidstars" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::rename(familyidstars = hh_id)
  }
  if ("person_id" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(participantid = paste(familyidstars, person_id, sep = "_"))
  }
  if ("role" %in% names(df_clean) && !"role_name" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::rename(role_name = role)
  }
  if ("pcr_sample" %in% names(df_clean) && !"ct_value" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::rename(ct_value = pcr_sample)
  }
  if ("test_result" %in% names(df_clean) && !"is_in_episode" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::rename(is_in_episode = test_result)
  }

  # Ensure date column exists
  if ("day_index" %in% names(df_clean) && !"date" %in% names(df_clean)) {
    df_clean <- df_clean %>% dplyr::mutate(date = study_start_date + (day_index - 1))
  }

  # =========================================================
  # 3. IMPUTATION STEP (Interval Sampling + Recovery Tail)
  # =========================================================

  if (!is.null(recovery_params)) {
    params_resolve <- recovery_params
  } else {
    params_resolve <- list(
      adult   = list(shape = 2, scale = 3),
      infant  = list(shape = 2, scale = 3),
      toddler = list(shape = 2, scale = 3),
      elderly = list(shape = 2, scale = 3)
    )
  }

  participants <- unique(df_clean$participantid)
  imputed_episodes <- list()

  for (pid in participants) {
    p_data <- df_clean %>%
      dplyr::filter(participantid == pid) %>%
      dplyr::arrange(date)

    # Get Role for Tail Parameters
    p_role <- p_data$role_name[1]
    if (is.na(p_role) || !p_role %in% role_levels) p_role <- "adult"
    cur_resolve <- params_resolve[[p_role]]
    if (is.null(cur_resolve)) cur_resolve <- list(shape = 2, scale = 3)

    # Get Episode IDs
    episode_ids <- if ("episode_id" %in% names(p_data)) unique(p_data$episode_id[p_data$episode_id > 0]) else integer(0)

    for (eid in episode_ids) {
      ep_rows <- p_data %>% dplyr::filter(episode_id == eid)
      if (nrow(ep_rows) == 0) next

      # --- 1. DETERMINE INFECTIOUS PERIOD (Data-Driven) ---
      first_pos_date <- min(ep_rows$date)
      last_pos_date  <- max(ep_rows$date)

      # Find Lower Bound (Last Negative)
      prev_neg_dates <- p_data %>% dplyr::filter(date < first_pos_date, episode_id == 0) %>% dplyr::pull(date)
      lower_bound <- if (length(prev_neg_dates) > 0) max(prev_neg_dates) else (first_pos_date - 7)

      # Find Upper Bound (Next Negative)
      next_neg_dates <- p_data %>% dplyr::filter(date > last_pos_date, episode_id == 0) %>% dplyr::pull(date)
      upper_bound <- if (length(next_neg_dates) > 0) min(next_neg_dates) else (last_pos_date + 7)

      # Sample Start (Uniformly)
      possible_starts <- seq(from = lower_bound + 1, to = first_pos_date, by = "day")
      T_inf_new <- if (length(possible_starts) > 1) sample(possible_starts, 1) else first_pos_date

      # Sample End of Infectiousness (Uniformly)
      possible_ends <- seq(from = last_pos_date, to = upper_bound - 1, by = "day")
      T_resolved_new <- if (length(possible_ends) > 1) sample(possible_ends, 1) else last_pos_date

      if (T_resolved_new < T_inf_new) T_resolved_new <- T_inf_new

      # --- 2. ADD RECOVERY TAIL (Model-Driven) ---
      resolve_draw <- max(1, round(stats::rgamma(1, shape = cur_resolve$shape, scale = cur_resolve$scale)))
      T_resolved_final <- T_resolved_new + resolve_draw

      # --- 3. PREVENT OVERLAP WITH NEXT EPISODE ---
      if (eid < max(episode_ids)) {
        next_ep_rows <- p_data %>% dplyr::filter(episode_id == (eid + 1))
        if (nrow(next_ep_rows) > 0) {
          next_ep_raw_start <- min(next_ep_rows$date)
          if (T_resolved_final >= next_ep_raw_start) {
            T_resolved_final <- next_ep_raw_start - 1
          }
        }
      }

      imputed_episodes[[length(imputed_episodes) + 1]] <- data.frame(
        participantid = pid, episode_id = eid,
        date_infection = T_inf_new,
        date_infectious_start = T_inf_new,
        date_infectious_end = T_resolved_new,
        date_resolved = T_resolved_new,
        date_resolved_final = T_resolved_final,
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine Imputed Episodes
  if (length(imputed_episodes) > 0) {
    df_imputed <- do.call(rbind, imputed_episodes)
  } else {
    df_imputed <- data.frame(
      participantid = character(), episode_id = integer(),
      date_infection = as.Date(character()), date_infectious_start = as.Date(character()),
      date_infectious_end = as.Date(character()), date_resolved = as.Date(character()),
      date_resolved_final = as.Date(character()), stringsAsFactors = FALSE
    )
  }

  # =========================================================
  # 4. METADATA & MATRIX BUILDING
  # =========================================================

  # Columns to preserve
  cols_to_keep <- c("participantid", "familyidstars", "role_name", "person_id")
  all_covs <- unique(c(covariates_susceptibility, covariates_infectivity))
  if (!is.null(all_covs) && length(all_covs) > 0) {
    missing <- setdiff(all_covs, names(df_clean))
    if (length(missing) > 0) stop(paste("Covariates missing in df_clean:", paste(missing, collapse = ", ")))
    cols_to_keep <- c(cols_to_keep, all_covs)
  }

  # Build Metadata
  df_meta_all <- df_clean %>%
    dplyr::select(dplyr::all_of(cols_to_keep)) %>%
    dplyr::distinct() %>%
    dplyr::rename(hh_id = familyidstars, role = role_name) %>%
    dplyr::filter(!is.na(role)) %>%
    dplyr::mutate(
      hh_num = as.numeric(gsub("\\D", "", hh_id)),
      p_num  = as.numeric(person_id)
    ) %>%
    dplyr::arrange(hh_num, p_num) %>%
    dplyr::select(-hh_num, -p_num)

  # Join Metadata with Imputed Episodes
  df_model_full <- df_meta_all %>%
    dplyr::left_join(df_imputed, by = "participantid")

  # Safety Check
  if (!"date_resolved_final" %in% names(df_model_full)) {
    df_model_full$date_resolved_final <- as.Date(NA)
    df_model_full$episode_id <- NA
    df_model_full$date_infection <- as.Date(NA)
    df_model_full$date_infectious_end <- as.Date(NA)
  }

  # Process Start/End Risk Dates
  df_model_full <- df_model_full %>%
    dplyr::group_by(participantid) %>%
    dplyr::mutate(
      prev_resolved_final = dplyr::lag(date_resolved_final, default = study_start_date - 1),
      start_risk_date = dplyr::case_when(
        is.na(episode_id) | episode_id == 1 ~ study_start_date,
        TRUE ~ prev_resolved_final + 1
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      start_risk = as.integer(start_risk_date - study_start_date) + 1,
      i_idx = dplyr::row_number()
    )

  N <- nrow(df_model_full)
  unique_hh <- unique(df_model_full$hh_id)
  H <- length(unique_hh)

  df_model_full$hh_id_int <- as.integer(factor(df_model_full$hh_id, levels = unique_hh))
  df_model_full$p_id_int  <- as.integer(as.factor(df_model_full$participantid))

  # --- COVARIATE MATRICES ---
  if (!is.null(covariates_susceptibility) && length(covariates_susceptibility) > 0) {
    X_susc <- as.matrix(df_model_full[, covariates_susceptibility, drop = FALSE])
    K_susc <- ncol(X_susc)
  } else {
    X_susc <- matrix(0, nrow = N, ncol = 0)
    K_susc <- 0
  }

  if (!is.null(covariates_infectivity) && length(covariates_infectivity) > 0) {
    X_inf <- as.matrix(df_model_full[, covariates_infectivity, drop = FALSE])
    K_inf <- ncol(X_inf)
  } else {
    X_inf <- matrix(0, nrow = N, ncol = 0)
    K_inf <- 0
  }

  # --- DATA MATRICES (I, Y, V) ---
  max_obs_val <- max(df_clean$ct_value, na.rm = TRUE)
  is_ct_data <- (is.finite(max_obs_val) && max_obs_val > 15)
  detected_vl_type <- if (is_ct_data) 0 else 1
  default_val <- if (is_ct_data) 45.0 else 0.0

  I <- matrix(0L, N, T_max)
  Y <- matrix(0L, N, T_max)
  V <- matrix(default_val, N, T_max)

  for (i in 1:N) {
    row <- df_model_full[i, ]

    if (!is.na(row$episode_id)) {
      idx_inf <- as.integer(row$date_infection - study_start_date) + 1
      idx_end <- as.integer(row$date_infectious_end - study_start_date) + 1

      # 1. Fill Infection Matrices
      if (idx_inf >= 1 && idx_inf <= T_max) I[i, idx_inf] <- 1L

      y_start <- max(1, idx_inf)
      y_end   <- min(T_max, idx_end)
      if (y_end >= y_start) Y[i, y_start:y_end] <- 1L

      # 2. MECHANISTIC IMPUTATION
      p_role <- row$role
      p_params <- imputation_params[[p_role]]

      # Fallback
      if (is.null(p_params)) {
        if (is_ct_data) p_params <- list(Cpeak = 33, r = 1.5, d = 1.2, t_peak = 5)
        else p_params <- list(v_p = 5, t_p = 4, lambda_g = 2.8, lambda_d = 1.0)
      }

      days_seq <- seq(y_start, y_end)

      if (length(days_seq) > 0) {
        t_vals <- days_seq - idx_inf
        mechanistic_vals <- sapply(t_vals, function(t) .calc_curve_val(t, p_params, is_ct_data))
        V[i, days_seq] <- mechanistic_vals

        # 3. RE-APPLY OBSERVED DATA
        obs <- df_clean %>%
          dplyr::filter(participantid == row$participantid, !is.na(ct_value))

        for (k in seq_len(nrow(obs))) {
          d_idx <- as.integer(obs$date[k] - study_start_date) + 1
          if (d_idx >= 1 && d_idx <= T_max) {
            V[i, d_idx] <- obs$ct_value[k]
          }
        }
      }
    }
  }

  # --- HOUSEHOLD STRUCTURE ---
  hh_members_list <- split(1:N, df_model_full$hh_id_int)
  hh_size_eps <- sapply(hh_members_list, length)
  hh_max_size <- max(hh_size_eps)
  hh_members <- matrix(0L, nrow = H, ncol = hh_max_size)
  for (h in 1:H) hh_members[h, 1:hh_size_eps[h]] <- hh_members_list[[h]]

  hh_size_df <- df_clean %>%
    dplyr::group_by(familyidstars) %>%
    dplyr::summarise(n = dplyr::n_distinct(participantid), .groups = "drop")
  hh_size_people <- hh_size_df$n[match(unique_hh, hh_size_df$familyidstars)]

  # --- SEASONALITY ---
  seasonal_forcing_mat <- matrix(1.0, nrow = T_max, ncol = 4)

  if (!is.null(surveillance_df)) {
    if (!all(c("date", "cases") %in% names(surveillance_df))) stop("surveillance_df must have 'date' and 'cases'")
    surveillance_df$date <- as.Date(surveillance_df$date)
    target_dates <- seq(from = study_start_date, to = study_end_date, by = "day")
    interp_res <- stats::approx(x = surveillance_df$date, y = surveillance_df$cases, xout = target_dates, method = "linear", rule = 2)
    daily_vec <- interp_res$y
    if (max(daily_vec, na.rm = TRUE) > 0) daily_vec <- daily_vec / max(daily_vec, na.rm = TRUE)
    daily_vec[is.na(daily_vec)] <- 0
    for (k in 1:4) seasonal_forcing_mat[, k] <- daily_vec
  }

  # =========================================================
  # 5. RETURN LIST
  # =========================================================

  list(
    N = N, T = T_max, H = H, R = length(role_levels), delta = 0,
    hh_id = df_model_full$hh_id_int,
    role_id = match(df_model_full$role, role_levels),
    I = I, Y = Y, V = V,
    start_risk = df_model_full$start_risk,
    p_id = df_model_full$p_id_int,
    hh_size_people = array(as.integer(hh_size_people), dim = H),
    hh_max_size = as.integer(hh_max_size),
    hh_members = hh_members,
    seasonal_forcing_mat = seasonal_forcing_mat,
    reference_phi = 1.0, reference_kappa = 1.0,
    use_vl_data = as.integer(use_vl_data),

    # Detect VL Type (High values > 15 imply Ct, Low values imply Log10)
    vl_type = as.integer(detected_vl_type),
    use_curve_logic = 1L,

    # Covariates
    K_susc = K_susc, X_susc = X_susc,
    K_inf  = K_inf,  X_inf  = X_inf,

    # Flexible Priors (Transmission)
    prior_beta1_type = p_beta1$type, prior_beta1_params = p_beta1$params,
    prior_beta2_type = p_beta2$type, prior_beta2_params = p_beta2$params,
    prior_alpha_type = p_alpha$type, prior_alpha_params = p_alpha$params,
    prior_cov_type   = p_cov$type,   prior_cov_params   = p_cov$params,

    # Flexible Priors (Viral Dynamics)
    prior_shape_type = p_shape$type, prior_shape_params = p_shape$params,
    prior_rate_type  = p_rate$type,  prior_rate_params  = p_rate$params,
    prior_ct50_type  = p_ct50$type,  prior_ct50_params  = p_ct50$params,
    prior_slope_type = p_slope$type, prior_slope_params = p_slope$params
  )
}


#' Build Stan Household Arrays (Per-Person Format)
#'
#' @param households List of household data frames.
#' @param T_max Integer; time horizon.
#' @param seasonal_forcing_list Named list of forcing vectors.
#' @param priors List of priors.
#' @param use_vl_data Logical.
#'
#' @return Stan data list.
#' @keywords internal
build_stan_household_arrays <- function(households,
                                        T_max,
                                        seasonal_forcing_list = NULL,
                                        priors = list(),
                                        use_vl_data = TRUE) {

  role_levels <- c("adult", "infant", "toddler", "elderly")

  # Parse priors
  p_beta1 <- .parse_prior(priors$beta1, 1, c(-5, 1))
  p_beta2 <- .parse_prior(priors$beta2, 1, c(-5, 1))
  p_alpha <- .parse_prior(priors$alpha, 1, c(-6, 2))
  p_cov   <- .parse_prior(priors$covariates, 1, c(0, 1))
  p_shape <- .parse_prior(priors$gen_shape, 1, c(5.0, 2.0))
  p_rate  <- .parse_prior(priors$gen_rate, 1, c(1.0, 0.5))
  p_ct50  <- .parse_prior(priors$ct50, 1, c(35.0, 3.0))
  p_slope <- .parse_prior(priors$slope, 1, c(1.5, 1.0))

  # Flatten households
  all_rows <- dplyr::bind_rows(households, .id = "hh_id_orig")

  N <- nrow(all_rows)
  H <- length(households)

  unique_hh <- unique(all_rows$hh_id)
  all_rows$hh_id_int <- as.integer(factor(all_rows$hh_id, levels = unique_hh))
  all_rows$p_id_int <- 1:N

  # Build matrices
  I <- matrix(0L, N, T_max)
  Y <- matrix(0L, N, T_max)
  V <- matrix(0.0, N, T_max)

  for (i in 1:N) {
    row <- all_rows[i, ]
    if (!is.na(row$infection_time)) {
      t_inf <- as.integer(row$infection_time)
      t_end <- as.integer(row$infectious_end %||% row$infection_resolved %||% (t_inf + 7))

      if (t_inf >= 1 && t_inf <= T_max) I[i, t_inf] <- 1L

      y_start <- max(1, t_inf)
      y_end   <- min(T_max, t_end)
      if (y_end >= y_start) Y[i, y_start:y_end] <- 1L
    }
  }

  # Household structure
  hh_members_list <- split(1:N, all_rows$hh_id_int)
  hh_size_eps <- sapply(hh_members_list, length)
  hh_max_size <- max(hh_size_eps)
  hh_members <- matrix(0L, nrow = H, ncol = hh_max_size)
  for (h in 1:H) hh_members[h, 1:hh_size_eps[h]] <- hh_members_list[[h]]

  # Seasonal forcing
  seasonal_forcing_mat <- matrix(1.0, nrow = T_max, ncol = 4)
  if (!is.null(seasonal_forcing_list)) {
    for (k in 1:4) {
      role_nm <- role_levels[k]
      if (!is.null(seasonal_forcing_list[[role_nm]])) {
        vec <- seasonal_forcing_list[[role_nm]]
        seasonal_forcing_mat[1:min(length(vec), T_max), k] <- vec[1:min(length(vec), T_max)]
      }
    }
  }

  list(
    N = N, T = T_max, H = H, R = 4, delta = 0,
    hh_id = all_rows$hh_id_int,
    role_id = match(all_rows$role, role_levels),
    I = I, Y = Y, V = V,
    start_risk = rep(1L, N),
    p_id = all_rows$p_id_int,
    hh_size_people = array(as.integer(hh_size_eps), dim = H),
    hh_max_size = as.integer(hh_max_size),
    hh_members = hh_members,
    seasonal_forcing_mat = seasonal_forcing_mat,
    reference_phi = 1.0, reference_kappa = 1.0,
    use_vl_data = as.integer(use_vl_data),
    vl_type = 1L,
    use_curve_logic = 1L,
    K_susc = 0, X_susc = matrix(0, nrow = N, ncol = 0),
    K_inf  = 0, X_inf  = matrix(0, nrow = N, ncol = 0),
    prior_beta1_type = p_beta1$type, prior_beta1_params = p_beta1$params,
    prior_beta2_type = p_beta2$type, prior_beta2_params = p_beta2$params,
    prior_alpha_type = p_alpha$type, prior_alpha_params = p_alpha$params,
    prior_cov_type   = p_cov$type,   prior_cov_params   = p_cov$params,
    prior_shape_type = p_shape$type, prior_shape_params = p_shape$params,
    prior_rate_type  = p_rate$type,  prior_rate_params  = p_rate$params,
    prior_ct50_type  = p_ct50$type,  prior_ct50_params  = p_ct50$params,
    prior_slope_type = p_slope$type, prior_slope_params = p_slope$params
  )
}
