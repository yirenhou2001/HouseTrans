#' Internal globals for NSE checks
#'
#' This file declares symbols used in non-standard evaluation (NSE),
#' e.g., in dplyr/data.table expressions, so R CMD check doesn't flag
#' "no visible binding for global variable" NOTES.
#'
#' @keywords internal
#' @noRd
utils::globalVariables(c(

  # Common NSE placeholders
  ".", "..cols", "..keep", ".SD", "pp",

  # Household/person identifiers and roles
  "HH", "hh_id", "person_id", "role", "role_name",
  "individual_ID", "participantid", "familyidstars",
  "role_id","value",

  # Time/day/test columns
  "day", "day_index", "test_date", "sample_date", "date",
  "day_idx", "bin_start_idx", "date_formatted", "date_bin",

  # Infection/detection timeline fields
  "infection_time", "infectious_start", "infectious_end",
  "infection_resolved", "infection_status", "is_in_episode",
  "detection_time", "first_pos_date", "last_pos_date",
  "date_infection", "date_infectious_start", "date_infectious_end",
  "date_resolved", "date_resolved_final", "prev_resolved_final",
  "start_risk_date", "start_risk", "i_idx",

  # Testing result fields
  "test_result", "pcr_sample", "ct_value",

  # Episode tracking (reinfection support)
  "episode_id", "type", "first_inf",

  # Stan-related / modeling fields
  "stan_id", "hh_factor", "hh_id_int", "p_id_int",
  "time_since_infection", "hh_num", "p_num",

  # Aggregation / SAR / summary fields
  "if_infection", "n_infections", "n_total", "n_infected", "sar",
  "n_primaries", "N_pop", "Primary_AR", "Primary_AR_Pct",
  "n_unique_people", "n_total_episodes", "n_reinfections",
  "Reinf_Rate_Among_Infected", "Reinf_Rate_Pct",

  # Viral load / plotting fields
  "index_vl", "vl_category", "event_type", "vl_test",
  "imputed_val", "episode_start", "days_since_inf",
  "val", "log_value", "group", "med",

  # Plotting fields
  "role_group", "cases", "valid", "total",
  "covariate", "y_plot", "label", "count_in_role", "rank_in_role",
  "target_x", "target_y", "source_x", "source_y", "source_role",
  "lbl_x", "lbl_y", "prob", "alpha_val",
  "x1", "y1", "x2", "y2", "ndx", "ndy", "ndist",
  "nmx", "nmy", "noffset", "nshift_x", "nshift_y",
  "nfinal_x", "nfinal_y", "lbl_x_num",
  "ny1", "ny2", "nx1", "nx2",

  # Transmission chain reconstruction
  "source", "target", "source_int", "n_id", "p_id",

  # Covariate fields
  "vacc_status",

  # Additional fields from collaborator's code
  "resolved_time", "infectious_end_history", "immunity_end_history"
))
