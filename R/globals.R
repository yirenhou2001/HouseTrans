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
  "HH", "hh_id", "person_id", "role",
  "ID_hh", "ID_indiv", "individual_ID", "indiv.index",

  # Time/day/test columns
  "day", "day_index", "test_date", "sample_date",

  # Infection/detection timeline fields
  "infection_time", "infectious_start", "infectious_end",
  "infection_resolved", "infection_status",
  "detection_time",
  "infected", "is_index",

  # Derived infection window / summary variables (from summarize_individuals)
  "inf_date", "inf_start_date", "inf_end_date",
  "inf_win_start", "inf_win_end",
  "T_FP_date", "T_LP_date",
  "infection.detected.start", "infection.detected.end",
  "infection.infectious.day",
  "last_neg_date", "last_negative",

  # Testing result fields
  "test_result", "pcr_sample",

  # Stan-related / modeling fields
  "stan_id", "hh_factor", "time_since_infection",

  # Aggregation / SAR / summary fields
  "if_infection", "n_infections", "n_total", "n_infected", "sar",

  # Viral load / plotting fields
  "index_vl", "vl_category", "event_type", "vl_test",

  # Covariate bins used in running_parameter_estimation (from check output)
  "age_cat", "agegrp2", "agegrp3", "agegrp4"
))
