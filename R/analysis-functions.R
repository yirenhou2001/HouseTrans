#' @title Analysis Functions
#' @description Functions for summarizing simulation results and attack rates
#' @name analysis_functions
#' @keywords internal
NULL


#' Summarize Attack Rates & Reinfections
#'
#' Calculates the Primary Attack Rate (proportion of people infected at least once)
#' and a separate summary of Reinfections (secondary episodes).
#'
#' @param sim_result A list object returned by \code{simulate_multiple_households_comm}.
#'
#' @return A list containing four dataframes:
#' \describe{
#'   \item{primary_overall}{Overall primary attack rate statistics.}
#'   \item{primary_by_role}{Primary attack rate stratified by age group.}
#'   \item{reinf_overall}{Overall summary of reinfection counts and rates.}
#'   \item{reinf_by_role}{Reinfection counts stratified by age group.}
#' }
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
#' rates <- summarize_attack_rates(sim_res)
#' }
#' @export
summarize_attack_rates <- function(sim_result) {

  if (!"hh_df" %in% names(sim_result)) stop("Input must contain 'hh_df'.")
  df <- sim_result$hh_df

  # ---------------------------------------------------------
  # 1. PREPARE DATA
  # ---------------------------------------------------------

  # A. Total Population (Denominator)
  pop_stats <- df %>%
    dplyr::select(hh_id, person_id, role) %>%
    dplyr::distinct() %>%
    dplyr::group_by(role) %>%
    dplyr::summarise(N_pop = dplyr::n(), .groups = "drop")

  total_pop <- sum(pop_stats$N_pop)

  # B. Infection Episodes (Numerator)
  infected_df <- df %>%
    dplyr::filter(!is.na(infection_time)) %>%
    dplyr::select(hh_id, person_id, role, infection_time)

  # ---------------------------------------------------------
  # 2. PRIMARY ATTACK RATE (First Infection Only)
  # ---------------------------------------------------------

  # Count UNIQUE people infected (at least once)
  primary_counts_by_role <- infected_df %>%
    dplyr::group_by(role) %>%
    dplyr::summarise(n_primaries = dplyr::n_distinct(paste(hh_id, person_id)), .groups = "drop")

  # Merge with Population
  primary_by_role <- pop_stats %>%
    dplyr::left_join(primary_counts_by_role, by = "role") %>%
    dplyr::mutate(
      n_primaries = ifelse(is.na(n_primaries), 0, n_primaries),
      Primary_AR  = n_primaries / N_pop,
      Primary_AR_Pct = scales::percent(Primary_AR, accuracy = 0.1)
    )

  # Overall Primary Stats
  total_primaries <- sum(primary_by_role$n_primaries)
  primary_overall <- data.frame(
    Total_Pop = total_pop,
    Total_Infected_People = total_primaries,
    Primary_Attack_Rate = total_primaries / total_pop,
    Primary_Attack_Rate_Pct = scales::percent(total_primaries / total_pop, accuracy = 0.1)
  )

  # ---------------------------------------------------------
  # 3. REINFECTION SUMMARY (Subsequent Infections)
  # ---------------------------------------------------------

  # Calculate Total Episodes vs Unique People
  reinf_counts_by_role <- infected_df %>%
    dplyr::group_by(role) %>%
    dplyr::summarise(
      n_total_episodes = dplyr::n(),
      n_unique_people  = dplyr::n_distinct(paste(hh_id, person_id)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      n_reinfections = n_total_episodes - n_unique_people,
      Reinf_Rate_Among_Infected = ifelse(n_unique_people > 0, n_reinfections / n_unique_people, 0)
    )

  # Merge to ensure all roles appear even if 0 reinfections
  reinf_by_role <- pop_stats %>%
    dplyr::select(role) %>%
    dplyr::left_join(reinf_counts_by_role, by = "role") %>%
    dplyr::mutate(
      n_reinfections = ifelse(is.na(n_reinfections), 0, n_reinfections),
      Reinf_Rate_Pct = scales::percent(ifelse(is.na(Reinf_Rate_Among_Infected), 0, Reinf_Rate_Among_Infected), accuracy = 0.1)
    ) %>%
    dplyr::select(role, n_unique_people, n_total_episodes, n_reinfections, Reinf_Rate_Pct)

  # Overall Reinfection Stats
  total_episodes_all <- nrow(infected_df)
  total_reinf_all    <- total_episodes_all - total_primaries

  reinf_overall <- data.frame(
    Total_Primary_Infections = total_primaries,
    Total_Reinfection_Events = total_reinf_all,
    Reinfection_Rate_Among_Infected = ifelse(total_primaries > 0, total_reinf_all / total_primaries, 0)
  )

  # ---------------------------------------------------------
  # 4. RETURN
  # ---------------------------------------------------------

  list(
    primary_overall = primary_overall,
    primary_by_role = primary_by_role,
    reinf_overall   = reinf_overall,
    reinf_by_role   = reinf_by_role
  )
}
