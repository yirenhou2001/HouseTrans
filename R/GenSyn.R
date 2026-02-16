#' Simulate Household Transmission and Estimate Parameters
#'
#' Simulates household infection dynamics with viral load trajectories,
#' reinfections, and covariate effects, then estimates transmission
#' parameters using Bayesian inference via Stan.
#'
#' @param n_households Integer; number of households to simulate.
#' @param start_date,end_date Character or Date; study period bounds.
#' @param surveillance_df Optional data frame with 'date' and 'cases' columns
#'   for seasonal forcing. If provided, overrides \code{seasonal_forcing_list}.
#' @param seasonal_forcing_list Optional named list with seasonal forcing vectors
#'   for each role (adult, infant, toddler, elderly).
#' @param covariates_config Optional list of covariate configurations for simulation.
#'   Each element should be a list with: \code{name} (column name), \code{efficacy}
#'   (effect size, 0-1), \code{effect_on} ("susceptibility", "infectivity", or "both"),
#'   and \code{coverage} (list of probabilities by role).
#' @param alpha_comm_by_role Numeric; baseline community acquisition rate.
#' @param beta1,beta2 Numeric; transmission coefficients.
#' @param delta Numeric; household size scaling exponent.
#' @param phi_by_role Named numeric vector; susceptibility multipliers by role.
#' @param kappa_by_role Named numeric vector; infectivity multipliers by role.
#' @param infectious_shape,infectious_scale Numeric; gamma parameters for infectious period.
#' @param waning_shape,waning_scale Numeric; gamma parameters for immunity waning period.
#' @param peak_day,width Numeric; infectivity profile parameters.
#' @param detect_threshold_log10 Numeric; viral load detection threshold (log10).
#' @param detect_threshold_Ct Numeric; Ct value detection threshold.
#' @param surveillance_interval Integer; days between scheduled tests.
#' @param test_daily Logical; switch to daily testing after first detection.
#' @param viral_testing Character; \code{"viral load"} or \code{"Ct"}.
#' @param V_ref,V_rho Numeric; viral load reference and power for transmission.
#' @param Ct_50,Ct_delta Numeric; Ct-based infectivity parameters.
#' @param VL_params_list Named list; role-specific viral load trajectory parameters.
#' @param Ct_params_list Named list; role-specific Ct trajectory parameters.
#' @param household_profile_list Named list; household composition probabilities.
#' @param max_infections Integer; maximum infections per person (for reinfection modeling).
#' @param use_vl_data Logical; whether to use viral load data in estimation.
#' @param priors Named list; prior specifications for Stan model. Each element
#'   should be a list with \code{dist} ("normal", "uniform", "lognormal") and
#'   \code{params} (parameter vector). Available priors: beta1, beta2, alpha,
#'   covariates, gen_shape, gen_rate, ct50, slope.
#' @param covariates_susceptibility Character vector; names of susceptibility covariates.
#' @param covariates_infectivity Character vector; names of infectivity covariates.
#' @param recovery_params Named list; Gamma parameters (shape, scale) for immunity
#'   tail duration by role.
#' @param imputation_params Named list; parameters for viral curve imputation by role.
#' @param stan_chains,stan_iter,stan_warmup Integers; Stan sampling controls.
#' @param stan_control List; Stan control parameters.
#' @param stan_cores Integer; number of CPU cores for Stan.
#' @param stan_refresh Integer; refresh rate for Stan output.
#' @param seed Integer; random seed for reproducibility.
#'
#' @return An object of class \code{"GenSynResult"} containing:
#' \describe{
#'   \item{$call}{The matched call}
#'   \item{$n_households}{Number of households simulated}
#'   \item{$simulation}{Raw simulation output with hh_df and diagnostic_df}
#'   \item{$stan_data}{Data prepared for Stan}
#'   \item{$fit}{The stanfit object}
#'   \item{$postprocessing}{Tidy posterior summary}
#'   \item{$attack_rates}{Attack rate and reinfection summaries}
#'   \item{$transmission_chains}{Reconstructed transmission links}
#' }
#' @examples
#' \dontrun{
#' # Basic simulation
#' result <- GenSyn(
#'   n_households = 50,
#'   start_date = "2024-07-01",
#'   end_date = "2025-06-30",
#'   stan_chains = 2,
#'   stan_iter = 1000,
#'   stan_warmup = 500,
#'   seed = 123
#' )
#'
#' # View results
#' print(result)
#' plot(result, which = "posterior")
#' }
#' @seealso \code{\link{TransmissionChainAnalysis}}, \code{\link{plot.GenSynResult}}
#' @export
GenSyn <- function(
    n_households = 50,
    start_date = "2024-07-01",
    end_date = "2025-06-30",
    surveillance_df = NULL,
    seasonal_forcing_list = NULL,
    covariates_config = NULL,
    alpha_comm_by_role = 5e-4,
    beta1 = 8e-3,
    beta2 = 8e-3,
    delta = 0,
    phi_by_role = c(adult = 1, infant = 4, toddler = 5, elderly = 1),
    kappa_by_role = c(adult = 1, infant = 1, toddler = 1.2, elderly = 1),
    infectious_shape = 3,
    infectious_scale = 1,
    waning_shape = 16,
    waning_scale = 10,
    peak_day = 1,
    width = 4,
    detect_threshold_log10 = 1e-6,
    detect_threshold_Ct = 99,
    surveillance_interval = 1,
    test_daily = FALSE,
    viral_testing = "viral load",
    V_ref = 3.0,
    V_rho = 2.5,
    Ct_50 = 40,
    Ct_delta = 2,
    VL_params_list = NULL,
    Ct_params_list = NULL,
    household_profile_list = NULL,
    max_infections = Inf,
    use_vl_data = TRUE,
    priors = NULL,
    covariates_susceptibility = NULL,
    covariates_infectivity = NULL,
    recovery_params = NULL,
    imputation_params = NULL,
    stan_chains = 4,
    stan_iter = 2000,
    stan_warmup = 1000,
    stan_control = list(adapt_delta = 0.95, max_treedepth = 15),
    stan_cores = 4,
    stan_refresh = 50,
    seed = NULL
) {

  call <- match.call()
  if (!is.null(seed)) set.seed(seed)
  if (is.null(priors)) priors <- default_priors

  # Simulate households
  sim_result <- simulate_multiple_households_comm(
    n_households = n_households,
    surveillance_df = surveillance_df,
    start_date = start_date,
    end_date = end_date,
    alpha_comm_by_role = alpha_comm_by_role,
    beta1 = beta1, beta2 = beta2, delta = delta,
    phi_by_role = phi_by_role, kappa_by_role = kappa_by_role,
    infectious_shape = infectious_shape, infectious_scale = infectious_scale,
    waning_shape = waning_shape, waning_scale = waning_scale,
    peak_day = peak_day, width = width,
    seasonal_forcing_list = seasonal_forcing_list,
    detect_threshold_log10 = detect_threshold_log10,
    detect_threshold_Ct = detect_threshold_Ct,
    surveillance_interval = surveillance_interval,
    test_daily = test_daily, viral_testing = viral_testing,
    V_ref = V_ref, V_rho = V_rho, Ct_50 = Ct_50, Ct_delta = Ct_delta,
    VL_params_list = VL_params_list, Ct_params_list = Ct_params_list,
    household_profile_list = household_profile_list,
    covariates_config = covariates_config,
    seed = seed, max_infections = max_infections
  )

  # Calculate attack rates
  attack_rates <- tryCatch(summarize_attack_rates(sim_result), error = function(e) NULL)

  # Prepare data for Stan
  df_for_stan <- sim_result$diagnostic_df
  if (!is.null(covariates_config)) {
    cov_names <- sapply(covariates_config, function(x) x$name)
    person_covariates <- sim_result$hh_df %>%
      dplyr::select(dplyr::any_of(c("hh_id", "person_id", cov_names))) %>%
      dplyr::distinct()
    df_for_stan <- df_for_stan %>%
      dplyr::left_join(person_covariates, by = c("hh_id", "person_id"))
  }

  stan_data <- prepare_stan_data(
    df_clean = df_for_stan,
    surveillance_df = surveillance_df,
    study_start_date = as.Date(start_date),
    study_end_date = as.Date(end_date),
    use_vl_data = use_vl_data,
    covariates_susceptibility = covariates_susceptibility,
    covariates_infectivity = covariates_infectivity,
    recovery_params = recovery_params,
    imputation_params = imputation_params,
    priors = priors,
    seed = seed
  )

  # Fit model
  fit <- fit_household_model(
    stan_data = stan_data,
    iter = stan_iter, chains = stan_chains, warmup = stan_warmup,
    control = stan_control, cores = stan_cores, refresh = stan_refresh
  )

  # Postprocess
  posterior_summary <- postprocess_stan_fit(fit)
  chains <- tryCatch(
    reconstruct_transmission_chains(fit, stan_data, min_prob_threshold = 0.01),
    error = function(e) NULL
  )

  out <- list(
    call = call, n_households = n_households,
    simulation = sim_result, surveillance_df = surveillance_df,
    start_date = start_date, end_date = end_date,
    stan_data = stan_data, fit = fit,
    postprocessing = posterior_summary, attack_rates = attack_rates,
    transmission_chains = chains
  )
  class(out) <- "GenSynResult"
  out
}


#' Print method for GenSynResult
#' @param x A \code{GenSynResult} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' \dontrun{
#' # Basic simulation
#' result <- GenSyn(
#'   n_households = 50,
#'   start_date = "2024-07-01",
#'   end_date = "2025-06-30",
#'   stan_chains = 2,
#'   stan_iter = 1000,
#'   stan_warmup = 500,
#'   seed = 123
#' )
#'
#' # View results
#' print(result)
#' }
#' @export
#' @method print GenSynResult
print.GenSynResult <- function(x, ...) {
  cat("GenSyn Result\n")
  cat("=============\n\n")
  cat("Households simulated:", x$n_households, "\n\n")

  if (!is.null(x$attack_rates)) {
    cat("--- Primary Attack Rate ---\n")
    print(x$attack_rates$primary_overall)
    cat("\n")
  }

  if (!is.null(x$postprocessing) && nrow(x$postprocessing) > 0) {
    cat("--- Posterior Summary (Key Parameters) ---\n")
    key_params <- c("beta1", "beta2", "alpha_comm", "phi_by_role", "kappa_by_role")
    subset <- x$postprocessing[grepl(paste(key_params, collapse = "|"), x$postprocessing$Parameter), ]
    if (nrow(subset) > 0) print(head(subset, 20)) else print(head(x$postprocessing, 10))
    cat("\n")
  }

  cat("Access components: $simulation, $fit, $postprocessing, $attack_rates, $transmission_chains\n")
  cat("Use plot(result, which = '...') for visualizations.\n")

  invisible(x)
}


#' Plot method for GenSynResult
#'
#' @param x A \code{GenSynResult} object.
#' @param which Character vector specifying which plots to generate. Options:
#'   \code{"posterior"}, \code{"covariate_effects"}, \code{"epidemic_curve"},
#'   \code{"transmission_chains"}, or \code{"all"}.
#' @param print Logical; whether to print plots immediately.
#' @param hh_id Integer; household ID for transmission chain plot (default: 1).
#' @param prob_cutoff Numeric; minimum probability threshold for chain links (default: 0).
#' @param ... Additional arguments (unused).
#' @return A named list of ggplot objects (invisibly if \code{print = TRUE}).
#' @examples
#' \dontrun{
#' # Basic simulation
#' result <- GenSyn(
#'   n_households = 50,
#'   start_date = "2024-07-01",
#'   end_date = "2025-06-30",
#'   stan_chains = 2,
#'   stan_iter = 1000,
#'   stan_warmup = 500,
#'   seed = 123
#' )
#'
#' # View results
#' plot(result, which = "posterior")
#' }
#' @export
#' @method plot GenSynResult
plot.GenSynResult <- function(x, which = "posterior", print = TRUE,
                              hh_id = 1, prob_cutoff = 0, ...) {

  if ("all" %in% which) {
    which <- c("posterior", "covariate_effects", "epidemic_curve", "transmission_chains")
  }

  plots <- list()

  if ("posterior" %in% which && !is.null(x$fit)) {
    plots$posterior <- tryCatch(
      plot_posterior_distributions(x$fit),
      error = function(e) .empty_plot(e$message)
    )
  }

  if ("covariate_effects" %in% which && !is.null(x$fit) && !is.null(x$stan_data)) {
    if (x$stan_data$K_susc > 0 || x$stan_data$K_inf > 0) {
      plots$covariate_effects <- tryCatch(
        plot_covariate_effects(x$fit, x$stan_data),
        error = function(e) .empty_plot(e$message)
      )
    }
  }

  if ("epidemic_curve" %in% which && !is.null(x$simulation) && !is.null(x$surveillance_df)) {
    plots$epidemic_curve <- tryCatch(
      plot_epidemic_curve(x$simulation, x$surveillance_df, x$start_date),
      error = function(e) .empty_plot(e$message)
    )
  }

  if ("transmission_chains" %in% which && !is.null(x$transmission_chains) && !is.null(x$stan_data)) {
    plots$transmission_chains <- tryCatch(
      plot_household_timeline(x$transmission_chains, x$stan_data,
                              target_hh_id = hh_id, start_date_str = x$start_date,
                              prob_cutoff = prob_cutoff),
      error = function(e) .empty_plot(e$message)
    )
  }

  if (print && length(plots) > 0) {
    for (nm in names(plots)) {
      cat("Plotting:", nm, "\n")
      print(plots[[nm]])
    }
  }

  invisible(plots)
}
