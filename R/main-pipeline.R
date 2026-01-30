#' Main parameter estimation pipeline
#'
#' End-to-end workflow for household transmission estimation via Bayesian Stan.
#' Works with simulated or user-provided data.
#'
#' @param user_data Optional; data.frame or list of data.frames.
#' @param synthetic_data Logical; if \code{TRUE}, simulate data internally.
#' @param n_households Integer; number of households to simulate.
#' @param start_date,end_date \code{Date}; analysis window.
#' @param seasonal_forcing_list Optional named list for seasonality.
#' @param max_days Integer; maximum days.
#' @param alpha_comm_by_role Numeric; community acquisition rate.
#' @param beta1,beta2 Numeric; transmission coefficients.
#' @param delta Numeric; household size scaling.
#' @param phi_by_role,kappa_by_role Named numeric vectors.
#' @param latent_shape,latent_scale Numeric; latent period parameters.
#' @param infectious_shape,infectious_scale Numeric; infectious period parameters.
#' @param resolve_shape,resolve_scale Numeric; resolution period parameters.
#' @param peak_day,width Numeric; infectivity profile.
#' @param ptrans_threshold,detect_threshold_log10,detect_threshold_Ct Numeric; testing thresholds.
#' @param surveillance_interval Integer; days between tests.
#' @param test_daily Logical; daily testing after detection.
#' @param viral_testing Character; \code{"viral load"} or \code{"Ct"}.
#' @param V_ref,V_rho Numeric; viral load parameters.
#' @param Ct_50,Ct_delta Numeric; Ct parameters.
#' @param VL_params_list,Ct_params_list,household_profile_list Named lists.
#' @param stan_file Path to Stan model file.
#' @param stan_chains,stan_iter,stan_warmup Integers; Stan sampling controls.
#' @param stan_control List; Stan \code{control} list.
#' @param stan_init Character or function; Stan initialization.
#' @param stan_refresh Integer; Stan refresh.
#' @param stan_cores Integer; CPU cores for Stan.
#'
#' @return A list with \code{raw_simulation}, \code{stan_data}, \code{fit},
#'   \code{posterior_summary}, and \code{diagnostic_df}.
#'
#' @seealso \code{\link{GenSyn}}
#' @export
main_parameter_estimation_pipeline <- function(
    user_data = NULL,
    synthetic_data = TRUE,
    n_households = 10,
    start_date = as.Date("2024-01-01"),
    end_date   = as.Date("2024-12-31"),

    # Seasonal forcing
    seasonal_forcing_list = NULL,
    max_days = 365,

    # Transmission parameters
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

    # Testing parameters
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

    # Stan controls
    stan_file = "model.stan",
    stan_chains = 4, stan_iter = 2000, stan_warmup = 1000,
    stan_control = list(adapt_delta = 0.99, max_treedepth = 20),
    stan_init = "random", stan_refresh = 50, stan_cores = 4
) {

  # Data generation / ingestion
  if (synthetic_data) {
    raw_dt <- simulate_households(
      n_households = n_households,
      start_date = start_date,
      end_date = end_date,
      max_days = max_days,
      seasonal_forcing_list = seasonal_forcing_list,
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
      VL_params_list = VL_params_list,
      Ct_params_list = Ct_params_list,
      household_profile_list = household_profile_list
    )
  } else {
    raw_dt <- user_data
  }

  # Extract households and diagnostic_df
  if (is.list(raw_dt) && !is.null(raw_dt$diagnostic_df)) {
    diagnostic_df <- raw_dt$diagnostic_df
    households <- raw_dt$households
  } else {
    households <- .normalize_households_input(raw_dt)
    diagnostic_df <- NULL
  }

  # Prepare Stan data
  if (!is.null(diagnostic_df) && nrow(diagnostic_df) > 0) {
    stan_data <- prepare_stan_data(
      raw_data = diagnostic_df,
      seasonal_forcing_list = seasonal_forcing_list,
      viral_testing_mode = viral_testing,
      T_max = max_days,
      V_ref = V_ref,
      V_rho = V_rho,
      peak_day = peak_day,
      width = width,
      alpha_comm_by_role = alpha_comm_by_role
    )
  } else {
    stan_data <- build_stan_household_arrays(
      households = households,
      T_max = max_days,
      seasonal_forcing_list = seasonal_forcing_list,
      alpha_comm_by_role = alpha_comm_by_role,
      beta1 = beta1,
      beta2 = beta2,
      V_ref = V_ref
    )
  }

  # Run Stan
  fit <- run_household_stan(
    stan_data,
    stan_file = stan_file,
    chains = stan_chains, iter = stan_iter, warmup = stan_warmup,
    control = stan_control, init = stan_init, refresh = stan_refresh, cores = stan_cores
  )

  posterior_summary <- postprocess_stan_fit(fit)

  list(
    raw_simulation    = raw_dt,
    households        = households,
    diagnostic_df     = diagnostic_df,
    stan_data         = stan_data,
    fit               = fit,
    posterior_summary = posterior_summary,
    postprocessing    = posterior_summary
  )
}
