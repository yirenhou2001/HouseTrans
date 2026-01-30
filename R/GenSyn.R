#' GenSyn: household simulation & estimation wrapper
#'
#' Runs \code{\link{main_parameter_estimation_pipeline}()} to simulate
#' household data, fit via Bayesian Stan, and assemble a \code{"GenSynResult"}.
#'
#' @param n_households Integer; number of households to simulate.
#' @param print_plots Logical; print plots if produced.
#' @param plots Character vector of plot names (\code{"daily"}, \code{"weekly"}, \code{"sar"}, \code{"timeline"}) or \code{"all"}.
#' @param index_vl_column Character; viral-load column name for plotting (default \code{"vl_test"}).
#' @param start_date,end_date \code{Date}; study window.
#' @param seasonal_forcing_list Optional named list of role vectors; seasonal forcing.
#' @param max_days Integer; maximum simulated days.
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
#' @param detect_threshold_log10 Numeric; VL detection threshold in log10.
#' @param detect_threshold_Ct Numeric; Ct detection threshold.
#' @param surveillance_interval Integer; days between scheduled tests.
#' @param test_daily Logical; switch to daily testing after first detection.
#' @param viral_testing Character; \code{"viral load"} or \code{"Ct"}.
#' @param V_ref,V_rho Numeric; viral load reference and power.
#' @param Ct_50,Ct_delta Numeric; Ct-based infectivity parameters.
#' @param VL_params_list Named list; role-specific VL trajectory parameters.
#' @param Ct_params_list Named list; role-specific Ct trajectory parameters.
#' @param household_profile_list Named list; household composition probabilities.
#' @param stan_file Path to a Stan model file.
#' @param stan_chains,stan_iter,stan_warmup Integers; Stan sampling controls.
#' @param stan_control List; Stan \code{control} list.
#' @param stan_init Character or function; Stan initialization.
#' @param stan_refresh Integer; Stan refresh rate.
#' @param stan_cores Integer; CPU cores for Stan.
#'
#' @return A \code{"GenSynResult"} list with elements: \code{$call},
#'   \code{$n_households}, \code{$results}, \code{$postprocessing}, \code{$plot_list}.
#'
#' @seealso \code{\link{main_parameter_estimation_pipeline}}, \code{\link{TransmissionChainAnalysis}}
#' @examples
#' \dontrun{
#' # Simulate and estimate
#' seasonal_forcing_list <- list(
#'   adult=rep(0.1,365), child=rep(0.1,365), elderly=rep(0.1,365), toddler=rep(0.1,365)
#' )
#' fit <- GenSyn(
#'   n_households=50,
#'   seasonal_forcing_list=seasonal_forcing_list, max_days=365,
#'   stan_chains=4, stan_iter=2000, stan_warmup=1000, stan_cores=4
#' )
#' }
#' @export
GenSyn <- function(
    n_households       = 50,
    print_plots        = TRUE,
    plots              = c("daily","weekly","timeline","sar"),
    index_vl_column    = "vl_test",

    # Study window
    start_date = as.Date("2024-01-01"),
    end_date   = as.Date("2024-12-31"),

    # Seasonal forcing
    seasonal_forcing_list = NULL,
    max_days              = 365,

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
    stan_file    = "model.stan",
    stan_chains  = 4,
    stan_iter    = 2000,
    stan_warmup  = 1000,
    stan_control = list(adapt_delta = 0.99, max_treedepth = 20),
    stan_init    = "random",
    stan_refresh = 50,
    stan_cores   = 4
) {

  results <- main_parameter_estimation_pipeline(
    user_data          = NULL,
    synthetic_data     = TRUE,
    n_households       = n_households,
    start_date         = start_date,
    end_date           = end_date,
    seasonal_forcing_list = seasonal_forcing_list,
    max_days           = max_days,
    alpha_comm_by_role = alpha_comm_by_role,
    beta1              = beta1,
    beta2              = beta2,
    delta              = delta,
    phi_by_role        = phi_by_role,
    kappa_by_role      = kappa_by_role,
    latent_shape       = latent_shape,
    latent_scale       = latent_scale,
    infectious_shape   = infectious_shape,
    infectious_scale   = infectious_scale,
    resolve_shape      = resolve_shape,
    resolve_scale      = resolve_scale,
    peak_day           = peak_day,
    width              = width,
    ptrans_threshold   = ptrans_threshold,
    detect_threshold_log10 = detect_threshold_log10,
    detect_threshold_Ct = detect_threshold_Ct,
    surveillance_interval = surveillance_interval,
    test_daily         = test_daily,
    viral_testing      = viral_testing,
    V_ref              = V_ref,
    V_rho              = V_rho,
    Ct_50              = Ct_50,
    Ct_delta           = Ct_delta,
    VL_params_list     = VL_params_list,
    Ct_params_list     = Ct_params_list,
    household_profile_list = household_profile_list,
    stan_file          = stan_file,
    stan_chains        = stan_chains,
    stan_iter          = stan_iter,
    stan_warmup        = stan_warmup,
    stan_control       = stan_control,
    stan_init          = stan_init,
    stan_refresh       = stan_refresh,
    stan_cores         = stan_cores
  )

  postprocessing <- results$posterior_summary

  plot_list <- list()
  if (length(plots)) {
    plot_list <- .make_transmission_plots(
      results,
      which = plots,
      print = print_plots,
      index_vl_column = index_vl_column
    )
  }

  out <- list(
    call              = match.call(),
    n_households      = n_households,
    plots             = plots,
    plot_list         = plot_list,
    results           = results,
    postprocessing    = postprocessing
  )
  class(out) <- "GenSynResult"
  out
}


#' Print method for \code{GenSynResult}
#'
#' @param x A \code{GenSynResult} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @exportS3Method print GenSynResult
print.GenSynResult <- function(x, ...) {
  cat("GenSyn result\n\n")
  cat("Households: ", x$n_households, "\n\n", sep = "")

  if (!is.null(x$postprocessing)) {
    cat("--- Posterior summary ---\n")
    print(x$postprocessing)
  } else {
    cat("(No posterior summary stored in `x$postprocessing`.)\n")
  }

  cat("\n(Full posterior summary stored in `x$postprocessing`)\n")
  cat("(Stan fit object stored in `x$results$fit`)\n")
  invisible(x)
}
