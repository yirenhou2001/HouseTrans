#' Analyze User-Supplied Household Transmission Data
#'
#' Runs Bayesian estimation on user-provided household testing or episode data.
#' This function does not simulate data; it expects input in tabular form.
#'
#' @param user_data A data.frame or list of data.frames with household data.
#'   Accepts two formats:
#'   \itemize{
#'     \item \strong{Long testing table}: columns \code{HH}, \code{individual_ID},
#'           \code{role}, \code{test_date}, \code{infection_status} (and optionally
#'           a viral load column)
#'     \item \strong{Per-person episode table}: columns \code{hh_id}, \code{person_id},
#'           \code{role}, \code{infection_time}, \code{infectious_start},
#'           \code{infectious_end}, \code{infection_resolved}
#'   }
#' @param start_date,end_date Date objects or character; study period bounds.
#' @param surveillance_df Optional data frame with 'date' and 'cases' columns for
#'   seasonal forcing.
#' @param seasonal_forcing_list Optional named list with seasonal forcing vectors
#'   for each role.
#' @param max_days Integer; maximum time horizon (days).
#' @param use_vl_data Logical; whether to use viral load data in estimation.
#' @param priors Named list; prior specifications for Stan model. Each element
#'   should be a list with \code{dist} ("normal", "uniform", "lognormal") and
#'   \code{params} (parameter vector).
#' @param covariates_susceptibility Character vector; column names for susceptibility
#'   covariates in the data.
#' @param covariates_infectivity Character vector; column names for infectivity
#'   covariates in the data.
#' @param recovery_params Named list; Gamma parameters (shape, scale) for immunity
#'   tail duration by role.
#' @param imputation_params Named list; parameters for viral curve imputation by role.
#' @param role_levels Character vector; canonical role levels for normalization.
#' @param stan_chains,stan_iter,stan_warmup Integers; Stan sampling controls.
#' @param stan_control List; Stan control parameters.
#' @param stan_cores Integer; number of CPU cores for Stan.
#' @param stan_refresh Integer; refresh rate for Stan output.
#' @param seed Integer; random seed for reproducibility.
#'
#' @return An object of class \code{"TransmissionChainResult"} containing:
#' \describe{
#'   \item{$call}{The matched call}
#'   \item{$user_data}{The processed user data}
#'   \item{$stan_data}{Data prepared for Stan}
#'   \item{$fit}{The stanfit object}
#'   \item{$postprocessing}{Tidy posterior summary}
#'   \item{$transmission_chains}{Reconstructed transmission links}
#' }
#' @examples
#' \dontrun{
#' # Per-person episode format
#' T_max <- 30
#' df_person <- data.frame(
#'   hh_id = c("HH1","HH1","HH1","HH2","HH2","HH2"),
#'   person_id = c(1, 2, 3, 1, 2, 3),
#'   role = c("adult","infant","elderly","adult","infant","elderly"),
#'   infection_time = c(2, 4, NA, 1, 3, NA),
#'   infectious_start = c(3, 6, NA, 2, 5, NA),
#'   infectious_end = c(8, 9, NA, 7, 9, NA),
#'   infection_resolved = c(9, 10, NA, 8, 10, NA)
#' )
#'
#' result <- TransmissionChainAnalysis(
#'   user_data = df_person,
#'   max_days = T_max,
#'   stan_chains = 2,
#'   stan_iter = 1000,
#'   stan_warmup = 500
#' )
#'
#' print(result)
#' plot(result, which = "posterior")
#' }
#' @seealso \code{\link{GenSyn}}, \code{\link{plot.TransmissionChainResult}}
#' @export
TransmissionChainAnalysis <- function(
    user_data,
    start_date = as.Date("2024-01-01"),
    end_date = as.Date("2024-12-31"),
    surveillance_df = NULL,
    seasonal_forcing_list = NULL,
    max_days = NULL,
    use_vl_data = TRUE,
    priors = NULL,
    covariates_susceptibility = NULL,
    covariates_infectivity = NULL,
    recovery_params = NULL,
    imputation_params = NULL,
    role_levels = c("adult", "infant", "toddler", "elderly"),
    stan_chains = 4,
    stan_iter = 2000,
    stan_warmup = 1000,
    stan_control = list(adapt_delta = 0.95, max_treedepth = 15),
    stan_cores = 4,
    stan_refresh = 50,
    seed = NULL
) {

  call <- match.call()

  # Validate input
  if (is.null(user_data)) {
    stop("`user_data` must be supplied to TransmissionChainAnalysis().", call. = FALSE)
  }

  # Handle list of data.frames
  if (is.list(user_data) && !is.data.frame(user_data)) {
    ok <- all(vapply(user_data, is.data.frame, logical(1)))
    if (!ok) {
      stop("If `user_data` is a list, it must contain only data.frames.", call. = FALSE)
    }
    user_data <- data.table::rbindlist(user_data, use.names = TRUE, fill = TRUE)
  }

  if (!is.data.frame(user_data)) {
    stop("`user_data` must be a data.frame (or a list of data.frames).", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  if (is.null(priors)) priors <- default_priors

  # Normalize roles
  if ("role" %in% names(user_data)) {
    user_data$role <- .norm_role(user_data$role)
  }

  # Check data format and prepare
  req_person_core <- c("hh_id", "person_id", "role", "infection_time",
                       "infectious_start", "infectious_end", "infection_resolved")
  req_long <- c("HH", "individual_ID", "role", "test_date", "infection_status")

  is_person_format <- all(req_person_core %in% names(user_data))
  is_long_format <- all(req_long %in% names(user_data))

  if (!is_person_format && !is_long_format) {
    stop(
      "TransmissionChainAnalysis could not interpret `user_data`.\n\n",
      "Please ensure that your data are either:\n",
      "  (a) a long household-testing table with columns: HH, individual_ID,\n",
      "      role, test_date, infection_status;\n",
      "  or\n",
      "  (b) a per-person table with columns: hh_id, person_id, role,\n",
      "      infection_time, infectious_start, infectious_end, infection_resolved.\n",
      call. = FALSE
    )
  }

  # Calculate max_days from dates if not provided
  if (is.null(max_days)) {
    max_days <- as.integer(as.Date(end_date) - as.Date(start_date)) + 1
  }

  # Prepare Stan data
  if (is_long_format) {
    # Long format - use full prepare_stan_data
    stan_data <- prepare_stan_data(
      df_clean = user_data,
      surveillance_df = surveillance_df,
      role_levels = role_levels,
      study_start_date = as.Date(start_date),
      study_end_date = as.Date(end_date),
      seasonal_forcing_list = seasonal_forcing_list,
      use_vl_data = use_vl_data,
      covariates_susceptibility = covariates_susceptibility,
      covariates_infectivity = covariates_infectivity,
      recovery_params = recovery_params,
      imputation_params = imputation_params,
      priors = priors,
      seed = seed
    )
  } else {
    # Per-person format - use simpler path
    # Add episode_id if missing
    if (!"episode_id" %in% names(user_data)) {
      user_data$episode_id <- ifelse(is.na(user_data$infection_time), NA_integer_, 1L)
    }

    # Split into households
    households <- split(user_data, user_data$hh_id, drop = TRUE)

    stan_data <- build_stan_household_arrays(
      households = households,
      T_max = max_days,
      seasonal_forcing_list = seasonal_forcing_list,
      priors = priors,
      use_vl_data = use_vl_data
    )
  }

  # Fit model
  fit <- fit_household_model(
    stan_data = stan_data,
    iter = stan_iter,
    chains = stan_chains,
    warmup = stan_warmup,
    control = stan_control,
    cores = stan_cores,
    refresh = stan_refresh
  )

  # Postprocess
  posterior_summary <- postprocess_stan_fit(fit)

  # Reconstruct transmission chains
  chains <- tryCatch(
    reconstruct_transmission_chains(fit, stan_data, min_prob_threshold = 0.01),
    error = function(e) NULL
  )

  out <- list(
    call = call,
    user_data = user_data,
    surveillance_df = surveillance_df,
    start_date = start_date,
    end_date = end_date,
    stan_data = stan_data,
    fit = fit,
    postprocessing = posterior_summary,
    transmission_chains = chains
  )

  class(out) <- "TransmissionChainResult"
  out
}


#' Print method for TransmissionChainResult
#'
#' @param x A \code{TransmissionChainResult} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' \dontrun{
#' # Per-person episode format
#' T_max <- 30
#' df_person <- data.frame(
#'   hh_id = c("HH1","HH1","HH1","HH2","HH2","HH2"),
#'   person_id = c(1, 2, 3, 1, 2, 3),
#'   role = c("adult","infant","elderly","adult","infant","elderly"),
#'   infection_time = c(2, 4, NA, 1, 3, NA),
#'   infectious_start = c(3, 6, NA, 2, 5, NA),
#'   infectious_end = c(8, 9, NA, 7, 9, NA),
#'   infection_resolved = c(9, 10, NA, 8, 10, NA)
#' )
#'
#' result <- TransmissionChainAnalysis(
#'   user_data = df_person,
#'   max_days = T_max,
#'   stan_chains = 2,
#'   stan_iter = 1000,
#'   stan_warmup = 500
#' )
#'
#' print(result)
#' }
#' @export
#' @method print TransmissionChainResult
print.TransmissionChainResult <- function(x, ...) {
  cat("TransmissionChainAnalysis Result\n")
  cat("================================\n\n")

  if (!is.null(x$user_data)) {
    n_hh <- length(unique(x$user_data$hh_id %||% x$user_data$HH))
    n_people <- nrow(x$user_data)
    cat("Households:", n_hh, "\n")
    cat("Person-episodes:", n_people, "\n\n")
  }

  if (!is.null(x$postprocessing) && nrow(x$postprocessing) > 0) {
    cat("--- Posterior Summary (Key Parameters) ---\n")
    key_params <- c("beta1", "beta2", "alpha_comm", "phi_by_role", "kappa_by_role")
    subset <- x$postprocessing[grepl(paste(key_params, collapse = "|"), x$postprocessing$Parameter), ]
    if (nrow(subset) > 0) {
      print(head(subset, 20))
    } else {
      print(head(x$postprocessing, 10))
    }
    cat("\n")
  }

  cat("Access components: $user_data, $fit, $postprocessing, $transmission_chains\n")
  cat("Use plot(result, which = '...') for visualizations.\n")

  invisible(x)
}


#' Plot method for TransmissionChainResult
#'
#' @param x A \code{TransmissionChainResult} object.
#' @param which Character vector specifying which plots to generate. Options:
#'   \code{"posterior"}, \code{"transmission_chains"}, \code{"covariate_effects"},
#'   \code{"epidemic_curve"}, or \code{"all"}.
#' @param print Logical; whether to print plots immediately.
#' @param hh_id Integer; household ID for transmission chain plot (default: 1).
#' @param prob_cutoff Numeric; minimum probability threshold for chain links (default: 0).
#' @param bin_width Integer; number of days per bin for epidemic curve (default: 7).
#' @param ... Additional arguments (unused).
#' @return A named list of ggplot objects (invisibly if \code{print = TRUE}).
#' @examples
#' \dontrun{
#' # Per-person episode format
#' T_max <- 30
#' df_person <- data.frame(
#'   hh_id = c("HH1","HH1","HH1","HH2","HH2","HH2"),
#'   person_id = c(1, 2, 3, 1, 2, 3),
#'   role = c("adult","infant","elderly","adult","infant","elderly"),
#'   infection_time = c(2, 4, NA, 1, 3, NA),
#'   infectious_start = c(3, 6, NA, 2, 5, NA),
#'   infectious_end = c(8, 9, NA, 7, 9, NA),
#'   infection_resolved = c(9, 10, NA, 8, 10, NA)
#' )
#'
#' result <- TransmissionChainAnalysis(
#'   user_data = df_person,
#'   max_days = T_max,
#'   stan_chains = 2,
#'   stan_iter = 1000,
#'   stan_warmup = 500
#' )
#'
#' plot(result, which = "posterior")
#' }
#' @export
#' @method plot TransmissionChainResult
plot.TransmissionChainResult <- function(x, which = "posterior", print = TRUE,
                                         hh_id = 1, prob_cutoff = 0,
                                         bin_width = 7, ...) {

  if ("all" %in% which) {
    which <- c("posterior", "transmission_chains", "covariate_effects", "epidemic_curve")
  }

  plots <- list()

  if ("posterior" %in% which && !is.null(x$fit)) {
    plots$posterior <- tryCatch(
      plot_posterior_distributions(x$fit),
      error = function(e) .empty_plot(e$message)
    )
  }

  if ("transmission_chains" %in% which && !is.null(x$transmission_chains) && !is.null(x$stan_data)) {
    plots$transmission_chains <- tryCatch(
      plot_household_timeline(x$transmission_chains, x$stan_data,
                              target_hh_id = hh_id,
                              start_date_str = as.character(x$start_date),
                              prob_cutoff = prob_cutoff),
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

  if ("epidemic_curve" %in% which && !is.null(x$user_data)) {
    plots$epidemic_curve <- tryCatch(
      plot_user_epidemic_curve(x$user_data, x$start_date, x$surveillance_df, bin_width),
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
