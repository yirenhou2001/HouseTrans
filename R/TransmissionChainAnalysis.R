#' TransmissionChainAnalysis: analyze user-supplied household data
#'
#' Runs the end-to-end workflow on \emph{user-provided} data using
#' Bayesian Stan estimation. This function does not simulate data;
#' it expects input in tabular form (see Details).
#'
#' @param user_data A data.frame or a list of data.frames describing household testing
#'   or episode data (see Details for accepted shapes). Lists are row-bound internally.
#' @param plots Character vector of plot names to build when compatible data are present:
#'   \code{c("daily","weekly","timeline","sar")} or \code{"all"}.
#' @param print_plots Logical; print plots if produced.
#' @param index_vl_column Optional character; name of the viral-load column used by SAR-by-VL plotting. Defaults to \code{"vl_test"} when present.
#' @param start_date,end_date \code{Date} study window.
#' @param seasonal_forcing_list Optional named list of role vectors for forcing.
#' @param max_days Integer; maximum modeled days.
#' @param stan_file Path to a Stan model file.
#' @param stan_chains,stan_iter,stan_warmup Integers; Stan sampling controls.
#' @param stan_control List; passed to \code{rstan::sampling(..., control = ...)}.
#' @param stan_init Character or function; Stan initialization.
#' @param stan_refresh Integer; Stan refresh rate.
#' @param stan_cores Integer; CPU cores for Stan.
#' @param vl_mode One of \code{c("auto","from_long")} indicating how to derive VL.
#' @param vl_source One of \code{c("column","simulate","none")} indicating the VL source.
#' @param vl_column Optional character; the VL column when \code{vl_source="column"}.
#' @param role_levels Character vector; canonical role levels for normalization.
#'
#' @details
#' \strong{Accepted input shapes}:
#' \itemize{
#'   \item Long testing table with columns like \code{HH}, \code{individual_ID}, \code{role},
#'         \code{test_date}, \code{infection_status}, and optionally a VL column.
#'   \item Per-person episode table with \code{hh_id}, \code{person_id}, \code{role},
#'         \code{infection_time}, \code{infectious_start}, \code{infectious_end},
#'         and optionally a viral-load trajectory.
#' }
#'
#' @return An object of class \code{"TransmissionChainResult"} containing:
#' \itemize{
#'   \item \code{$results}: Contains \code{fit} (the Stan object) and \code{posterior_summary}.
#'   \item \code{$postprocessing}: Stan posterior summary.
#'   \item \code{$plot_list}: named list of ggplot objects when built.
#' }
#'
#' @seealso \code{\link{main_parameter_estimation_pipeline}}, \code{\link{GenSyn}}
#'
#' @examples
#' \dontrun{
#' T_max <- 12
#' df_person <- data.frame(
#'   hh_id = c("HH1","HH1","HH1","HH2","HH2","HH2"),
#'   person_id = c(1,2,3,1,2,3),
#'   role = c("adult","child","elderly","adult","child","elderly"),
#'   infection_time = c(2, 4, NA, 1, 3, NA),
#'   infectious_start = c(3, 6, NA, 2, 5, NA),
#'   infectious_end = c(8, 9, NA, 7, 9, NA),
#'   infection_resolved = c(9, 10, NA, 8, 10, NA)
#' )
#' seasonal_forcing_list <- list(
#'   adult = rep(1, T_max), child = rep(1, T_max),
#'   elderly = rep(1, T_max), toddler = rep(1, T_max)
#' )
#' result <- TransmissionChainAnalysis(
#'   user_data = df_person,
#'   seasonal_forcing_list = seasonal_forcing_list,
#'   max_days = T_max,
#'   stan_chains = 1, stan_iter = 300, stan_warmup = 150
#' )
#' }
#' @export
TransmissionChainAnalysis <- function(
    user_data,

    plots          = c("daily","weekly","timeline","sar"),
    print_plots    = FALSE,
    index_vl_column = "vl_test",

    # Analysis window
    start_date = as.Date("2024-01-01"),
    end_date   = as.Date("2024-12-31"),

    # Seasonal forcing
    seasonal_forcing_list = NULL,
    max_days = 365,

    # Stan controls
    stan_file    = "model.stan",
    stan_chains = 4, stan_iter = 2000, stan_warmup = 1000,
    stan_control = list(adapt_delta = 0.99, max_treedepth = 20),
    stan_init = "random", stan_refresh = 50, stan_cores = 4,

    # VL options
    vl_mode   = c("auto","from_long"),
    vl_source = c("column","simulate","none"),
    vl_column = NULL,
    role_levels = c("adult","child","elderly","toddler")
) {
  if (is.null(user_data)) {
    stop("`user_data` must be supplied to TransmissionChainAnalysis().", call. = FALSE)
  }

  # Allow list of data.frames but not other list types
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

  vl_mode   <- match.arg(vl_mode)
  vl_source <- match.arg(vl_source)

  # Normalize roles
  if ("role" %in% names(user_data)) {
    user_data$role <- .norm_role(user_data$role)
  }

  households <- tryCatch(
    prepare_stan_households_from_user_data(
      user_data   = user_data,
      role_levels = role_levels,
      vl_mode     = vl_mode,
      vl_source   = vl_source,
      vl_column   = vl_column,
      start_date  = start_date,
      end_date    = end_date
    ),
    error = function(e) {
      msg <- paste0(
        "TransmissionChainAnalysis could not interpret `user_data`.\n\n",
        "Please ensure that your data are either:\n",
        "  (a) a long household-testing table with columns such as HH, individual_ID,\n",
        "      role, test_date, infection_status (and optionally a viral-load column);\n",
        "  or\n",
        "  (b) a per-person table with columns such as hh_id, person_id, role,\n",
        "      infection_time, infectious_start, infectious_end, and optionally\n",
        "      vl_full_trajectory.\n\n",
        "Original error from prepare_stan_households_from_user_data():\n  ",
        conditionMessage(e)
      )
      stop(msg, call. = FALSE)
    }
  )

  stan_data <- build_stan_household_arrays(
    households = households,
    T_max = max_days,
    seasonal_forcing_list = seasonal_forcing_list
  )

  fit <- run_household_stan(
    stan_data,
    stan_file = stan_file,
    chains    = stan_chains,
    iter      = stan_iter,
    warmup    = stan_warmup,
    control   = stan_control,
    init      = stan_init,
    refresh   = stan_refresh,
    cores     = stan_cores
  )

  posterior_summary <- postprocess_stan_fit(fit)
  postprocessing_tbl <- posterior_summary

  out <- list(
    call              = match.call(),
    results = list(
      raw_simulation    = households,
      stan_data         = stan_data,
      fit               = fit,
      posterior_summary = posterior_summary
    ),
    postprocessing = postprocessing_tbl,
    plots          = plots,
    plot_list      = list()
  )
  class(out) <- "TransmissionChainResult"
  return(out)
}


#' Print a TransmissionChainResult
#'
#' Nicely prints sections available in a \code{TransmissionChainResult} returned by
#' \code{\link{TransmissionChainAnalysis}}.
#'
#' @param x A \code{TransmissionChainResult} object.
#' @param ... Passed to or from other methods (unused).
#'
#' @return \code{x}, returned invisibly.
#'
#' @seealso \code{\link{TransmissionChainAnalysis}}
#' @exportS3Method print TransmissionChainResult
print.TransmissionChainResult <- function(x, ...) {
  cat("TransmissionChainAnalysis result\n\n")

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
