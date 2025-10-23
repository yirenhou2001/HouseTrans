#' Simulate households, estimate transmission, and (optionally) summarize
#'
#' Generates synthetic household transmission data and runs the full parameter
#' estimation pipeline. By default, it prints a post-processing table comparing
#' the mean estimates to “true” values (bias and relative bias). Optionally,
#' it prints a per-individual data summary produced by the pipeline.
#'
#' @section Workflow:
#' \enumerate{
#'   \item Simulate households and test records.
#'   \item Summarize individuals and impute infection timelines.
#'   \item Build a person–day table and run repeated ML estimation.
#'   \item (Optional) Print data summary and post-processing table.
#' }
#'
#' @param synthetic_data Logical. Must be \code{TRUE}. If \code{FALSE}, an error is thrown with
#'   guidance to use \code{\link{TransmissionChainAnalysis}}.
#' @param data_summary Logical. If \code{TRUE}, prints \code{results$summarized_data}. Default \code{FALSE}.
#' @param postprocessing Logical. If \code{TRUE}, prints \code{\link{postprocessing_estimates}}. Default \code{TRUE}.
#'
#' @param n_households Integer. Number of households to simulate. Default \code{10}.
#' @param n_runs Integer. Repeated estimation runs. Default \code{10}.
#'
#' @param hh.size Integer. Household size (constant across HHs unless a random draw is passed).
#'   Default \code{sample(3:7, 1)}.
#' @param tests.per.week Integer. Tests per person per week (1–3). Default \code{1}.
#'
#' @param Covariates Logical. Generate synthetic covariates. Default \code{FALSE}.
#' @param Covariates_list Character vector of covariate names to generate when \code{Covariates=TRUE}.
#'   Default \code{c("Vaccination status", "Antibody Level")}.
#' @param Covariate_specifications List with per-covariate generation specs (type, dist, time_varying, params).
#'
#' @param day_series_covariates Logical. If \code{TRUE}, builds day-series list-columns for covariates in the person–day table. Default \code{TRUE}.
#' @param series_cols Character vector of covariate names (pre-normalization) to build as day-series; if \code{NULL}, all detected covariates are used.
#'
#' @param comm_covariate_cols Character vector of community-level covariates (no intercept) for the likelihood.
#' @param hh_covariate_cols Character vector of household covariates shared across roles.
#' @param hh_by_role Logical. If \code{TRUE}, fit role-specific household covariates.
#' @param hh_role_covariate_cols Named list with elements \code{infant}, \code{sibling}, \code{adult}, \code{elder}
#'   providing covariate names per role when \code{hh_by_role=TRUE}.
#' @param standardize_covariates Logical. Z-score non-binary numeric columns in model matrices. Default \code{TRUE}.
#'
#' @param lambda_comm,lambda_hh Numeric. L2 penalties for community and household covariate coefficients.
#'
#' @param p.comm.base.infant.fix,p.comm.multiplier.sibling,p.comm.multiplier.parent,p.comm.multiplier.elder
#'   Community infection baseline and role multipliers used in simulation.
#' @param p.hh.base.infant,p.hh.multiplier.sibling,p.hh.multiplier.parent,p.hh.multiplier.elder
#'   Household transmission baseline and role multipliers used in simulation.
#' @param p.imm.base.sibling,p.imm.base.parent,p.imm.base.elder Baseline immunity probabilities.
#' @param partial.immunity.infant,partial.immunity.sibling,partial.immunity.parent,partial.immunity.elder
#'   Breakthrough modifiers for partially immune individuals.
#' @param duration.latent,duration.infect.inf,multiplier.dur.sibpar Natural-history durations used in simulation.
#' @param p.detect Detection probability in testing algorithm.
#' @param amplitude,phase Seasonal modulation of community risk.
#' @param start_date,end_date Simulation window as \code{Date}.
#'
#' @param latent_par,report_par,infect_par Lists with Gamma(\code{shape}, \code{scale}) parameters for imputation of latent delay,
#'   reporting delay, and infectious period.
#'
#' @param start_par Numeric vector of initial parameters (auto-expanded if length mismatches model).
#' @param lambda Base L2 penalty for slope-like parameters (age multipliers and household role offsets).
#' @param lambda0 Penalty anchoring \code{delta0} near \code{delta0_true}.
#' @param lambda_alpha Penalty anchoring \code{alpha0} near \code{alpha0_true}.
#' @param delta0_true,alpha0_true Numeric anchors (logit scale) for community baseline and household baseline.
#'
#' @param true_values Optional named numeric vector of reference values. Names
#'   should match columns in \code{theta_mat}. Unmatched names are ignored; missing
#'   references yield \code{NA} bias/relative bias.
#'
#' @return (Invisibly) a list with elements:
#' \itemize{
#'   \item \code{results}: Output of \code{\link{main_parameter_estimation_pipeline}} (raw simulations, summaries, person-day, estimates).
#'   \item \code{postprocessing}: The table returned by \code{\link{postprocessing_estimates}} (or \code{NULL} if \code{postprocessing=FALSE}).
#' }
#'
#' @details This function is a thin wrapper around \code{\link{main_parameter_estimation_pipeline}}
#'   with \code{synthetic_data=TRUE}. Any extra columns generated (e.g., covariates) are automatically
#'   carried into the summarization and may be used in the likelihood via the covariate mapping arguments.
#'
#' @seealso \code{\link{TransmissionChainAnalysis}}, \code{\link{main_parameter_estimation_pipeline}},
#'   \code{\link{postprocessing_estimates}}
#'
#' @examples
#' \dontrun{
#' out <- GenSyn(n_households = 50, n_runs = 20, Covariates = TRUE)
#' names(out)
#' head(out$results$person_day)
#' }
#' @export
GenSyn <- function(
    synthetic_data = TRUE,
    data_summary   = FALSE,
    postprocessing = TRUE,
    n_households = 10,
    n_runs       = 10,
    hh.size        = sample(3:7, 1),
    tests.per.week = 1,
    Covariates               = FALSE,
    Covariates_list          = c("Vaccination status", "Antibody Level"),
    Covariate_specifications = NULL,
    day_series_covariates = TRUE,
    series_cols           = NULL,
    comm_covariate_cols   = NULL,
    hh_covariate_cols     = NULL,
    hh_by_role            = FALSE,
    hh_role_covariate_cols = NULL,
    standardize_covariates = TRUE,
    lambda_comm = 0.01,
    lambda_hh   = 0.01,
    p.comm.base.infant.fix   = 0.002,
    p.comm.multiplier.sibling = 1,
    p.comm.multiplier.parent  = 1,
    p.comm.multiplier.elder   = 1,
    p.hh.base.infant        = 0.2,
    p.hh.multiplier.sibling = 5.267686e-01,
    p.hh.multiplier.parent  = 8.008933e-01,
    p.hh.multiplier.elder   = 6.008933e-01,
    p.imm.base.sibling = 1e-10,
    p.imm.base.parent  = 1e-10,
    p.imm.base.elder   = 1e-10,
    partial.immunity.infant  = 1e-10,
    partial.immunity.sibling = 1e-10,
    partial.immunity.parent  = 1e-10,
    partial.immunity.elder   = 1e-10,
    duration.latent     = 1,
    duration.infect.inf = 2,
    multiplier.dur.sibpar = 0.5,
    p.detect            = 0.999,
    amplitude  = 2.65810 * 0,
    phase      = -0.408,
    start_date = as.Date("2024-09-21"),
    end_date   = as.Date("2025-04-17"),
    latent_par = list(shape = 2, scale = 1),
    report_par = list(shape = 1, scale = 1.5),
    infect_par = list(shape = 3, scale = 2),
    start_par    = c(-6, 0.02, -2, rep(0, 6)),
    lambda       = 0.01,
    lambda0      = 0.2,
    lambda_alpha = 5,
    delta0_true  = qlogis(0.002),
    alpha0_true  = qlogis(0.2),
    true_values = c(
      # Community infection
      delta0 = log(7.148217e-05),
      gamma2 = log(7.148217e-05 * 4.331956e+00) - log(7.148217e-05),
      gamma3 = log(7.148217e-05 * 1.835466e+00) - log(7.148217e-05),
      gamma4 = log(7.148217e-05 * 2) - log(7.148217e-05),
      # Household transmission
      alpha0 = log(0.2888953),
      beta2  = log(0.2888953 * 0.5267686) - log(0.2888953),
      beta3  = log(0.2888953 * 0.8008933) - log(0.2888953),
      beta4  = log(0.2888953 * 0.6008933) - log(0.2888953)
    )
) {
  if (isFALSE(synthetic_data)) {
    stop("Please use TransmissionChainAnalysis if not doing analysis on synthetic data!", call. = FALSE)
  }

  # Run the full pipeline for synthetic data
  results <- main_parameter_estimation_pipeline(
    user_data      = NULL,
    synthetic_data = TRUE,
    n_households   = n_households,
    n_runs         = n_runs,
    hh.size        = hh.size,
    tests.per.week = tests.per.week,
    Covariates               = Covariates,
    Covariates_list          = Covariates_list,
    Covariate_specifications = Covariate_specifications,
    day_series_covariates = day_series_covariates,
    series_cols           = series_cols,
    comm_covariate_cols   = comm_covariate_cols,
    hh_covariate_cols     = hh_covariate_cols,
    hh_by_role            = hh_by_role,
    hh_role_covariate_cols = hh_role_covariate_cols,
    standardize_covariates = standardize_covariates,
    lambda_comm = lambda_comm,
    lambda_hh   = lambda_hh,
    p.comm.base.infant.fix   = p.comm.base.infant.fix,
    p.comm.multiplier.sibling = p.comm.multiplier.sibling,
    p.comm.multiplier.parent  = p.comm.multiplier.parent,
    p.comm.multiplier.elder   = p.comm.multiplier.elder,
    p.hh.base.infant        = p.hh.base.infant,
    p.hh.multiplier.sibling = p.hh.multiplier.sibling,
    p.hh.multiplier.parent  = p.hh.multiplier.parent,
    p.hh.multiplier.elder   = p.hh.multiplier.elder,
    p.imm.base.sibling = p.imm.base.sibling,
    p.imm.base.parent  = p.imm.base.parent,
    p.imm.base.elder   = p.imm.base.elder,
    partial.immunity.infant  = partial.immunity.infant,
    partial.immunity.sibling = partial.immunity.sibling,
    partial.immunity.parent  = partial.immunity.parent,
    partial.immunity.elder   = partial.immunity.elder,
    duration.latent        = duration.latent,
    duration.infect.inf    = duration.infect.inf,
    multiplier.dur.sibpar  = multiplier.dur.sibpar,
    p.detect               = p.detect,
    amplitude              = amplitude,
    phase                  = phase,
    start_date             = start_date,
    end_date               = end_date,
    latent_par = latent_par,
    report_par = report_par,
    infect_par = infect_par,
    start_par    = start_par,
    lambda       = lambda,
    lambda0      = lambda0,
    lambda_alpha = lambda_alpha,
    delta0_true  = delta0_true,
    alpha0_true  = alpha0_true
  )

  # Prepare optional outputs without printing
  summarized_data <- if (isTRUE(data_summary)) results$summarized_data else NULL
  postprocessing_est <- if (isTRUE(postprocessing)) {
    postprocessing_estimates(results$estimates, true_values = true_values)
  } else {
    NULL
  }

  # Return visibly with a custom class so print() can control console output
  structure(
    list(
      results          = results,
      summarized_data  = summarized_data,   # convenience copy; also lives in results$summarized_data
      postprocessing   = postprocessing_est
    ),
    class = "GenSynResult"
  )
}

#' Print a GenSynResult
#'
#' Nicely prints sections available in a \code{GenSynResult} returned by
#' \code{\link{GenSyn}}. If present, the per-individual data summary and the
#' post-processing table are shown; otherwise concise guidance is printed.
#'
#' @param x A \code{GenSynResult} object, typically the result of \code{\link{GenSyn}}.
#' @param ... Passed to or from other methods (unused).
#'
#' @details
#' This S3 method is invoked when a \code{GenSynResult} is printed (e.g., typing
#' the object at the console or calling \code{print(x)}). The method does not
#' perform any computation and does not modify global state; it only formats and
#' prints components that were included in the returned object.
#'
#' @return \code{x}, returned invisibly.
#'
#' @seealso \code{\link{GenSyn}}
#' @exportS3Method print GenSynResult
print.GenSynResult <- function(x, ...) {
  cat("GenSyn result\n")
  if (!is.null(x$summarized_data)) {
    cat("\n--- Data summary ---\n")
    print(x$summarized_data)
  } else {
    cat("\n(Data summary not requested; set data_summary = TRUE to include.)\n")
  }
  if (!is.null(x$postprocessing)) {
    cat("\n--- Post-processing of estimates ---\n")
    print(x$postprocessing)
  } else {
    cat("\n(Post-processing not requested; set postprocessing = TRUE to include.)\n")
  }
  invisible(x)
}
