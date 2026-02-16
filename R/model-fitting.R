#' @title Model Fitting Functions
#' @description Functions for fitting Stan models to household data
#' @name model_fitting
#' @keywords internal
NULL


#' Fit Household Transmission Model
#'
#' Fits the compiled Stan model to the prepared household data.
#'
#' @param stan_data A list of data formatted by \code{prepare_stan_data}.
#' @param iter Integer; number of iterations per chain (including warmup).
#' @param chains Integer; number of Markov chains.
#' @param warmup Integer; number of warmup iterations.
#' @param init_fun Function or List; initial values for the sampler.
#' @param control List; Stan control parameters.
#' @param cores Integer; number of CPU cores.
#' @param refresh Integer; refresh rate for output.
#' @param ... Additional arguments passed to \code{rstan::sampling}.
#'
#' @return A \code{stanfit} object containing the posterior samples.
#' @keywords internal
fit_household_model <- function(stan_data,
                                iter = 2000,
                                chains = 4,
                                warmup = 1000,
                                init_fun = NULL,
                                control = list(adapt_delta = 0.95, max_treedepth = 15),
                                cores = 4,
                                refresh = 50,
                                ...) {

  # 1. Define Default Initial Values
  # These defaults are crucial for convergence in transmission models
  if (is.null(init_fun)) {
    init_fun <- function() {
      list(
        log_beta1 = -5.3,  # Approx log(0.005)
        log_beta2 = -5.3,

        # Initialize multipliers near 1.0 (log scale 0) to start neutral
        log_phi_by_role_raw = rep(0.1, max(1, stan_data$R - 1)),
        log_kappa_by_role_raw = rep(0.1, max(1, stan_data$R - 1)),

        # Viral Load parameters (if used)
        V_ref = 3.0,
        V_rho = 2.5,

        # Community infection rate placeholders
        log_beta3 = 0,
        log_beta4 = 0
      )
    }
  }

  # 2. Run Sampler
  # Use pre-compiled Stan model from the package
  out <- rstan::sampling(
    object = stanmodels$household_transmission,
    data = stan_data,
    iter = iter,
    chains = chains,
    warmup = warmup,
    init = init_fun,
    control = control,
    cores = cores,
    refresh = refresh,
    ...
  )

  return(out)
}


#' Postprocess Stan Fit
#'
#' Extracts and summarizes posterior distributions from a stanfit object.
#'
#' @param fit A \code{stanfit} object.
#'
#' @return A data frame with parameter summaries.
#' @keywords internal
postprocess_stan_fit <- function(fit) {
  if (is.null(fit)) return(data.frame())

  tryCatch({
    summ <- rstan::summary(fit)$summary
    df <- as.data.frame(summ)
    df$Parameter <- rownames(summ)
    rownames(df) <- NULL

    # Reorder columns
    cols <- c("Parameter", "mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")
    cols <- cols[cols %in% names(df)]
    df <- df[, cols]

    return(df)
  }, error = function(e) {
    warning("Failed to postprocess Stan fit: ", e$message)
    return(data.frame())
  })
}
