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
#' person_covariates <- sim_res$hh_df %>%
#'   dplyr::select(hh_id, person_id, vacc_status) %>%   dplyr::distinct()
#'
#' df_for_stan <- sim_res$diagnostic_df %>%
#'   dplyr::left_join(person_covariates, by = c("hh_id", "person_id"))
#'
#' stan_input <- prepare_stan_data(
#'   df_clean = df_for_stan,
#'   use_vl_data      = 1,
#'   study_start_date = as.Date(study_start),
#'   study_end_date   = as.Date(study_end),
#'   surveillance_df  = surveillance_data,
#'   study_start_date = as.Date(study_start),
#'   study_end_date   = as.Date(study_end),
#'   imputation_params = VL_params_list)
#' fit <- fit_household_model(stan_input)
#' }
#' @export
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

#' Reconstruct Transmission Chains (Covariate-Aware)
#'
#' Identifies potential infectors for each infection event based on posterior
#' estimates, accounting for covariates if present.
#'
#' @param fit A stanfit object.
#' @param stan_data The list data passed to Stan.
#' @param min_prob_threshold Minimum probability threshold for links.
#'
#' @return A data frame of transmission links with probabilities.
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
#' person_covariates <- sim_res$hh_df %>%
#'   dplyr::select(hh_id, person_id, vacc_status) %>%   dplyr::distinct()
#'
#' df_for_stan <- sim_res$diagnostic_df %>%
#'   dplyr::left_join(person_covariates, by = c("hh_id", "person_id"))
#'
#' stan_input <- prepare_stan_data(
#'   df_clean = df_for_stan,
#'   use_vl_data      = 1,
#'   study_start_date = as.Date(study_start),
#'   study_end_date   = as.Date(study_end),
#'   surveillance_df  = surveillance_data,
#'   study_start_date = as.Date(study_start),
#'   study_end_date   = as.Date(study_end),
#'   imputation_params = VL_params_list)
#' fit <- fit_household_model(stan_input)
#' chains <- reconstruct_transmission_chains(fit, stan_input, min_prob_threshold = 0.001)
#' }
#' @export
reconstruct_transmission_chains <- function(fit, stan_data, min_prob_threshold = 0.01) {

  # 1. Extract Posterior Medians
  post <- rstan::extract(fit)
  p_beta1 <- stats::median(post$beta1)
  p_beta2 <- stats::median(post$beta2)
  p_alpha <- mean(post$alpha_comm)

  p_gen_shape <- stats::median(post$gen_shape)
  p_gen_rate  <- stats::median(post$gen_rate)
  p_ct50      <- stats::median(post$Ct50)
  p_slope     <- stats::median(post$slope_ct)

  p_phi   <- apply(post$phi_by_role, 2, stats::median)
  p_kappa <- apply(post$kappa_by_role, 2, stats::median)

  # Extract Covariate Coefficients
  p_beta_susc <- if (!is.null(post$beta_susc)) apply(as.matrix(post$beta_susc), 2, stats::median) else numeric(0)
  p_beta_inf  <- if (!is.null(post$beta_inf))  apply(as.matrix(post$beta_inf), 2, stats::median) else numeric(0)

  # 2. Setup Data
  N <- stan_data$N
  T_max <- stan_data$T
  I <- stan_data$I
  Y <- stan_data$Y
  V <- stan_data$V

  X_susc <- if (stan_data$K_susc > 0) stan_data$X_susc else matrix(0, nrow = N, ncol = 0)
  X_inf  <- if (stan_data$K_inf > 0)  stan_data$X_inf  else matrix(0, nrow = N, ncol = 0)

  infection_day <- rep(0, N)
  for (n in 1:N) {
    idx <- which(I[n, ] == 1)
    if (length(idx) > 0) infection_day[n] <- idx[1]
  }

  # Pre-calculate Gamma Curve
  g_curve <- numeric(T_max)
  for (d in 1:T_max) {
    g_curve[d] <- stats::dgamma(d, shape = p_gen_shape, rate = p_gen_rate)
  }
  g_curve <- g_curve / max(g_curve)

  # 3. Iterate through Infections
  results_list <- list()
  infected_episodes <- which(infection_day > 0)

  for (target_idx in infected_episodes) {

    t_inf <- infection_day[target_idx]
    hh    <- stan_data$hh_id[target_idx]
    role_t <- stan_data$role_id[target_idx]

    # Calculate Target Susceptibility (Phi * Covariates)
    log_susc_mod <- 0
    if (length(p_beta_susc) > 0) {
      log_susc_mod <- sum(X_susc[target_idx, ] * p_beta_susc)
    }
    phi_eff <- p_phi[role_t] * exp(log_susc_mod)

    # Community Hazard
    season_val <- stan_data$seasonal_forcing_mat[t_inf, role_t]
    lambda_comm <- phi_eff * p_alpha * season_val

    # Household Hazards
    hh_members <- which(stan_data$hh_id == hh)
    source_probs <- numeric(length(hh_members))
    source_ids   <- integer(length(hh_members))

    for (k in seq_along(hh_members)) {
      source_idx <- hh_members[k]
      source_ids[k] <- source_idx

      if (source_idx == target_idx) next
      if (stan_data$p_id[source_idx] == stan_data$p_id[target_idx]) next
      if (Y[source_idx, t_inf] == 0) next

      src_inf_day <- infection_day[source_idx]
      if (src_inf_day == 0 || src_inf_day > t_inf) next

      dt <- t_inf - src_inf_day + 1
      if (dt > T_max || dt < 1) next

      val_g <- g_curve[dt]
      val_v <- 0
      if (stan_data$use_vl_data == 1) {
        raw_v <- V[source_idx, t_inf]
        if (stan_data$vl_type == 1) {
          val_v <- (max(0, raw_v) / p_ct50)^p_slope
        } else {
          exponent <- -(raw_v - p_ct50) * p_slope
          val_v <- 1 / (1 + exp(-exponent))
        }
      }

      term_combined <- 0
      if (stan_data$use_vl_data == 0) {
        term_combined <- p_beta1 + (p_beta2 * val_g)
      } else {
        term_combined <- (p_beta1 * val_g) + (p_beta2 * val_v)
      }

      # Calculate Source Infectivity (Kappa * Covariates)
      log_inf_mod <- 0
      if (length(p_beta_inf) > 0) {
        log_inf_mod <- sum(X_inf[source_idx, ] * p_beta_inf)
      }
      kappa_eff <- p_kappa[stan_data$role_id[source_idx]] * exp(log_inf_mod)

      scaling <- (1 / max(stan_data$hh_size_people[hh], 1))^stan_data$delta
      h_source <- scaling * kappa_eff * term_combined

      lambda_source <- phi_eff * h_source
      source_probs[k] <- lambda_source
    }

    # 4. Save Probabilities
    total_hazard <- lambda_comm + sum(source_probs)

    if (total_hazard > 0) {
      # A. Community Link
      prob_comm <- lambda_comm / total_hazard
      if (prob_comm >= min_prob_threshold) {
        results_list[[length(results_list) + 1]] <- data.frame(
          target = target_idx, hh_id = hh, day = t_inf,
          source = "Community", prob = prob_comm, stringsAsFactors = FALSE
        )
      }

      # B. Household Links
      prob_hh_vec <- source_probs / total_hazard
      for (k in seq_along(prob_hh_vec)) {
        if (prob_hh_vec[k] >= min_prob_threshold) {
          results_list[[length(results_list) + 1]] <- data.frame(
            target = target_idx, hh_id = hh, day = t_inf,
            source = as.character(source_ids[k]),
            prob = prob_hh_vec[k], stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(results_list) > 0) do.call(rbind, results_list) else data.frame()
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
