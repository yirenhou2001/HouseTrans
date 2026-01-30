#' Prepare per-household inputs for the Stan RSV/VL model
#'
#' Converts user data to a list of per-household data frames required by the
#' Stan RSV/VL pipeline. Accepts either a long test-day table or a per-person
#' episodes table.
#'
#' @param user_data Either:
#'   \itemize{
#'     \item \emph{Long format:} columns \code{HH}, \code{individual_ID},
#'           \code{role}, \code{test_date}, \code{infection_status}
#'           (optionally a VL column); or
#'     \item \emph{Per-person format:} columns \code{hh_id}, \code{person_id},
#'           \code{role}, \code{infection_time}, \code{infectious_start},
#'           \code{infectious_end} (optionally \code{vl_full_trajectory}).
#'   }
#' @param role_levels Character vector of allowed roles after normalization.
#' @param vl_mode Character; currently used only to trigger trajectory building.
#' @param vl_source Character; one of \code{"none"}, \code{"column"}, \code{"simulate"}.
#' @param vl_column Optional name of the VL column when \code{vl_source = "column"}.
#' @param start_date,end_date Optional \code{Date} bounds.
#'
#' @return A named list of data frames (one per household).
prepare_stan_households_from_user_data <- function(
    user_data,
    role_levels = c("adult", "child", "elderly", "toddler"),
    vl_mode   = c("from_long", "auto"),
    vl_source = c("none", "column", "simulate"),
    vl_column = NULL,
    start_date = NULL,
    end_date   = NULL
) {
  vl_mode   <- match.arg(vl_mode)
  vl_source <- match.arg(vl_source)

  norm_role_vec <- function(x) {
    x <- tolower(as.character(x))
    x[x == "infant"]  <- "toddler"; x[x == "baby"] <- "toddler"
    x[x == "sibling"] <- "child"; x[x == "kid"] <- "child"
    x[x == "parent"]  <- "adult"; x[x == "mother"] <- "adult"; x[x == "father"] <- "adult"
    x[x == "elder"]   <- "elderly"; x[x == "grandparent"] <- "elderly"
    x
  }

  req_person_core <- c("hh_id", "person_id", "role", "infection_time", "infectious_start", "infectious_end", "infection_resolved")
  if (inherits(user_data, "data.frame") && all(req_person_core %in% names(user_data))) {
    df <- data.frame(user_data, stringsAsFactors = FALSE)
    df$role <- norm_role_vec(df$role)
    bad_roles <- setdiff(unique(df$role), role_levels)
    if (length(bad_roles)) stop("Unknown role labels: ", paste(bad_roles, collapse = ", "), call. = FALSE)
    if (!"detection_time" %in% names(df)) {
      df$detection_time <- ifelse(is.finite(df$infectious_start), df$infectious_start, df$infection_time)
    }
    if (!"vl_full_trajectory" %in% names(df)) {
      vl_list <- vector("list", nrow(df))
      for (i in seq_len(nrow(df))) {
        it <- df$infection_time[i]; ir <- df$infection_resolved[i]
        vl_list[[i]] <- if (is.na(it) || is.na(ir) || ir < it) numeric(0) else rep(-5, ir - it + 1L)
      }
      df$vl_full_trajectory <- vl_list
    }
    df$hh_id <- as.character(df$hh_id); df$person_id <- as.integer(df$person_id)
    df <- df[order(df$hh_id, df$person_id), , drop = FALSE]
    return(split(df, df$hh_id, drop = TRUE))
  }

  required_long <- c("HH", "individual_ID", "role", "test_date", "infection_status")
  if (!all(required_long %in% names(user_data))) {
    stop("Long-format user_data missing columns: ", paste(setdiff(required_long, names(user_data)), collapse = ", "), call. = FALSE)
  }

  dt <- data.table::as.data.table(user_data)
  dt[, HH := as.character(HH)]; dt[, individual_ID := as.integer(individual_ID)]
  dt[, role := norm_role_vec(role)]; dt[, infection_status := as.integer(infection_status)]

  if (inherits(dt$test_date, "Date")) {
    if (is.null(start_date)) start_date <- min(dt$test_date, na.rm = TRUE)
    dt[, day := as.integer(test_date - start_date) + 1L]
  } else { dt[, day := as.integer(test_date)] }

  bad_roles <- setdiff(unique(dt$role), role_levels)
  if (length(bad_roles)) stop("Unknown role labels: ", paste(bad_roles, collapse = ", "), call. = FALSE)

  person <- dt[, {
    infected_days <- day[infection_status == 1L]
    if (length(infected_days) == 0L) {
      list(infection_time = NA_integer_, infectious_start = NA_integer_, infectious_end = NA_integer_,
           detection_time = NA_integer_, infection_resolved = NA_integer_)
    } else {
      it <- min(infected_days); ie <- max(infected_days)
      list(infection_time = it, infectious_start = it, infectious_end = ie, detection_time = it, infection_resolved = ie + 1L)
    }
  }, by = .(hh_id = HH, person_id = individual_ID, role)]
  person <- as.data.frame(person, stringsAsFactors = FALSE)

  vl_list <- vector("list", nrow(person))
  if (vl_source == "none") {
    for (i in seq_len(nrow(person))) {
      it <- person$infection_time[i]; ir <- person$infection_resolved[i]
      vl_list[[i]] <- if (is.na(it) || is.na(ir) || ir < it) numeric(0) else rep(-5, ir - it + 1L)
    }
  } else if (vl_source == "column") {
    if (is.null(vl_column) || !vl_column %in% names(dt)) stop("vl_column not found.", call. = FALSE)
    max_day <- max(dt$day, na.rm = TRUE)
    for (i in seq_len(nrow(person))) {
      sub <- dt[dt$HH == person$hh_id[i] & dt$individual_ID == person$person_id[i], ]
      traj <- rep(-5, max_day)
      if (nrow(sub)) { idx <- sub$day; vals <- sub[[vl_column]]; ok <- !is.na(idx) & idx >= 1 & idx <= max_day; traj[idx[ok]] <- vals[ok] }
      vl_list[[i]] <- traj
    }
  } else if (vl_source == "simulate") {
    for (i in seq_len(nrow(person))) {
      it <- person$infection_time[i]; ir <- person$infection_resolved[i]
      if (is.na(it) || is.na(ir) || ir < it) { vl_list[[i]] <- numeric(0) }
      else { t_rel <- 0:as.integer(ir - it); rp <- draw_random_VL_params(person$role[i])
      vl_list[[i]] <- simulate_viral_load_trajectory(t_rel, rp$v_p, rp$t_p, rp$lambda_g, rp$lambda_d) }
    }
  }
  person$vl_full_trajectory <- vl_list
  person$hh_id <- as.character(person$hh_id); person$person_id <- as.integer(person$person_id)
  person$role <- factor(person$role, levels = role_levels)
  person <- person[order(person$hh_id, person$person_id), , drop = FALSE]
  split(person, person$hh_id, drop = TRUE)
}


#' Prepare Stan data from diagnostic testing records
#'
#' Converts raw diagnostic testing data into Stan-ready format with the
#' \code{time_since_infection} matrix for soft-gate latent period modeling.
#'
#' @param raw_data Data frame with columns \code{hh_id}, \code{person_id}, \code{role},
#'   \code{day_index}, \code{test_result}, and optionally \code{pcr_sample}.
#' @param seasonal_forcing_list Named list with seasonal forcing vectors by role.
#' @param viral_testing_mode Character; \code{"viral load"} or \code{"Ct"}.
#' @param T_max Integer; maximum time horizon (days).
#' @param V_ref,V_rho Numeric; viral load transmission parameters.
#' @param peak_day,width Numeric; infectivity profile parameters.
#' @param alpha_comm_by_role Numeric; community acquisition rate.
#' @param reference_phi,reference_kappa Numeric; reference parameters.
#'
#' @return Named list suitable for Stan.
prepare_stan_data <- function(raw_data, seasonal_forcing_list = NULL, viral_testing_mode = "viral load",
                              T_max = 365, V_ref = 3, V_rho = 2.5, peak_day = 2, width = 4,
                              alpha_comm_by_role = 5e-3, reference_phi = 1, reference_kappa = 1) {

  raw_data <- raw_data %>% dplyr::arrange(hh_id, person_id, day_index)
  all_hh <- unique(raw_data$hh_id)
  proper_order <- all_hh[order(as.integer(gsub("\\D", "", all_hh)))]

  individuals <- raw_data %>%
    dplyr::group_by(hh_id, person_id, role) %>% dplyr::slice(1) %>% dplyr::ungroup() %>%
    dplyr::mutate(stan_id = dplyr::row_number(), hh_factor = as.integer(factor(hh_id, levels = proper_order))) %>%
    dplyr::arrange(hh_factor, stan_id)

  N <- nrow(individuals); H <- max(individuals$hh_factor)
  role_levels <- c("adult", "child", "toddler", "elderly")

  I <- matrix(0L, N, T_max); V <- matrix(0, N, T_max); V_term <- matrix(0, N, T_max)
  time_since_infection <- matrix(-1L, N, T_max)
  if (viral_testing_mode == "Ct") V[] <- 45

  for (i in 1:N) {
    p_id <- individuals$person_id[i]; h_id <- individuals$hh_id[i]; role <- individuals$role[i]
    history <- raw_data %>% dplyr::filter(hh_id == h_id, person_id == p_id)
    pos_days <- history$day_index[history$test_result == 1]
    neg_days <- history$day_index[history$test_result == 0]

    if (length(pos_days) > 0) {
      first_pos <- min(pos_days); last_pos <- max(pos_days)
      prev_negs <- neg_days[neg_days < first_pos]
      if (length(prev_negs) > 0) {
        last_neg <- max(prev_negs); possible_starts <- seq(last_neg + 1, first_pos)
        t_start <- if (length(possible_starts) == 1) possible_starts else sample(possible_starts, 1)
      } else {
        possible_starts <- seq(max(1, first_pos - 7), first_pos)
        t_start <- sample(possible_starts, 1)
      }
      next_negs <- neg_days[neg_days > last_pos]
      if (length(next_negs) > 0) {
        first_neg_after <- min(next_negs); possible_ends <- seq(last_pos, first_neg_after - 1)
        t_end <- if (length(possible_ends) == 1) possible_ends else sample(possible_ends, 1)
      } else {
        possible_ends <- seq(last_pos, min(T_max, last_pos + 14)); t_end <- sample(possible_ends, 1)
      }

      if (t_start <= T_max) I[i, t_start] <- 1L
      if (t_start <= T_max) for (t_idx in t_start:T_max) time_since_infection[i, t_idx] <- t_idx - t_start

      vs <- max(1, t_start); ve <- min(T_max, t_end)
      if (vs <= ve) {
        rel_days <- (vs:ve) - t_start
        if (viral_testing_mode == "viral load") {
          p <- default_VL_params[[role]]
          traj_vals <- simulate_viral_load_trajectory(rel_days, p$v_p, p$t_p, p$lambda_g, p$lambda_d)
          traj_vals[traj_vals < 0] <- 0
        } else {
          p <- default_Ct_params[[role]]
          traj_vals <- simulate_Ct_trajectory(rel_days, p$Cpeak, p$r, p$d, p$t_peak)
          traj_vals[traj_vals > 45] <- 45
        }
        V[i, vs:ve] <- traj_vals
      }
    }
  }

  for (i in 1:N) for (t in 1:T_max) {
    val <- V[i, t]
    V_term[i, t] <- if (viral_testing_mode == "viral load" && val > 0) (val / V_ref)^V_rho else 0
  }

  hh_members_list <- split(individuals$stan_id, individuals$hh_factor)
  hh_size <- sapply(hh_members_list, length); hh_max_size <- max(hh_size)
  hh_members <- matrix(1L, nrow = H, ncol = hh_max_size)
  for (h in 1:H) hh_members[h, 1:hh_size[h]] <- hh_members_list[[h]]

  if (is.null(seasonal_forcing_list)) {
    seasonal_forcing_list <- list(adult = rep(1, T_max), child = rep(1, T_max), toddler = rep(1, T_max), elderly = rep(1, T_max))
  }
  seasonal_forcing_mat <- matrix(0, nrow = T_max, ncol = 4)
  seasonal_forcing_mat[, 1] <- seasonal_forcing_list$adult[1:T_max]
  seasonal_forcing_mat[, 2] <- seasonal_forcing_list$child[1:T_max]
  seasonal_forcing_mat[, 3] <- seasonal_forcing_list$toddler[1:T_max]
  seasonal_forcing_mat[, 4] <- seasonal_forcing_list$elderly[1:T_max]

  g_profile <- g_rescaled(1:T_max, peak_day = peak_day, width = width)

  list(N = N, T = T_max, H = H, R = length(role_levels), delta = 0,
       hh_id = individuals$hh_factor, role_id = match(individuals$role, role_levels),
       I = I, time_since_infection = time_since_infection, V_term = V_term,
       hh_size = as.integer(hh_size), hh_max_size = as.integer(hh_max_size), hh_members = hh_members,
       alpha_comm_by_role = alpha_comm_by_role, g_profile = as.numeric(g_profile),
       seasonal_forcing_mat = seasonal_forcing_mat, reference_phi = reference_phi, reference_kappa = reference_kappa)
}


#' Build Stan data arrays from household list
#'
#' Converts a list of per-household data frames into Stan-ready arrays.
#'
#' @param households List of per-household data frames.
#' @param T_max Integer; maximum time horizon.
#' @param seasonal_forcing_list Named list of seasonal forcing vectors.
#' @param alpha_comm_by_role,beta1,beta2,V_ref Numeric; model parameters.
#' @param reference_phi,reference_kappa Numeric; reference parameters.
#' @param g_peak_day,g_width Numeric; infectivity profile parameters.
#'
#' @return Named list for Stan.
#' @export
build_stan_household_arrays <- function(households, T_max = 365L, seasonal_forcing_list = NULL,
                                        alpha_comm_by_role = 5e-3, beta1 = 0.2, beta2 = 0.6, V_ref = 1e3,
                                        reference_phi = 1, reference_kappa = 1, g_peak_day = 2, g_width = 4) {

  .clamp_int <- function(x, lo, hi) as.integer(pmax(lo, pmin(hi, x)))

  .ensure_sf_mat <- function(sf_list, T_max) {
    role_levels <- c("adult", "child", "elderly", "toddler")
    if (is.null(sf_list)) sf_list <- setNames(lapply(role_levels, function(r) rep(1, T_max)), role_levels)
    mat <- matrix(1, nrow = T_max, ncol = 4); colnames(mat) <- role_levels
    for (r in role_levels) if (r %in% names(sf_list)) mat[, r] <- rep_len(sf_list[[r]], T_max)
    mat
  }

  all_individuals <- do.call(rbind, lapply(households, function(hh) if (is.null(hh) || nrow(hh) == 0) NULL else hh))
  if (is.null(all_individuals) || nrow(all_individuals) == 0) stop("Stan data builder: no individuals found.")

  if (!("hh_id" %in% names(all_individuals))) {
    if ("HH" %in% names(all_individuals)) all_individuals$hh_id <- all_individuals$HH
    else stop("Stan data builder: need `hh_id` or `HH`.")
  }
  if (!("role" %in% names(all_individuals))) stop("Stan data builder: need `role` column.")
  all_individuals$role <- .norm_role(all_individuals$role)
  role_levels <- c("adult", "child", "elderly", "toddler")

  for (nm in c("infection_time", "infectious_start", "infectious_end", "infection_resolved")) {
    if (!nm %in% names(all_individuals)) all_individuals[[nm]] <- NA_integer_
    else all_individuals[[nm]] <- suppressWarnings(as.integer(all_individuals[[nm]]))
  }

  all_individuals$hh_id <- factor(all_individuals$hh_id)
  hh_id <- as.integer(all_individuals$hh_id)
  N <- nrow(all_individuals); H <- length(levels(all_individuals$hh_id)); T <- as.integer(T_max)
  role_id <- match(all_individuals$role, role_levels)

  I <- matrix(0L, nrow = N, ncol = T); time_since_infection <- matrix(-1L, nrow = N, ncol = T)
  for (i in seq_len(N)) {
    s <- all_individuals$infection_time[i]
    if (!is.na(s) && s >= 1L && s <= T) {
      I[i, as.integer(s)] <- 1L
      for (t_idx in s:T) time_since_infection[i, t_idx] <- t_idx - s
    }
  }

  Y <- matrix(0L, nrow = N, ncol = T)
  for (i in seq_len(N)) {
    s <- all_individuals$infectious_start[i]; e <- all_individuals$infectious_end[i]
    if (!is.na(s) && !is.na(e)) {
      s_idx <- .clamp_int(s, 1L, T); e_idx <- .clamp_int(e, 1L, T)
      if (e_idx >= s_idx) Y[i, s_idx:e_idx] <- 1L
    }
  }

  V <- matrix(0, nrow = N, ncol = T); V_term <- matrix(0, nrow = N, ncol = T)
  if ("vl_full_trajectory" %in% names(all_individuals)) {
    for (i in seq_len(N)) {
      traj <- all_individuals$vl_full_trajectory[[i]]
      if (is.null(traj) || !length(traj)) next
      inf_time <- all_individuals$infection_time[i]; inf_end <- all_individuals$infection_resolved[i]
      if (is.na(inf_time) || is.na(inf_end)) next
      days_seq <- seq.int(inf_time, inf_end)
      valid_len <- min(length(traj), length(days_seq), T - as.integer(inf_time) + 1L)
      if (valid_len < 1L) next
      days_fill <- as.integer(days_seq[seq_len(valid_len)]); vals_fill <- as.numeric(traj[seq_len(valid_len)])
      vals_fill[!is.finite(vals_fill) | vals_fill < 0] <- 0
      keep <- days_fill >= 1L & days_fill <= T
      if (any(keep)) V[i, days_fill[keep]] <- vals_fill[keep]
    }
  }
  if (any(V > 0)) {
    n_exp <- 3.91; p_v <- 0.83e-2; K <- 5.39; Vm <- 6.59; sai <- 4.57
    powV <- V^n_exp; denom <- powV + K^n_exp
    idx <- denom > 0; V_term[idx] <- 1 - exp(-sai * Vm * p_v * powV[idx] / denom[idx])
  }

  hh_members_list <- split(seq_len(N), hh_id); hh_size <- vapply(hh_members_list, length, integer(1)); hh_max_size <- max(hh_size)
  hh_members <- matrix(1L, nrow = H, ncol = hh_max_size)
  for (h in seq_len(H)) hh_members[h, 1:hh_size[h]] <- hh_members_list[[h]]

  g_profile <- exp(-0.5 * ((seq_len(T) - g_peak_day) / g_width)^2); g_profile <- g_profile / max(g_profile)
  sf_mat <- .ensure_sf_mat(seasonal_forcing_list, T_max = T)
  sf_mat <- sf_mat[, c("adult", "child", "elderly", "toddler"), drop = FALSE]

  list(N = N, T = T, H = H, R = length(role_levels), delta = 0, hh_id = hh_id, role_id = as.integer(role_id),
       Y = Y, I = I, V = V, V_term = V_term, time_since_infection = time_since_infection,
       alpha_comm_by_role = alpha_comm_by_role,
       hh_size = as.integer(hh_size), hh_max_size = as.integer(hh_max_size), hh_members = hh_members,
       g_profile = as.numeric(g_profile), V_ref = V_ref, beta1 = beta1, beta2 = beta2,
       reference_phi = reference_phi, reference_kappa = reference_kappa, seasonal_forcing_mat = sf_mat)
}


#' Run the Stan household model
#'
#' Wrapper around \code{rstan::sampling()} with reasonable defaults.
#'
#' @param stan_data Named list from \code{build_stan_household_arrays} or \code{prepare_stan_data}.
#' @param stan_file Path to Stan model file.
#' @param chains,iter,warmup Stan MCMC settings.
#' @param control List passed to \code{rstan::sampling()}.
#' @param init,refresh,cores Stan parameters.
#' @param stan_code Optional character string with Stan program code.
#' @param package_name Optional package name for \code{system.file} lookup.
#'
#' @return A \code{stanfit} object.
run_household_stan <- function(stan_data, stan_file = "model.stan", chains = 4, iter = 2000, warmup = 1000,
                               control = list(adapt_delta = 0.99, max_treedepth = 20),
                               init = "random", refresh = 50, cores = 4, stan_code = NULL, package_name = NULL) {

  res <- .resolve_stan_model(stan_file = stan_file, stan_code = stan_code, package_name = package_name)

  if (requireNamespace("rstan", quietly = TRUE)) {
    sm <- if (!is.null(res$code)) rstan::stan_model(model_code = res$code, verbose = FALSE)
    else rstan::stan_model(file = res$file, verbose = FALSE)
    op <- options(mc.cores = cores); on.exit(options(op), add = TRUE)
    fit <- rstan::sampling(sm, data = stan_data, chains = chains, iter = iter, warmup = warmup,
                           control = control, init = init, refresh = refresh)
    return(fit)
  } else { stop("Package 'rstan' not available.") }
}


#' Resolve Stan model file path
#' @keywords internal
.resolve_stan_model <- function(stan_file, stan_code = NULL, package_name = NULL) {
  if (!is.null(stan_code)) return(list(file = NULL, code = stan_code))
  if (!is.null(package_name)) {
    pkg_path <- system.file("stan", stan_file, package = package_name)
    if (nzchar(pkg_path) && file.exists(pkg_path)) return(list(file = pkg_path, code = NULL))
  }
  if (file.exists(stan_file)) return(list(file = stan_file, code = NULL))
  inst_path <- system.file("stan", stan_file, package = utils::packageName())
  if (nzchar(inst_path) && file.exists(inst_path)) return(list(file = inst_path, code = NULL))
  stop("Stan file not found: ", stan_file, call. = FALSE)
}


#' Tidy posterior summary for role multipliers
#'
#' Extracts posterior summaries from a \code{stanfit}.
#'
#' @param fit A \code{stanfit} object.
#' @return A data frame with posterior summaries.
#' @export
postprocess_stan_fit <- function(fit) {
  s <- tryCatch(rstan::summary(fit)$summary, error = function(e) NULL)
  if (is.null(s) || !nrow(s)) return(data.frame(Parameter = character()))
  tab <- as.data.frame(s); tab$Parameter <- rownames(tab); rownames(tab) <- NULL
  names(tab) <- sub("^X2[.]5[.]$", "2.5%", names(tab))
  names(tab) <- sub("^X50[.]$", "50%", names(tab))
  names(tab) <- sub("^X97[.]5[.]$", "97.5%", names(tab))
  keep <- intersect(c("Parameter", "mean", "sd", "2.5%", "50%", "97.5%"), names(tab))
  tab[, keep, drop = FALSE]
}
