#' Convert a user data frame to a list of household tables
#'
#' Splits a long-format test table into a list of per-household \code{data.table}s,
#' preserving the required columns and (optionally) any extra covariates.
#'
#' @param df Data frame with at least household ID, individual ID, role, test date,
#'   infection status, and community risk columns.
#' @param hh_col,id_col,role_col,date_col,inf_col,comm_col Character. Column names
#'   for household ID (default \code{"HH"}), individual ID (\code{"individual_ID"}),
#'   role (\code{"role"}), test date (\code{"test_date"}), infection status
#'   (\code{"infection_status"}), and community risk (\code{"community_risk"}).
#' @param keep_extra_cols Logical; keep additional user columns (default \code{TRUE}).
#'
#' @return List of \code{data.table}s, one per household (household ID placed first).
dataframe_to_household_list <- function(df,
                                        hh_col   = "HH",
                                        id_col   = "individual_ID",
                                        role_col = "role",
                                        date_col = "test_date",
                                        inf_col  = "infection_status",
                                        comm_col = "community_risk",

                                        # keep any additional user covariates by default
                                        keep_extra_cols = TRUE) {

  stopifnot(is.data.frame(df))
  DT <- data.table::as.data.table(df)

  req <- c(hh_col, id_col, role_col, date_col, inf_col, comm_col)
  miss <- setdiff(req, names(DT))
  if (length(miss)) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  DT[[hh_col]]   <- as.integer(DT[[hh_col]])
  DT[[id_col]]   <- as.integer(DT[[id_col]])
  DT[[role_col]] <- as.character(DT[[role_col]])

  if (inherits(DT[[date_col]], "Date")) {
    origin <- min(DT[[date_col]], na.rm = TRUE)
    DT[[date_col]] <- as.integer(DT[[date_col]] - origin) + 1L
  } else {
    DT[[date_col]] <- as.integer(DT[[date_col]])
  }
  DT[[inf_col]]  <- as.integer(DT[[inf_col]])
  DT[[comm_col]] <- as.numeric(DT[[comm_col]])

  #Trim columns to keep
  base_cols <- req
  extra_cols <- setdiff(names(DT), base_cols)
  cols <- if (keep_extra_cols) c(base_cols, extra_cols) else base_cols

  data.table::setorderv(DT, c(hh_col, id_col, date_col))
  split_list <- split(DT[, ..cols], by = hh_col, keep.by = FALSE, drop = TRUE)

  hh_ids <- as.integer(names(split_list))
  for (i in seq_along(split_list)) {
    split_list[[i]][, (hh_col) := hh_ids[i]]
    data.table::setcolorder(split_list[[i]], c(hh_col, setdiff(names(split_list[[i]]), hh_col)))
  }
  return(split_list)
}



#' Summarize individual-level infection data
#'
#' Produces one row per individual with detection windows, inferred infectious
#' windows (relative days), index-case flags, observation bounds, role/age
#' classification, and aggregated covariates. Optionally builds day-series
#' list-columns for selected covariates.
#'
#' @param raw_dt List of household-level data frames/tables.
#' @param study_start,study_end \code{Date}. Analysis window defining day indices.
#' @param day_series_covariates Logical; add day-series list-columns (default \code{TRUE}).
#' @param series_cols Character or \code{NULL}; covariates to series-encode (default \code{NULL} = all).
#'
#' @return \code{data.table} with one row per individual containing detection and
#'   infectious windows, flags, observation bounds, role/age, and covariate summaries.
summarize_individuals <- function(raw_dt, study_start, study_end,
                                  day_series_covariates = TRUE,
                                  series_cols = NULL
) {
  stopifnot(inherits(study_start, "Date"), inherits(study_end, "Date"))
  time.steps <- as.integer(study_end - study_start)  # tmax used for day-series length

  #Bind all households into long format rows (test records)
  raw_long <- data.table::rbindlist(raw_dt, use.names = TRUE, fill = TRUE)

  # inside summarize_individuals(), just after you set raw_long:
  raw_long[, role := tolower(role)]
  # include common aliases
  raw_long[role %in% c("parent","adult"), role := "adult"]
  raw_long[role %in% c("elder","elderly"), role := "elder"]
  raw_long[role %in% c("child","sibling"), role := "sibling"]
  raw_long[role %in% c("infant","baby","toddler"), role := "infant"]

  data.table::setDT(raw_long)

  #Core columns that are NOT covariates
  core_cols <- c("HH","individual_ID","role","test_date","infection_status","community_risk")

  #Normalize covariate column names
  .normalize <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_|_$", "", x)
    tolower(x)
  }
  covar_raw <- setdiff(names(raw_long), core_cols)
  if (length(covar_raw)) {
    data.table::setnames(raw_long, covar_raw, .normalize(covar_raw))
  }

  #Candidate covariates are everything outside core columns with at least one non NA
  covar_cols <- setdiff(names(raw_long), core_cols)
  covar_cols <- covar_cols[vapply(covar_cols, function(x) any(!is.na(raw_long[[x]])), logical(1))]

  #If the user specified a subset for day-series, normalize those names and intersect
  series_cols_norm <- NULL
  if (!is.null(series_cols)) {
    series_cols_norm <- .normalize(series_cols)
    series_cols_norm <- intersect(series_cols_norm, covar_cols)
  }

  #Helpers
  .mode_chr <- function(v) {
    v <- v[!is.na(v)]
    if (!length(v)) return(NA_character_)
    names(sort(table(v), decreasing = TRUE))[1]
  }
  .mode_num <- function(v) {
    v <- v[!is.na(v)]
    if (!length(v)) return(NA_real_)
    as.numeric(names(sort(table(v), decreasing = TRUE))[1])
  }

  #Aggregate covariates per individual (mean/mode/first/last + timevarying)
  data.table::setorder(raw_long, HH, individual_ID, test_date)

  covar_summary <- NULL
  if (length(covar_cols) > 0) {
    covar_summary <- raw_long[, {
      out <- vector("list", 0L)
      for (cn in covar_cols) {
        v <- .SD[[cn]]
        is_num   <- is.numeric(v)
        timevary <- length(unique(v[!is.na(v)])) > 1

        if (is_num) {
          nz <- which(!is.na(v))
          out[[paste0(cn, "_first")]] <- if (length(nz)) v[nz[1]] else NA_real_
          out[[paste0(cn, "_mean")]]  <- mean(v, na.rm = TRUE)
          out[[paste0(cn, "_last")]]  <- if (length(nz)) v[nz[length(nz)]] else NA_real_
        } else {
          if (is.logical(v)) {
            out[[paste0(cn, "_mode")]] <- as.integer(.mode_num(as.numeric(v)))
          } else {
            out[[paste0(cn, "_mode")]] <- .mode_chr(as.character(v))
          }
        }
        out[[paste0(cn, "_timevarying")]] <- timevary
      }
      out
    }, by = .(HH, individual_ID, role)]
  }

  ##########################################################################
  #                                                                        #
  #            Original summarization of infection windows                 #
  #                                                                        #
  ##########################################################################
  raw_sum <- data_summarization(raw_long)  # uses normalized names; produces per-person rows
  data.table::setDT(raw_sum)
  data.table::setorder(raw_sum, HH)

  dt <- data.table::copy(raw_sum)
  data.table::setnames(dt, "individual_ID", "indiv.index")

  conv <- function(x) data.table::fifelse(is.na(x), as.Date(NA), study_start + as.integer(x))
  dt[, T_FP_date      := conv(infection.detected.start)]
  dt[, T_LP_date      := conv(infection.detected.end)]
  dt[, last_neg_date  := conv(last_negative)]
  dt[, inf_date       := as.Date(NA)]
  dt[, inf_start_date := T_FP_date]
  dt[, inf_end_date   := T_LP_date]

  parse_vec <- function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_integer_)
    vals <- suppressWarnings(as.integer(strsplit(s, ",")[[1]]))
    vals <- vals[!is.na(vals)]
    if (!length(vals)) return(NA_integer_)
    vals
  }
  dt[, `:=`(
    inf_win_start = sapply(infection.infectious.day, function(s) {
      v <- parse_vec(s); if (all(is.na(v))) NA_integer_ else min(v, na.rm = TRUE)
    }),
    inf_win_end = sapply(infection.infectious.day, function(s) {
      v <- parse_vec(s); if (all(is.na(v))) NA_integer_ else max(v, na.rm = TRUE)
    })
  )]
  dt[is.finite(inf_win_start), inf_start_date := study_start + inf_win_start]
  dt[is.finite(inf_win_end),   inf_end_date   := study_start + inf_win_end]

  #Infection flags and index
  dt[, infected := !is.na(T_FP_date)]
  dt[, is_index := FALSE]
  dt[infected == TRUE, is_index := T_FP_date == min(T_FP_date, na.rm = TRUE), by = HH]

  #Role -> age category
  dt[, age_cat :=
       data.table::fifelse(role == "infant",  1L,
                           data.table::fifelse(role == "sibling", 2L,
                                               data.table::fifelse(role == "adult",   3L,
                                                                   data.table::fifelse(role == "elder",   4L, NA_integer_))))]
  if (any(is.na(dt$age_cat))) stop("Unknown role labels remain, extend the role to age mapping.")

  #Relative day columns
  dt[, `:=`(
    inf_day_rl            = as.integer(inf_date       - study_start),
    infectious_day_rl     = as.integer(inf_start_date - study_start),
    infectious_end_day_rl = as.integer(inf_end_date   - study_start)
  )]
  dt[infected == FALSE,
     `:=`(inf_day_rl = NA_integer_,
          infectious_day_rl = NA_integer_,
          infectious_end_day_rl = NA_integer_)]

  #IDs & index recompute on new IDs
  data.table::setnames(dt, "HH", "ID_hh")
  dt[, ID_indiv := sprintf("HH%03d_%02d", ID_hh, indiv.index)]
  dt[, is_index := FALSE]
  dt[infected == TRUE, is_index := T_FP_date == min(T_FP_date, na.rm = TRUE), by = ID_hh]

  #Observation window
  dt[, `:=`(obs_start_date = study_start, obs_end_date = study_end)]

  #Merge per-person aggregate covariates (mean/mode/first/last/timevarying)
  if (!is.null(covar_summary)) {
    data.table::setDT(covar_summary)
    key_cols <- c("HH","individual_ID","role")
    cov_cols <- setdiff(names(covar_summary), key_cols)
    add_cols <- setdiff(cov_cols, names(dt))
    if (length(add_cols)) {
      cs_reduced <- covar_summary[, c(key_cols, add_cols), with = FALSE]
      dt <- merge(
        dt, cs_reduced,
        by.x = c("ID_hh","indiv.index","role"),
        by.y = c("HH","individual_ID","role"),
        all.x = TRUE,
        sort  = FALSE
      )
    }
  }

  ##########################################################################
  #                                                                        #
  #            Build day-series list-columns for covariates                #
  #                                                                        #
  ##########################################################################
  if (day_series_covariates && length(covar_cols) > 0) {
    #Choose which columns to build series for
    build_cols <- if (!is.null(series_cols_norm)) series_cols_norm else covar_cols

    #For faster access:
    key_cols <- c("ID_hh","indiv.index")
    #We need the raw_long rows to map day indices for each individual
    raw_long[, day := as.integer(test_date)]  # already relative days in your setup

    #LOCF helper (works for numeric or character vectors)
    locf_fill <- function(x) {
      if (all(is.na(x))) return(x)
      # forward fill
      for (k in seq_along(x)) if (is.na(x[k]) && k > 1L) x[k] <- x[k-1L]
      # backfill head if still NA
      if (is.na(x[1L])) x[1L] <- x[which(!is.na(x))[1L]]
      # if still any gap (all NA originally), just return x
      x
    }

    #Build series per individual
    data.table::setkeyv(raw_long, c("HH","individual_ID","role","day"))

    #Map HH/individual_ID -> ID_hh/indiv.index using dt
    id_map <- dt[, .(ID_hh, indiv.index, role)]
    setnames(id_map, c("ID_hh","indiv.index"), c("HH","individual_ID"))

    #Join to get consistent (HH, individual_ID, role) combos
    raw_long2 <- merge(raw_long, id_map, by = c("HH","individual_ID","role"), all.x = TRUE)

    #For each individual, construct a tmax+1 vector for each build_cols
    series_dt <- raw_long2[, {
      out <- vector("list", 0L)
      days <- day
      for (cn in build_cols) {
        #Initialize vector of length tmax+1 with NA (numeric or character, depending on observed type)
        vals_obs <- .SD[[cn]]
        is_num <- is.numeric(vals_obs)
        if (is_num) {
          s <- rep(NA_real_, time.steps + 1L)
        } else {
          #Use character for categorical; factors can be coerced later
          s <- rep(NA_character_, time.steps + 1L)
        }
        #Place observed values on their observed days (indices are day+1)
        good <- !is.na(days) & !is.na(vals_obs) &
          (days >= 0L) & (days <= time.steps)
        if (any(good)) {
          s[days[good] + 1L] <- vals_obs[good]
        }
        #LOCF + backfill to avoid leading NA
        s <- locf_fill(s)
        out[[cn]] <- list(s)  # list-column cell
      }
      out
    }, by = .(HH, individual_ID, role), .SDcols = build_cols]

    #Merge series back onto dt by (ID_hh, indiv.index, role)
    if (nrow(series_dt)) {
      dt <- merge(
        dt, series_dt,
        by.x = c("ID_hh","indiv.index","role"),
        by.y = c("HH","individual_ID","role"),
        all.x = TRUE,
        sort = FALSE
      )
    }
  }
  data.table::setorder(dt, ID_hh, indiv.index)
  dt[]
}



#' Impute infection timelines from delay distributions
#'
#' Imputes infection date, infectious start/end, and component delays using
#' gamma distributions, optionally scaled by covariate functions.
#'
#' @param dt \code{data.table} from \code{\link{summarize_individuals}}.
#' @param study_start \code{Date} origin for relative day indices.
#' @param latent_par,report_par,infect_par Lists with \code{shape} and \code{scale}
#'   for latent, reporting, and infectious periods.
#' @param latent_scale_fn,report_scale_fn,infect_scale_fn Optional functions taking
#'   \code{dt[idx]} (infected rows) and returning numeric scale multipliers.
#'
#' @return The input \code{dt} with columns \code{latent_delay}, \code{report_delay},
#'   \code{infect_period}, \code{inf_date}, \code{inf_start_date}, \code{inf_end_date},
#'   and relative-day variants.
infectious_time_imputation <- function(dt, study_start,
                                       latent_par, report_par, infect_par,

                                       # For imputing delays to infection depending on the covariates
                                       # For example,shorter infectious period if vaccinated, or delays depending on antibody-level
                                       latent_scale_fn = NULL,
                                       report_scale_fn = NULL,
                                       infect_scale_fn = NULL) {
  dt <- data.table::copy(dt)

  #Ensure columns exist
  for (nm in c("latent_delay","report_delay","infect_period",
               "inf_date","inf_start_date","inf_end_date",
               "T_FP_date","T_LP_date","last_neg_date","infected")) {
    if (!nm %in% names(dt)) dt[, (nm) := NA]
  }

  #Vectorized truncated-gamma sampler (allows vector upper & scale)
  rtrunc_gamma_vec <- function(shape, scale, upper){
    n <- length(upper)
    up <- upper; up[!is.finite(up) | up <= 0] <- 1e-8
    u  <- runif(n)
    qgamma(u * pgamma(up, shape = shape, rate = 1/scale),
                  shape = shape, rate = 1/scale)
  }

  #Indices of infected rows
  idx <- which(isTRUE(dt$infected))
  if (length(idx) == 0) return(dt)

  #Max look-back window to first positive (or 14 if no prior negative recorded)
  max_back <- ifelse(is.na(dt$last_neg_date[idx]),
                     14,
                     as.integer(dt$T_FP_date[idx] - dt$last_neg_date[idx]))

  #Per-row scale multipliers (default 1)
  m_lat <- if (is.null(latent_scale_fn)) rep(1, length(idx)) else as.numeric(latent_scale_fn(dt[idx]))
  m_rep <- if (is.null(report_scale_fn)) rep(1, length(idx)) else as.numeric(report_scale_fn(dt[idx]))
  m_inf <- if (is.null(infect_scale_fn)) rep(1, length(idx)) else as.numeric(infect_scale_fn(dt[idx]))

  #Draw latent delays
  lat <- rtrunc_gamma_vec(shape = latent_par$shape,
                          scale = latent_par$scale * m_lat,
                          upper = pmax(max_back, 1e-8))

  #Draw reporting delays with remaining headroom
  rep_upper <- pmax(max_back - lat, 1e-8)
  repd <- rtrunc_gamma_vec(shape = report_par$shape,
                           scale = report_par$scale * m_rep,
                           upper = rep_upper)

  #Draw infectious period
  infd <- rgamma(length(idx), shape = infect_par$shape, rate = 1/infect_par$scale) * m_inf
  infd[infd < 0] <- 0

  dt$latent_delay[idx] <- lat
  dt$report_delay[idx] <- repd
  dt$infect_period[idx] <- infd

  #Infection date = first positive - (latent + report), but not before last negative + 1
  inf_date_i <- dt$T_FP_date[idx] - (lat + repd)
  has_neg    <- !is.na(dt$last_neg_date[idx])
  inf_date_i[has_neg] <- pmax(inf_date_i[has_neg], dt$last_neg_date[idx][has_neg] + 1)

  dt$inf_date[idx] <- inf_date_i
  dt$inf_start_date[idx] <- dt$inf_date[idx] + lat

  #Cap end by last positive date if available, otherwise by observed window (if present)
  if ("obs_end_date" %in% names(dt)) {
    cap_end <- ifelse(is.na(dt$T_LP_date[idx]), dt$obs_end_date[idx], dt$T_LP_date[idx])
  } else {
    cap_end <- dt$T_LP_date[idx]
  }
  dt$inf_end_date[idx] <- pmin(dt$inf_start_date[idx] + infd, cap_end, na.rm = TRUE)

  #Relative-day indices
  dt$inf_day_rl[idx]            <- as.integer(dt$inf_date[idx]       - study_start)
  dt$infectious_day_rl[idx]     <- as.integer(dt$inf_start_date[idx] - study_start)
  dt$infectious_end_day_rl[idx] <- as.integer(dt$inf_end_date[idx]   - study_start)

  return(dt)
}



#' Construct person-day long data
#'
#' Expands individual timelines into daily rows and computes within-household
#' infectious counts by infector role for likelihood-based estimation.
#'
#' @param dt \code{data.table} from \code{\link{infectious_time_imputation}}.
#' @param tmax Integer; maximum day index.
#' @param cases_t Numeric of length \code{tmax + 1}; community intensity for days 0..\code{tmax}.
#' @param covariate_cols Character; names of covariates to carry (scalars or day series).
#'
#' @return \code{data.table} with columns:
#'   \code{agegrp2}, \code{agegrp3}, \code{agegrp4}, \code{n_inf},
#'   \code{n_inf_infant}, \code{n_inf_sibling}, \code{n_inf_adult}, \code{n_inf_elder},
#'   \code{cases}, \code{event}, \code{ID_indiv}, \code{ID_hh}, \code{day}, and requested covariates.
build_person_day_table <- function(dt, tmax, cases_t, covariate_cols = character(0)) {
  dt <- data.table::copy(dt)

  #Sanity checks
  if (length(cases_t) != (tmax + 1L))
    stop("`cases_t` must have length tmax + 1.")
  req <- c("ID_hh","ID_indiv","role","age_cat","infected","is_index",
           "inf_day_rl","infectious_day_rl","infectious_end_day_rl")
  miss <- setdiff(req, names(dt))
  if (length(miss)) stop("Missing required columns in `dt`: ", paste(miss, collapse = ", "))

  #Helper: convert span [a,b] to 1-based indices for daily arrays of length tmax+1
  safe_span_idx <- function(a, b, tmax) {
    if (length(a) == 0L || length(b) == 0L) return(integer(0))
    if (is.na(a) || is.na(b)) return(integer(0))
    a <- max(as.integer(a), 0L); b <- min(as.integer(b), tmax)
    if (a > b) return(integer(0))
    (a:b) + 1L
  }

  #Helper: extract the value of covariate `col` for a given record `rec` on day `d`
  #Supports:
  #  - scalar columns (recycled across days)
  #  - list-columns whose cell is a vector indexed by day+1 (numeric or character)
  cov_value_for_day <- function(rec, col, d, default = NA) {
    if (!col %in% names(rec)) return(default)
    x <- rec[[col]]
    # If it's a list-column, pull the first (and only) cell content
    if (is.list(x)) x <- x[[1]]
    # If it's a factor, drop to character to keep day-specific levels stable
    if (is.factor(x)) x <- as.character(x)
    # If it's a vector (length > 1), treat as a day-series and index day+1
    if (length(x) > 1L) {
      idx <- d + 1L
      if (idx >= 1L && idx <= length(x)) {
        return(x[[idx]])
      } else {
        return(default)
      }
    }
    #Otherwise, scalar: recycle across days
    x[1]
  }

  rows <- list()
  for (hh in unique(dt$ID_hh)) {
    hhdat <- dt[ID_hh == hh]

    #Per-day infectious counts by infector role (people counts)
    n_inf_total   <- integer(tmax + 1L)
    n_inf_infant  <- integer(tmax + 1L)
    n_inf_sibling <- integer(tmax + 1L)
    n_inf_adult   <- integer(tmax + 1L)
    n_inf_elder   <- integer(tmax + 1L)

    #Mark infectious spans for each infected household member
    for (j in hhdat[infected == TRUE]$ID_indiv) {
      rec <- hhdat[ID_indiv == j]
      idx <- safe_span_idx(rec$infectious_day_rl, rec$infectious_end_day_rl, tmax)
      if (length(idx)) {
        n_inf_total[idx] <- n_inf_total[idx] + 1L
        ac <- as.integer(rec$age_cat[1L])
        if (ac == 1L) n_inf_infant[idx]  <- n_inf_infant[idx]  + 1L
        if (ac == 2L) n_inf_sibling[idx] <- n_inf_sibling[idx] + 1L
        if (ac == 3L) n_inf_adult[idx]   <- n_inf_adult[idx]   + 1L
        if (ac == 4L) n_inf_elder[idx]   <- n_inf_elder[idx]   + 1L
      }
    }

    #Generate person-day rows for non-index susceptibles
    for (i in hhdat[is_index == FALSE]$ID_indiv) {
      rec <- hhdat[ID_indiv == i]
      inf_d <- rec$inf_day_rl

      #Self infectious span
      idx_self <- safe_span_idx(rec$infectious_day_rl, rec$infectious_end_day_rl, tmax)
      self_inf <- integer(tmax + 1L); if (length(idx_self)) self_inf[idx_self] <- 1L

      #Age info (1,2,3,4); dummies with infant as baseline
      ac <- as.integer(rec$age_cat[1L])
      agegrp2 <- as.integer(ac == 2L)
      agegrp3 <- as.integer(ac == 3L)
      agegrp4 <- as.integer(ac == 4L)

      #Self-by-role flags (exclude own infectiousness from HH counts)
      self_infant  <- as.integer(ac == 1L) * self_inf
      self_sibling <- as.integer(ac == 2L) * self_inf
      self_adult   <- as.integer(ac == 3L) * self_inf
      self_elder   <- as.integer(ac == 4L) * self_inf

      for (d in 0:tmax) {
        if (!is.na(inf_d) && d > inf_d) break
        idx <- d + 1L

        #Subtract self only on days when self is infectious
        n_tot_d <- n_inf_total[idx]   - self_inf[idx]
        n_inf_d <- n_inf_infant[idx]  - self_infant[idx]
        n_sib_d <- n_inf_sibling[idx] - self_sibling[idx]
        n_ad_d  <- n_inf_adult[idx]   - self_adult[idx]
        n_el_d  <- n_inf_elder[idx]   - self_elder[idx]

        #Defensive clamp
        if (n_tot_d < 0L || n_inf_d < 0L || n_sib_d < 0L || n_ad_d < 0L || n_el_d < 0L) {
          n_tot_d <- max(0L, n_tot_d); n_inf_d <- max(0L, n_inf_d)
          n_sib_d <- max(0L, n_sib_d); n_ad_d  <- max(0L, n_ad_d); n_el_d <- max(0L, n_el_d)
        }

        # Build per-day covariate list: day-specific if list-column vector, else scalar
        cov_list <- list()
        if (length(covariate_cols)) {
          for (col in covariate_cols) {
            cov_list[[col]] <- cov_value_for_day(rec, col, d, default = NA)
          }
        }

        rows[[length(rows) + 1L]] <- c(list(
          #Age dummies (infant is reference)
          agegrp2   = agegrp2,
          agegrp3   = agegrp3,
          agegrp4   = agegrp4,

          #Household infectious counts (excluding self)
          n_inf         = n_tot_d,
          n_inf_infant  = n_inf_d,
          n_inf_sibling = n_sib_d,
          n_inf_adult   = n_ad_d,
          n_inf_elder   = n_el_d,

          cases     = cases_t[idx],
          event     = as.integer(!is.na(inf_d) && d == inf_d),

          ID_indiv  = i,
          ID_hh     = hh,
          day       = d
        ), cov_list)
      }
    }
  }
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}




#' Run penalized maximum likelihood parameter estimation
#'
#' Repeats optimization \code{n_runs} times from jittered starts to estimate
#' community and household transmission parameters with optional covariate effects
#' and quadratic penalties.
#'
#' @param long_dt data.table or data.frame. Person day data from \code{\link{build_person_day_table}}.
#' @param n_runs Integer. Number of optimization runs for multi start.
#' @param start_par Numeric vector or \code{NULL}. Initial parameters. If \code{NULL},
#'   a vector of the correct length is initialized with \code{delta0_true} and \code{alpha0_true}
#'   in the intercept slots and zeros elsewhere.
#' @param lambda Numeric. Base L2 penalty on \code{gamma} age effects and \code{z_*} role offsets.
#' @param lambda0,lambda_alpha Numeric. Penalty strengths that pull \code{delta0} and \code{alpha0}
#'   toward \code{delta0_true} and \code{alpha0_true}.
#' @param delta0_true,alpha0_true Numeric anchors for the intercept penalties on the logit scale.
#' @param comm_covariate_cols Character vector of community risk covariate names with no intercept.
#' @param hh_covariate_cols Character vector of household covariates shared across roles.
#' @param hh_by_role Logical. If \code{TRUE}, allow role specific household covariates through
#'   \code{hh_role_covariate_cols}.
#' @param hh_role_covariate_cols Named list with elements \code{infant}, \code{sibling}, \code{adult}, \code{elder}
#'   that give covariate names per role. Falls back to \code{hh_covariate_cols} if a role list is missing.
#' @param standardize_covariates Logical. Z score non binary numeric columns in model matrices.
#' @param lambda_comm,lambda_hh Numeric. L2 penalties for community and household covariate coefficients.
#' @param verbose Logical. If \code{TRUE}, print notes on dropped or unknown covariates and initialization.
#'
#' @return A numeric matrix of dimension \code{n_runs} by \code{n_parameters} with column names that match
#'   the parameter layout, for example \code{delta0}, \code{gamma2}, \code{alpha0}, \code{z_sib}, and \code{theta_*}.
running_parameter_estimation <- function(long_dt, n_runs, start_par = NULL,
                                         lambda,
                                         lambda0,
                                         lambda_alpha,
                                         delta0_true,
                                         alpha0_true,

                                         comm_covariate_cols = NULL,
                                         hh_covariate_cols   = NULL,
                                         hh_by_role = FALSE,
                                         hh_role_covariate_cols = NULL,
                                         standardize_covariates = FALSE,

                                         lambda_comm = lambda,
                                         lambda_hh   = lambda,
                                         verbose = TRUE) {

  dat <- data.table::as.data.table(long_dt)

  log1mexp <- function(logp) ifelse(logp < log(0.5), log1p(-exp(logp)), log(-expm1(logp)))
  mmprod   <- function(X, b) if (is.null(X) || NCOL(X) == 0L || is.null(b)) 0 else drop(X %*% b)
  `%||%`   <- function(a,b) if (is.null(a)) b else a

  ##########################################################################
  #                                                                        #
  #                 Building covariate matrix X_comm                       #
  #                                                                        #
  ##########################################################################
  build_X <- function(D, cols, standardize = FALSE, label = "X", verbose = TRUE) {
    req <- cols %||% character(0)
    keep <- intersect(req, names(D))
    drop <- setdiff(req, keep)
    if (verbose && length(drop))
      warning(sprintf("%s: dropping unknown covariates: %s",
                      label, paste(drop, collapse = ", ")), call. = FALSE)
    if (length(keep) == 0L) return(NULL)
    M <- data.frame(D[, ..keep], check.names = FALSE)
    for (nm in names(M)) if (is.character(M[[nm]])) M[[nm]] <- factor(M[[nm]])
    X <- model.matrix(~ 0 + ., data = M)  # no intercept
    if (standardize && NCOL(X) > 0L) {
      for (j in seq_len(NCOL(X))) {
        uj <- sort(unique(X[, j]))
        is_binary <- length(uj) <= 2 && all(uj %in% c(0,1))
        if (!is_binary) {
          m <- mean(X[, j]); s <- sd(X[, j])
          if (is.finite(s) && s > 0) X[, j] <- (X[, j] - m)/s
        }
      }
    }
    X
  }

  need <- c("agegrp2","agegrp3","agegrp4",
            "n_inf_infant","n_inf_sibling","n_inf_adult","n_inf_elder","event")
  miss <- setdiff(need, names(dat))
  if (length(miss)) stop("long_dt missing columns: ", paste(miss, collapse=", "))

  age_mat <- as.matrix(dat[, .(agegrp2, agegrp3, agegrp4)])
  X_comm  <- build_X(dat, comm_covariate_cols, standardize_covariates,
                     label = "Community", verbose = verbose)

  if (!isTRUE(hh_by_role)) {
    X_hh <- build_X(dat, hh_covariate_cols, standardize_covariates,
                    label = "Household", verbose = verbose)
    X_hh_inf <- X_hh_sib <- X_hh_ad <- X_hh_el <- NULL
  } else {
    role_cols <- list(
      infant  = (hh_role_covariate_cols$infant  %||% hh_covariate_cols),
      sibling = (hh_role_covariate_cols$sibling %||% hh_covariate_cols),
      adult   = (hh_role_covariate_cols$adult   %||% hh_covariate_cols),
      elder   = (hh_role_covariate_cols$elder   %||% hh_covariate_cols)
    )
    X_hh_inf <- build_X(dat, role_cols$infant,  standardize_covariates,
                        label = "HH-infant",  verbose = verbose)
    X_hh_sib <- build_X(dat, role_cols$sibling, standardize_covariates,
                        label = "HH-sibling", verbose = verbose)
    X_hh_ad  <- build_X(dat, role_cols$adult,   standardize_covariates,
                        label = "HH-adult",    verbose = verbose)
    X_hh_el  <- build_X(dat, role_cols$elder,   standardize_covariates,
                        label = "HH-elder",    verbose = verbose)
    X_hh <- NULL
  }

  # Parameter layout and names
  n_comm <- if (is.null(X_comm)) 0L else NCOL(X_comm)

  ##########################################################################
  #                                                                        #
  #                  Building covariate matrix X_hh                        #
  #                                                                        #
  ##########################################################################
  if (!isTRUE(hh_by_role)) {
    n_hh <- if (is.null(X_hh)) 0L else NCOL(X_hh)
    n_par <- 8 + n_comm + n_hh
    par_names <- c("delta0","gamma2","gamma3","gamma4",
                   "alpha0","z_sib","z_ad","z_el",
                   if (n_comm > 0L) paste0("theta_comm_", colnames(X_comm)) else NULL,
                   if (n_hh   > 0L) paste0("theta_hh_",   colnames(X_hh))   else NULL)
  } else {
    n_inf <- if (is.null(X_hh_inf)) 0L else NCOL(X_hh_inf)
    n_sib <- if (is.null(X_hh_sib)) 0L else NCOL(X_hh_sib)
    n_ad  <- if (is.null(X_hh_ad))  0L else NCOL(X_hh_ad)
    n_el  <- if (is.null(X_hh_el))  0L else NCOL(X_hh_el)
    n_par <- 8 + n_comm + (n_inf + n_sib + n_ad + n_el)
    par_names <- c("delta0","gamma2","gamma3","gamma4",
                   "alpha0","z_sib","z_ad","z_el",
                   if (n_comm > 0L) paste0("theta_comm_",   colnames(X_comm))   else NULL,
                   if (n_inf  > 0L) paste0("theta_hh_inf_", colnames(X_hh_inf)) else NULL,
                   if (n_sib  > 0L) paste0("theta_hh_sib_", colnames(X_hh_sib)) else NULL,
                   if (n_ad   > 0L) paste0("theta_hh_ad_",  colnames(X_hh_ad))  else NULL,
                   if (n_el   > 0L) paste0("theta_hh_el_",  colnames(X_hh_el))  else NULL)
  }

  if (is.null(start_par) || length(start_par) != n_par) {
    start_par <- numeric(n_par)
    start_par[1] <- delta0_true  # delta0
    start_par[5] <- alpha0_true  # alpha0
    if (verbose)
      message("Initialized start_par of length ", n_par, " (based on available covariates).")
  }

  # Negative log-likelihood
  negll <- function(par, dat, eps = 1e-12,
                    lambda = 0.01, lambda0 = 0.2, lambda_alpha = 5,
                    delta0_true = qlogis(0.002),
                    alpha0_true = qlogis(0.2),
                    lambda_comm = lambda, lambda_hh = lambda) {

    delta0 <- par[1]
    gamma  <- par[2:4]
    alpha0 <- par[5]
    z_sib  <- par[6]; z_ad <- par[7]; z_el <- par[8]

    k <- 8
    theta_comm <- if (n_comm > 0L) { th <- par[(k+1):(k+n_comm)]; k <- k + n_comm; th } else NULL

    if (!isTRUE(hh_by_role)) {
      n_hh <- if (is.null(X_hh)) 0L else NCOL(X_hh)
      theta_hh <- if (n_hh > 0L) { th <- par[(k+1):(k+n_hh)]; k <- k + n_hh; th } else NULL

      eta_comm <- delta0 + drop(age_mat %*% gamma) + mmprod(X_comm, theta_comm)
      eta_inf  <- alpha0 + mmprod(X_hh, theta_hh)
      eta_sib  <- alpha0 + z_sib + mmprod(X_hh, theta_hh)
      eta_ad   <- alpha0 + z_ad  + mmprod(X_hh, theta_hh)
      eta_el   <- alpha0 + z_el  + mmprod(X_hh, theta_hh)
    } else {
      n_inf <- if (is.null(X_hh_inf)) 0L else NCOL(X_hh_inf)
      n_sib <- if (is.null(X_hh_sib)) 0L else NCOL(X_hh_sib)
      n_ad  <- if (is.null(X_hh_ad))  0L else NCOL(X_hh_ad)
      n_el  <- if (is.null(X_hh_el))  0L else NCOL(X_hh_el)

      theta_inf <- if (n_inf > 0L) { th <- par[(k+1):(k+n_inf)]; k <- k + n_inf; th } else NULL
      theta_sib <- if (n_sib > 0L) { th <- par[(k+1):(k+n_sib)]; k <- k + n_sib; th } else NULL
      theta_ad  <- if (n_ad  > 0L) { th <- par[(k+1):(k+n_ad )]; k <- k + n_ad ; th } else NULL
      theta_el  <- if (n_el  > 0L) { th <- par[(k+1):(k+n_el )]; k <- k + n_el ; th } else NULL

      eta_comm <- delta0 + drop(age_mat %*% gamma) + mmprod(X_comm, theta_comm)
      eta_inf  <- alpha0 + mmprod(X_hh_inf, theta_inf)
      eta_sib  <- alpha0 + z_sib + mmprod(X_hh_sib, theta_sib)
      eta_ad   <- alpha0 + z_ad  + mmprod(X_hh_ad,  theta_ad)
      eta_el   <- alpha0 + z_el  + mmprod(X_hh_el,  theta_el)
    }

    # Probabilities
    p_comm_prob <- pmin(pmax(plogis(eta_comm), eps), 1 - eps)
    p_hh_infant  <- pmin(pmax(plogis(eta_inf), eps), 1 - eps)
    p_hh_sibling <- pmin(pmax(plogis(eta_sib), eps), 1 - eps)
    p_hh_adult   <- pmin(pmax(plogis(eta_ad), eps), 1 - eps)
    p_hh_elder   <- pmin(pmax(plogis(eta_el), eps), 1 - eps)

    # Escape probability and total infection probability
    log_escape <- log1mexp(log(p_comm_prob)) +
      dat$n_inf_infant  * log1mexp(log(p_hh_infant))  +
      dat$n_inf_sibling * log1mexp(log(p_hh_sibling)) +
      dat$n_inf_adult   * log1mexp(log(p_hh_adult))   +
      dat$n_inf_elder   * log1mexp(log(p_hh_elder))

    p_tot <- pmin(pmax(-expm1(log_escape), eps), 1 - eps)
    nll <- -sum(dat$event * log(p_tot) + (1 - dat$event) * log1mexp(log(p_tot)))

    # Penalties
    pen <- 0
    pen <- pen + lambda0      * (delta0 - delta0_true)^2
    pen <- pen + lambda_alpha * (alpha0 - alpha0_true)^2
    pen <- pen + lambda       * (sum(gamma^2) + z_sib^2 + z_ad^2 + z_el^2)

    if (!is.null(theta_comm)) pen <- pen + lambda_comm * sum(theta_comm^2)
    if (!isTRUE(hh_by_role)) {
      if (!is.null(X_hh))     pen <- pen + lambda_hh   * sum(theta_hh^2)
    } else {
      if (!is.null(X_hh_inf)) pen <- pen + lambda_hh   * sum(theta_inf^2)
      if (!is.null(X_hh_sib)) pen <- pen + lambda_hh   * sum(theta_sib^2)
      if (!is.null(X_hh_ad))  pen <- pen + lambda_hh   * sum(theta_ad^2)
      if (!is.null(X_hh_el))  pen <- pen + lambda_hh   * sum(theta_el^2)
    }

    nll + pen
  }

  # Multi-start optimization
  multi_start_optim <- function(start_par, fn, dat, n_start = 5, ...) {
    best_fit <- NULL; best_val <- Inf
    for (i in seq_len(n_start)) {
      trial <- start_par + rnorm(length(start_par), 0, 1)
      fit <- optim(trial, fn = fn, dat = dat, ..., method = "BFGS",
                   control = list(maxit = 2e4))
      if (fit$value < best_val) { best_val <- fit$value; best_fit <- fit }
    }
    best_fit
  }

  theta_mat <- matrix(NA_real_, n_runs, length(start_par))
  colnames(theta_mat) <- par_names

  for (m in seq_len(n_runs)) {
    fit <- multi_start_optim(
      start_par, fn = negll, dat = dat, n_start = 5,
      lambda = lambda, lambda0 = lambda0, lambda_alpha = lambda_alpha,
      delta0_true = delta0_true, alpha0_true = alpha0_true,
      lambda_comm = lambda_comm, lambda_hh = lambda_hh
    )
    theta_mat[m, ] <- fit$par
  }
  theta_mat
}



#' Generate standardized synthetic data for one household
#'
#' Simulates test-day observations for one household over a date window with
#' community seasonality, within-household transmission, adaptive testing,
#' baseline/partial immunity, and optional covariates.
#'
#' @param household_id Integer; household identifier written to \code{HH}.
#' @param hh.size Integer; household size.
#' @param tests.per.week Integer; tests per person per week.
#' @param p.comm.base.infant.fix Numeric; baseline community infection prob/day (infant).
#' @param p.comm.multiplier.sibling,p.comm.multiplier.parent,p.comm.multiplier.elder
#'   Numeric; community multipliers by role.
#' @param p.hh.base.infant Numeric; baseline within-household infection prob/day (infant source).
#' @param p.hh.multiplier.sibling,p.hh.multiplier.parent,p.hh.multiplier.elder
#'   Numeric; within-household multipliers by source role.
#' @param p.imm.base.sibling,p.imm.base.parent,p.imm.base.elder Numeric; baseline immunity at day 1.
#' @param partial.immunity.infant,partial.immunity.sibling,partial.immunity.parent,partial.immunity.elder
#'   Numeric; partial-immunity modifiers by role.
#' @param duration.latent Integer; mean latent period (days).
#' @param duration.infect.inf Integer; mean infectious duration for infants (days).
#' @param multiplier.dur.sibpar Numeric; infectious-duration multiplier for non-infants.
#' @param p.detect Numeric; detection probability if infected on a test day.
#' @param amplitude,phase Numeric; seasonality parameters for community risk.
#' @param start_date,end_date \code{Date}; simulation window.
#' @param Covariates Logical; generate additional covariates.
#' @param Covariates_list Character; covariate names.
#' @param Covariate_specifications List; optional per-covariate specs.
#'
#' @return Data frame with columns \code{HH}, \code{individual_ID}, \code{role},
#'   \code{test_date} (1 = \code{start_date}), \code{infection_status},
#'   \code{community_risk}, and optional covariates.
generate_synthetic_data_one <- function(
    household_id,
    hh.size = sample(3:7,1),
    tests.per.week = 2,

    p.comm.base.infant.fix = 0.001,
    p.comm.multiplier.sibling = 1,
    p.comm.multiplier.parent = 1,
    p.comm.multiplier.elder = 1,

    p.hh.base.infant = 0.1,
    p.hh.multiplier.sibling = 1,
    p.hh.multiplier.parent = 1,
    p.hh.multiplier.elder = 1,

    p.imm.base.sibling = 1e-10,
    p.imm.base.parent = 1e-10,
    p.imm.base.elder = 1e-10,

    partial.immunity.infant = 1e-10,
    partial.immunity.sibling = 1e-10,
    partial.immunity.parent = 1e-10,
    partial.immunity.elder = 1e-10,

    duration.latent = 2,
    duration.infect.inf = 3,
    multiplier.dur.sibpar = 0.5,
    p.detect = 0.999,

    amplitude = 0,
    phase = -0.408,
    start_date = as.Date("2024-09-21"),
    end_date = as.Date("2025-04-17"),

    #Generate synthetic covariates
    Covariates = FALSE,
    Covariates_list = c("Vaccination status", "Antibody Level"),
    Covariate_specifications = NULL
)
{

  ##########################################################################
  #                                                                        #
  #              Helper function for cleaning covariate names              #
  #                                                                        #
  ##########################################################################
  sanitize_name <- function(x) {
    x2 <- tolower(trimws(x)); x2 <- gsub("[^a-z0-9]+", "_", x2); x2 <- gsub("^_|_$","",x2)
    if (x2 == "") x2 <- "covariate"; x2
  }

  ##########################################################################
  #                                                                        #
  #                    Distributions covariates may take                   #
  #                                                                        #
  ##########################################################################
  #AR(1)
  ar1_normal <- function(T, mu=0, sigma=1, rho=0.9) {
    x <- numeric(T); x[1] <- rnorm(1, mu, sigma)
    for (t in 2:T) x[t] <- mu + rho*(x[t-1]-mu) + rnorm(1, 0, sigma*sqrt(1-rho^2))
    x
  }

  #AR(1) on log-scale
  ar1_lognormal <- function(T, mu_log=log(100), sigma_log=0.5, rho=0.9) {
    z <- ar1_normal(T, mu=mu_log, sigma=sigma_log, rho=rho)
    exp(z)
  }

  #Random walk (Normal steps), with optional bounds
  clip <- function(x, lo=-Inf, hi=Inf) pmin(pmax(x, lo), hi)
  rw_normal <- function(T, mu0=0, sigma_step=0.05, lo=-Inf, hi=Inf) {
    x <- numeric(T); x[1] <- mu0
    for (t in 2:T) x[t] <- clip(x[t-1] + rnorm(1, 0, sigma_step), lo, hi)
    x
  }

  #Exponential decay on log scale with noise
  decay_lognormal <- function(T, x0=100, lambda=0.01, noise_sdlog=0.2, floor=1) {
    t <- 0:(T-1)
    mean_traj <- pmax(floor, x0 * exp(-lambda * t))
    noise <- rlnorm(T, 0, noise_sdlog)
    mean_traj * noise
  }


  # -------- timeline & testing -------------------------------------------
  time.steps <- as.integer(end_date - start_date) + 1

  test.days1 <- seq(1, time.steps, by=7)
  test.days2 <- seq(5, time.steps, by=7)
  test.days3 <- seq(3, time.steps, by=7)
  if (tests.per.week == 1)      test.days <- test.days1
  else if (tests.per.week == 2) test.days <- c(test.days1, test.days2)
  else if (tests.per.week == 3) test.days <- c(test.days1, test.days2, test.days3)
  else                          test.days <- test.days1
  test.days <- sort(test.days)
  baseline.test.days <- test.days

  # -------- household composition ----------------------------------------
  n.adult.base <- 2
  remaining <- hh.size - (1 + n.adult.base)
  n.elder <- ifelse(remaining > 0, sample(0:min(2, remaining), 1), 0)
  remaining <- remaining - n.elder
  n.sib <- max(remaining, 0)

  hh.roles <- c("infant", rep("adult", n.adult.base), rep("sibling", n.sib), rep("elder", n.elder))

  latent <- array(0, dim=c(4, time.steps, hh.size))
  infectious <- array(0, dim=c(4, time.steps, hh.size))
  immune <- matrix(NA, nrow=time.steps, ncol=hh.size)

  immune[1, hh.roles=="infant"] <- 0
  immune[1, hh.roles %in% "adult"] <- rbinom(sum(hh.roles=="adult"), 1, p.imm.base.parent)
  immune[1, hh.roles %in% "elder"] <- rbinom(sum(hh.roles=="elder"), 1, p.imm.base.elder)
  immune[1, hh.roles=="sibling"]   <- rbinom(sum(hh.roles=="sibling"), 1, p.imm.base.sibling)

  p.comm.init <- p.comm.base.infant.fix * exp(amplitude * cos((2*pi*(1+40*7)/365.25) + phase))
  p.comm.init.vec <- ifelse(hh.roles=="infant", p.comm.init,
                            ifelse(hh.roles=="adult", p.comm.init * p.comm.multiplier.parent,
                                   ifelse(hh.roles=="elder", p.comm.init * p.comm.multiplier.elder,
                                          p.comm.init * p.comm.multiplier.sibling)))
  latent[1,1,] <- ifelse(immune[1,] > 0, 0, rbinom(hh.size, 1, p.comm.init.vec))
  infectious[1,1,] <- 0

  p_comm_mat <- matrix(NA_real_, nrow=time.steps, ncol=hh.size)
  p_comm_mat[1, ] <- p.comm.init.vec

  for (i in 2:time.steps) {
    p.comm.base.infant <- p.comm.base.infant.fix * exp(amplitude * cos((2*pi*(i+40*7)/365.25) + phase))
    for (j in 1:hh.size) {
      if (hh.roles[j] == "infant") {
        p.comm <- p.comm.base.infant; partial.immunity <- partial.immunity.infant; duration.infect <- duration.infect.inf
      } else if (hh.roles[j] == "adult") {
        p.comm <- p.comm.base.infant * p.comm.multiplier.parent; partial.immunity <- partial.immunity.parent; duration.infect <- duration.infect.inf * multiplier.dur.sibpar
      } else if (hh.roles[j] == "elder") {
        p.comm <- p.comm.base.infant * p.comm.multiplier.elder; partial.immunity <- partial.immunity.elder; duration.infect <- duration.infect.inf * multiplier.dur.sibpar
      } else {
        p.comm <- p.comm.base.infant * p.comm.multiplier.sibling; partial.immunity <- partial.immunity.sibling; duration.infect <- duration.infect.inf * multiplier.dur.sibpar
      }
      p_comm_mat[i, j] <- p.comm

      N.inf.infant  <- sum(infectious[,(i-1),hh.roles=="infant"])
      N.inf.sibling <- sum(infectious[,(i-1),hh.roles=="sibling"])
      N.inf.parent  <- sum(infectious[,(i-1),hh.roles=="adult"])
      N.inf.elder   <- sum(infectious[,(i-1),hh.roles=="elder"])

      log.qi <- log(1 - p.comm) +
        N.inf.infant  * log(1 - p.hh.base.infant) +
        N.inf.sibling * log(1 - p.hh.base.infant * p.hh.multiplier.sibling) +
        N.inf.parent  * log(1 - p.hh.base.infant * p.hh.multiplier.parent) +
        N.inf.elder   * log(1 - p.hh.base.infant * p.hh.multiplier.elder)

      qi  <- exp(log.qi); pri <- 1 - qi

      was_inf <- as.integer(max(infectious[,(i-1),j]) > 0)
      was_lat <- as.integer(max(latent[,(i-1),j]) > 0)
      susceptible <- (1 - was_inf) * (1 - was_lat) * (1 - as.integer(immune[(i-1),j] > 0))

      new.inf <- rbinom(1,1,pri) * susceptible +
        rbinom(1,1,pri*partial.immunity) * immune[(i-1),j] *
        as.numeric(((i-1) - max(c(0, which(immune[1:(i-1), j] == 0)))) > 25 | ((i-1) < 26))

      stay.latent     <- max(latent[,(i-1),j])    * rbinom(1,1,(1 - 1/duration.latent)^(1/4))
      stay.infectious <- max(infectious[,(i-1),j]) * rbinom(1,1,(1 - 1/duration.infect)^(1/4))

      latent[1,i,j] <- new.inf + latent[1,(i-1),j] * stay.latent
      latent[2,i,j] <- latent[1,(i-1),j]*(1-stay.latent) + latent[2,(i-1),j]*stay.latent
      latent[3,i,j] <- latent[2,(i-1),j]*(1-stay.latent) + latent[3,(i-1),j]*stay.latent
      latent[4,i,j] <- latent[3,(i-1),j]*(1-stay.latent) + latent[4,(i-1),j]*stay.latent

      infectious[1,i,j] <- latent[4,(i-1),j]*(1-stay.latent) + infectious[1,(i-1),j]*stay.infectious
      infectious[2,i,j] <- infectious[1,(i-1),j]*(1-stay.infectious) + infectious[2,(i-1),j]*stay.infectious
      infectious[3,i,j] <- infectious[2,(i-1),j]*(1-stay.infectious) + infectious[3,(i-1),j]*stay.infectious
      infectious[4,i,j] <- infectious[3,(i-1),j]*(1-stay.infectious) + infectious[4,(i-1),j]*stay.infectious

      immune[i,j] <- infectious[4,(i-1),j]*(1-stay.infectious) + immune[i-1,j] - immune[i-1,j]*new.inf
    }
  }

  latent2 <- apply(latent, c(2,3), sum)
  infectious2 <- apply(infectious, c(2,3), sum)
  infected <- (latent2 + infectious2) > 0

  test.schedule <- rep(FALSE, time.steps); test.schedule[baseline.test.days] <- TRUE
  detect.inf_full <- matrix(0, nrow=time.steps, ncol=hh.size)
  daily.testing <- FALSE; neg.streak <- rep(0, hh.size)

  for (d in 1:time.steps) {
    if (!test.schedule[d]) next
    detect.inf_full[d,] <- rbinom(hh.size, 1, infected[d,] * p.detect)
    if (!daily.testing && any(detect.inf_full[d,] > 0)) {
      daily.testing <- TRUE
      neg.streak <- ifelse(detect.inf_full[d,]==0, 1, 0)
      if (d < time.steps) test.schedule[(d+1):time.steps] <- TRUE
    } else if (daily.testing) {
      neg.streak <- neg.streak + (detect.inf_full[d,]==0)
      neg.streak[detect.inf_full[d,]>0] <- 0
      if (all(neg.streak >= 2)) {
        daily.testing <- FALSE
        if (d < time.steps) test.schedule[(d+1):time.steps] <- FALSE
        remaining <- baseline.test.days[baseline.test.days > d]
        test.schedule[remaining] <- TRUE
      }
    }
  }

  test.days <- which(test.schedule)
  test_df <- data.frame(
    HH = household_id,
    individual_ID = rep(1:hh.size, each = length(test.days)),
    role = rep(hh.roles, each = length(test.days)),
    test_date = rep(test.days, times = hh.size),
    infection_status = as.integer(as.vector(infected[test.days, ])),
    community_risk = as.vector(p_comm_mat[test.days, ])
  )

  ##########################################################################
  #                                                                        #
  #                   When Covariates parameter is TRUE                    #
  #                                                                        #
  ##########################################################################
  if (isTRUE(Covariates)) {

    #Normalize requested names
    req <- unique(tolower(trimws(Covariates_list)))
    req <- req[req != ""]

    #Default specification if none is given for vaccination status
    default_specifications <- list(
      vaccination_status = list(
        type="binary", dist="bernoulli", time_varying=FALSE,
        params=list(probability_by_role=function(role){
          ifelse(role=="infant", 0.05,
                 ifelse(role=="sibling",0.30,
                        ifelse(role=="adult",  0.70, 0.85)))
        })
      ),

      #Default specification if none is given for antibody level
      antibody_level = list(
        type="continuous", dist="ar1-lognormal", time_varying=TRUE,
        params=list(mu_log=log(100), sigma_log=0.5, rho=0.9)
      )
    )

    #If user passed specifications, merge on top of defaults
    specifications <- default_specifications
    if (!is.null(Covariate_specifications) && length(Covariate_specifications)>0) {
      # user specifications can add new or override existing
      for (nm in names(Covariate_specifications)) specifications[[sanitize_name(nm)]] <- Covariate_specifications[[nm]]
    }

    #If user did not pass specifications, will consider covariate comes from a Normal distribution
    for (nm in req) {
      key <- sanitize_name(nm)
      if (is.null(specifications[[key]])) {
        specifications[[key]] <- list(type="continuous", dist="normal", time_varying=FALSE,
                                      params=list(mean=0, sd=1))
      }
    }

    #Generator per specification returns matrix [time.steps x hh.size]
    gen_series <- function(spec) {
      type <- spec$type; dist <- spec$dist; tv <- isTRUE(spec$time_varying)
      par <- spec$params %||% list()
      out <- matrix(NA_real_, nrow=time.steps, ncol=hh.size)

      for (j in 1:hh.size) {
        role <- hh.roles[j]

        #Binary covariates
        if (identical(type,"binary")) {
          if (!tv || is.null(dist) || dist=="bernoulli") {
            get_p_for_role <- function(par, role) {
              if (!is.null(par$probability_by_role)) {
                v <- par$probability_by_role
                if (is.function(v)) {
                  return(v(role))
                } else if (is.numeric(v)) {
                  if (!is.null(names(v)) && role %in% names(v)) return(unname(v[role]))  # Named vector
                  if (length(v) == 1L) return(v)                                       # Single number
                  stop("Probability_by_role must be a function, a single number, or a named vector with roles.")
                } else {
                  stop("Unsupported probability_by_role type. Use function, single number, or named vector.")
                }
              }
              return(par$p %||% 0.2)
            }
            p <- get_p_for_role(par, role)
            x <- rep(rbinom(1,1,p), time.steps)
          } else stop("Unsupported binary dist: ", dist)
        } else if (identical(type,"categorical")) {
          levels <- par$levels %||% c("A","B","C")
          probs  <- par$probs  %||% rep(1/length(levels), length(levels))
          if (!tv) {
            val <- sample(levels, 1, prob=probs, replace=TRUE)
            x <- rep(val, time.steps)
          } else {
            x <- sample(levels, time.steps, prob=probs, replace=TRUE)
          }

          #Continuous covariates
        } else if (identical(type,"continuous")) {
          if (!tv) {
            if (dist=="normal") {
              mu <- par$mean %||% 0; sd <- par$sd %||% 1; x <- rep(rnorm(1,mu,sd), time.steps)
            } else if (dist=="lognormal") {
              mu <- par$meanlog %||% 0; sd <- par$sdlog %||% 1; x <- rep(rlnorm(1,mu,sd), time.steps)
            } else if (dist=="uniform") {
              a <- par$min %||% 0; b <- par$max %||% 1; x <- rep(runif(1,a,b), time.steps)
            } else stop("Unsupported continuous dist (static): ", dist)
          } else { # time-varying
            if (dist=="ar1-normal") {
              mu <- par$mu %||% 0; sd <- par$sd %||% 1; rho <- par$rho %||% 0.9
              x <- ar1_normal(time.steps, mu, sd, rho)
            } else if (dist=="ar1-lognormal") {
              mu <- par$mu_log %||% log(100); sd <- par$sigma_log %||% 0.5; rho <- par$rho %||% 0.9
              x <- ar1_lognormal(time.steps, mu, sd, rho)
            } else if (dist=="rw-normal") {
              mu0 <- par$mu0 %||% 0; s <- par$sigma_step %||% 0.05
              lo <- par$min %||% -Inf; hi <- par$max %||% Inf
              x <- rw_normal(time.steps, mu0, s, lo, hi)
            } else if (dist=="decay-lognormal") {
              x0 <- par$x0 %||% 100; lam <- par$lambda %||% 0.01
              s  <- par$noise_sdlog %||% 0.2; fl <- par$floor %||% 1
              x <- decay_lognormal(time.steps, x0, lam, s, fl)
            } else stop("Unsupported continuous dist (time-varying): ", dist)
          }
        } else stop("Unsupported type: ", type)

        out[,j] <- x
      }
      out
    }

    `%||%` <- function(a,b) if (is.null(a)) b else a

    #Build only for the requested covariates (or all if Covariates_list empty)
    wanted_keys <- if (length(req)>0) unique(sapply(req, sanitize_name)) else names(specifications)

    for (key in wanted_keys) {
      ser <- gen_series(specifications[[key]])  # [time.steps x hh.size]
      col_name <- key
      vals <- as.vector(ser[test.days, ])  # Subset to observation days
      test_df[[col_name]] <- vals
    }
  }
  return(test_df)
}
