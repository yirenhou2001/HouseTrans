#' Summarize infection episodes (one row per individual)
#'
#' Collapses long testing records to episode-level features and optional covariate
#' summaries. Supply either a single data frame \code{df} with the six core columns,
#' or provide all six vectors individually. An optional covariate frame can be merged.
#'
#' @param df Data frame with columns \code{HH}, \code{individual_ID}, \code{role},
#'   \code{test_date}, \code{infection_status}, \code{community_risk}. Ignored when all six
#'   vector inputs are provided.
#' @param Household_ID,Individual_ID,Household_role,Sample_test_days,Infectious_status,Community_rate_infection Vectors used only when \code{df} is \code{NULL}.
#' @param Covariate_DataFrame Optional data frame to merge before summarizing; joined by available keys among
#'   \code{c("HH","individual_ID","test_date")} (at least \code{individual_ID} required).
#' @param covariate_cols Optional character vector restricting which covariates from
#'   \code{Covariate_DataFrame} are summarized.
#'
#' @return A data frame with one row per individual including:
#'   \code{HH}, \code{individual_ID}, episode counts and timing fields
#'   (e.g., \code{infection.detected.start}, \code{infection.detected.end},
#'   \code{infection.true.duration}, \code{last_negative}),
#'   \code{infection.infectious.day} (comma-separated days),
#'   \code{community.risk}, \code{role}, and (if present) basic summaries for covariates.
#'
#' @details Records are ordered by \code{HH}, \code{individual_ID}, \code{test_date}.
#' Episodes are runs of \code{infection_status == 1}; first-episode fields refer to the earliest run.
#'
#' @seealso \code{\link{summarize_individuals}}, \code{\link{main_parameter_estimation_pipeline}}
data_summarization <- function(
    df = NULL,
    Household_ID = NULL, Individual_ID = NULL, Household_role = NULL,
    Sample_test_days = NULL, Infectious_status = NULL, Community_rate_infection = NULL,

    # User covariate data frame
    Covariate_DataFrame = NULL,
    #Allow user to specify which covariates to summarize
    covariate_cols = NULL
) {
  tryCatch({

    ##########################################################################
    #                                                                        #
    #                    Build and check data frame                          #
    #                                                                        #
    ##########################################################################
    core_cols <- c("HH","individual_ID","role","test_date","infection_status","community_risk")

    if (!is.null(Household_ID) && !is.null(Individual_ID) && !is.null(Household_role) &&
        !is.null(Sample_test_days) && !is.null(Infectious_status) && !is.null(Community_rate_infection)) {
      newdf <- data.frame(
        HH = Household_ID,
        individual_ID = Individual_ID,
        role = Household_role,
        test_date = Sample_test_days,
        infection_status = Infectious_status,
        community_risk = Community_rate_infection,
        stringsAsFactors = FALSE
      )
      df <- newdf[order(newdf$HH, newdf$individual_ID, newdf$test_date), ]
    } else if (!is.null(df)) {
      missing_core <- setdiff(core_cols, names(df))
      if (length(missing_core))
        stop("`df` is missing required columns: ", paste(missing_core, collapse=", "))
      df <- df[order(df$HH, df$individual_ID, df$test_date), ]
    } else {
      stop("Either a full data frame `df` or all six vector inputs must be provided.")
    }

    ##########################################################################
    #                                                                        #
    #                   Merge user covariate data frame                      #
    #                                                                        #
    ##########################################################################
    if (!is.null(Covariate_DataFrame)) {
      cv <- Covariate_DataFrame
      # Choose best join keys available
      join_keys <- intersect(c("HH","individual_ID","test_date"), names(cv))
      if (all(c("HH","individual_ID","test_date") %in% join_keys)) {
        keys <- c("HH","individual_ID","test_date")
      } else if (all(c("HH","individual_ID") %in% join_keys)) {
        keys <- c("HH","individual_ID")
      } else if ("individual_ID" %in% join_keys) {
        keys <- "individual_ID"
      } else {
        stop("`Covariate_DataFrame` must contain at least `individual_ID` (and ideally `HH`, and/or `test_date`).")
      }

      # Avoid overwriting base df columns
      keep_cols <- setdiff(names(cv), keys)
      if (!is.null(covariate_cols)) keep_cols <- intersect(keep_cols, covariate_cols)
      if (length(keep_cols) > 0) {
        collide <- intersect(keep_cols, names(df))
        if (length(collide)) {
          ren <- setNames(paste0(collide, "_cv"), collide)
          names(cv)[match(collide, names(cv))] <- ren
          keep_cols <- setdiff(names(cv), keys)  # refresh
        }
        df <- merge(df, cv[, c(keys, keep_cols), drop = FALSE], by = keys, all.x = TRUE, sort = FALSE)
      }
    }

    ##########################################################################
    #                                                                        #
    #                  Decide which covariates to summarize                  #
    #                                                                        #
    ##########################################################################
    # Anything beyond the 6 core columns counts as a covariate
    extra_cols <- setdiff(names(df), core_cols)
    if (!is.null(covariate_cols)) extra_cols <- intersect(extra_cols, covariate_cols)

    sanitize <- function(x) {
      x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
      gsub("_+", "_", gsub("^_|_$", "", x))
    }

    mode_val <- function(v) {
      v <- v[!is.na(v)]
      if (!length(v)) return(NA)
      tt <- sort(table(v), decreasing = TRUE)
      names(tt)[1]
    }

    ##########################################################################
    #                                                                        #
    #                     Split by person and summarize                      #
    #                                                                        #
    ##########################################################################
    split_df <- split(df, list(df$HH, df$individual_ID), drop = TRUE)

    res_list <- lapply(split_df, function(ind_df) {
      status <- ind_df$infection_status
      dates  <- ind_df$test_date
      role   <- unique(ind_df$role)
      comm_risk <- mean(ind_df$community_risk, na.rm = TRUE)

      r   <- rle(status)
      idx_end   <- cumsum(r$lengths)
      idx_start <- idx_end - r$lengths + 1

      starts <- dates[idx_start[r$values == 1]]
      ends   <- dates[idx_end[r$values == 1]]
      n_inf  <- sum(r$values == 1)
      dur    <- ends - starts + 1

      first_pos <- if (length(starts)) min(starts) else Inf
      has_pre_neg <- ind_df$infection_status == 0 & ind_df$test_date < first_pos
      last_negative <- if (any(has_pre_neg)) max(ind_df$test_date[has_pre_neg]) else NA_integer_

      infectious_day <- if (length(starts) > 0) paste(starts[1]:ends[1], collapse = ",") else NA_character_

      out <- list(
        HH = ind_df$HH[1],
        individual_ID = ind_df$individual_ID[1],
        n.true.infection = n_inf,
        n.detected.infection = n_inf,
        infection.detected.start = if (length(starts) > 0) starts[1] else NA_integer_,
        infection.detected.end   = if (length(ends) > 0)   ends[1]   else NA_integer_,
        infection.true.duration  = if (length(dur) > 0)    dur[1]    else NA_integer_,
        last_negative = last_negative,
        infection.infectious.day = infectious_day,
        community.risk = comm_risk,
        role = role
      )

      # Covariate summaries, if there are any
      if (length(extra_cols) > 0) {
        for (cn in extra_cols) {
          v <- ind_df[[cn]]
          sn <- sanitize(cn)

          if (all(is.na(v))) {
            out[[paste0(sn, "_timevarying")]] <- NA
            out[[paste0(sn, "_mode")]]  <- NA
            out[[paste0(sn, "_first")]] <- NA
            out[[paste0(sn, "_mean")]]  <- NA
            out[[paste0(sn, "_last")]]  <- NA
            next
          }

          is_num  <- is.numeric(v)
          is_logi <- is.logical(v)
          is_cat  <- is.factor(v) || is.character(v) || is_logi

          is_bin_num <- FALSE
          if (is_num) {
            u <- unique(v[!is.na(v)])
            is_bin_num <- length(u) <= 2 && all(u %in% c(0,1))
          }

          if (is_cat || is_bin_num) {
            mv <- mode_val(if (is_logi) as.integer(v) else v)
            out[[paste0(sn, "_mode")]] <- if (is_logi || is_bin_num) as.numeric(mv) else as.character(mv)
          } else if (is_num) {
            nz <- which(!is.na(v))
            out[[paste0(sn, "_first")]] <- v[nz[1]]
            out[[paste0(sn, "_mean")]]  <- mean(v, na.rm = TRUE)
            out[[paste0(sn, "_last")]]  <- v[nz[length(nz)]]
          }

          out[[paste0(sn, "_timevarying")]] <- length(unique(v[!is.na(v)])) > 1
        }
      }

      as.data.frame(out, stringsAsFactors = FALSE)
    })

    out_df <- do.call(rbind, res_list)
    rownames(out_df) <- NULL
    return(out_df)

  }, error = function(e) {
    message("An error occurred during summarization: ", conditionMessage(e))
    return(NULL)
  })
}



#' Try episode-level summarization (safe helper)
#'
#' Runs \code{data_summarization()} only when core columns are present;
#' otherwise returns \code{NULL}.
#'
#' @param long_df Data frame expected to contain `HH`, `individual_ID`, `role`,
#'   `test_date`, `infection_status`, `community_risk`.
#' @param Covariate_DataFrame Optional covariate frame to merge beforehand.
#' @param covariate_cols Optional character vector restricting covariate summary.
#'
#' @return A data frame (as from \code{data_summarization()}) or \code{NULL}
#'   when requirements are not met.
#' @keywords internal
.try_episode_summary <- function(long_df, Covariate_DataFrame = NULL, covariate_cols = NULL) {
  core <- c("HH","individual_ID","role","test_date","infection_status","community_risk")
  if (!is.null(long_df) && is.data.frame(long_df) && all(core %in% names(long_df))) {
    return(
      data_summarization(
        df = long_df,
        Covariate_DataFrame = Covariate_DataFrame,
        covariate_cols = covariate_cols
      )
    )
  }
  NULL
}



#' Minimal legacy-like summaries from RSV/VL households (helper)
#'
#' Converts RSV/VL simulator households to (i) individual summaries and
#' (ii) a minimal person-day table with viral-load test values for plotting.
#'
#' @param households List of per-household data frames from the RSV/VL simulator;
#'   each should include `hh_id`, `role`, `infectious_start`, `infectious_end`,
#'   `vl_full_trajectory`, and attribute `test_days`.
#' @param seasonal_forcing_list Named list of role vectors (`adult`, `child`,
#'   `elderly`, `toddler`) passed to \code{households_to_long_tests()}.
#' @param start_date,end_date \code{Date} range used internally by legacy
#'   summarization/imputation.
#'
#' @return List with:
#'   \itemize{
#'     \item \code{summarized_data}: output of legacy summarization/imputation.
#'     \item \code{person_day}: data frame with `HH`, `individual_ID`, `day`,
#'       and `vl_test` (for SAR-by-VL plotting).
#'   }
#' @details Internally converts households to long test days, then reuses
#'   legacy summarization and imputation.
#' @keywords internal
.stan_make_summary_from_households <- function(households,
                                               seasonal_forcing_list,
                                               start_date, end_date) {
  # 1) Convert the simulator-style list to a long test-day table
  long_tests <- households_to_long_tests(households)

  # 2) Reuse your legacy summarization + imputation (gives HH/ID/roles/inf_* cols)
  dt <- summarize_individuals(raw_dt = list(long_tests),
                              study_start = start_date, study_end = end_date)
  dt <- infectious_time_imputation(dt = dt, study_start = start_date)

  # 3) Build a *minimal* person-day table for SAR-by-VL (just what the plot needs)
  #    Expected by .plot_sar_by_index_vl(): HH, individual_ID, day, and a VL column (vl_test)
  pd_min <- long_tests |>
    dplyr::transmute(HH, individual_ID, day = test_date, vl_test)

  list(summarized_data = dt, person_day = pd_min)
}



#' Flatten simulator households to a long test-day table
#'
#' Expands simulator output to per-person per-day rows suitable for legacy
#' processing and plotting.
#'
#' @param households List of per-household data frames (legacy simulator).
#'   Each should include `individual_ID`, `role`, `infection_status`,
#'   `infectious_start`, `infectious_end`. Optional: `infection_time`,
#'   `infection_resolved`, scalar `community_risk`, and attribute `test_days`.
#'
#' @return Data frame with columns:
#'   `HH`, `individual_ID`, `role`, `test_date`, `infection_status`
#'   and (if available) `community_risk`.
#'
#' @details If `test_days` is absent, the horizon is inferred from
#' `infectious_end`/`infection_resolved` (fallback to 1). Missing `infection_status`
#' is treated as 0. Household index in the list is used for `HH`.
households_to_long_tests <- function(households) {
  if (is.null(households) || length(households) == 0L) {
    return(data.frame(
      HH = integer(0), individual_ID = integer(0), role = character(0),
      test_date = integer(0), infection_status = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  out_list <- vector("list", length(households))

  for (h in seq_along(households)) {
    hh <- households[[h]]
    if (is.null(hh) || nrow(hh) == 0L) next

    # --- Minimal normalization / safety ---
    # Ensure required columns exist
    need_cols <- c("individual_ID", "role", "infection_status",
                   "infectious_start", "infectious_end")
    for (nm in need_cols) if (!nm %in% names(hh)) hh[[nm]] <- NA

    # Coerce infection_status to {0,1} with NAs -> 0 (prevents if(NA) later)
    hh$infection_status <- as.integer(ifelse(is.na(hh$infection_status), 0L, hh$infection_status))

    # Compute a safe max_days:
    # try (a) attr test_days, (b) infectious_end, (c) infection_resolved
    test_days_attr <- attr(hh, "test_days", exact = TRUE)
    candidates <- c(
      suppressWarnings(as.integer(test_days_attr)),
      suppressWarnings(as.integer(hh$infectious_end)),
      suppressWarnings(as.integer(hh$infection_resolved))
    )
    max_days <- suppressWarnings(max(candidates, na.rm = TRUE))
    if (!is.finite(max_days) || is.na(max_days) || max_days < 1L) max_days <- 1L
    max_days <- as.integer(max_days)

    # Optional per-row scalar community risk (carry through if present)
    has_comm_risk <- "community_risk" %in% names(hh)

    # Pre-allocate a list of rows
    rows <- vector("list", nrow(hh) * max_days)
    k <- 0L

    for (i in seq_len(nrow(hh))) {
      id_i   <- hh$individual_ID[i]
      role_i <- hh$role[i]

      # Safe bounds
      s_i <- suppressWarnings(as.integer(hh$infectious_start[i]))
      e_i <- suppressWarnings(as.integer(hh$infectious_end[i]))

      # Status flag (explicit TRUE/FALSE, never NA)
      status1 <- isTRUE(hh$infection_status[i] == 1L)

      # Optional scalar community_risk for this person
      comm_i <- if (has_comm_risk) hh$community_risk[i] else NA

      # Build day rows
      for (d in seq_len(max_days)) {
        # In-range infectious if both start/end present; otherwise FALSE
        in_window <- (!is.na(s_i) && !is.na(e_i) && d >= s_i && d <= e_i)

        infect_today <- as.integer(status1 && in_window)

        k <- k + 1L
        if (has_comm_risk) {
          rows[[k]] <- list(HH = as.integer(h),
                            individual_ID = id_i,
                            role = role_i,
                            test_date = as.integer(d),
                            infection_status = infect_today,
                            community_risk = comm_i)
        } else {
          rows[[k]] <- list(HH = as.integer(h),
                            individual_ID = id_i,
                            role = role_i,
                            test_date = as.integer(d),
                            infection_status = infect_today)
        }
      }
    }

    hh_long <- do.call(rbind.data.frame, rows)
    if (!has_comm_risk) {
      hh_long$community_risk <- NULL
    }
    out_list[[h]] <- hh_long
  }

  out <- do.call(rbind.data.frame, out_list)
  rownames(out) <- NULL
  out
}
