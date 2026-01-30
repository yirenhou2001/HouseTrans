#' Daily infections plot (RSV/VL simulator)
#'
#' Aggregates simulated household data to daily new infections by role.
#'
#' @param households List of data frames from the RSV/VL simulator; each must
#'   contain \code{role} and \code{infection_time}.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.rsv_plot_daily <- function(households) {
  df <- do.call(rbind, households)
  df %>%
    dplyr::mutate(if_infection = !is.na(infection_time)) %>%
    dplyr::group_by(role, day = infection_time) %>%
    dplyr::summarise(n_infections = sum(if_infection), .groups = "drop") %>%
    ggplot2::ggplot(ggplot2::aes(x = day, y = n_infections, fill = role)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::labs(x = "Day of infection", y = "New infections", fill = "Age group",
                  title = "Daily new infections by age group") +
    ggplot2::theme_classic(base_size = 14)
}


#' Weekly infections plot (RSV/VL simulator)
#'
#' Aggregates simulated household data to weekly new infections by role.
#'
#' @param households List of data frames from the RSV/VL simulator; each must
#'   contain \code{role} and \code{infection_time}.
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.rsv_plot_weekly <- function(households) {
  df <- do.call(rbind, households)
  df %>%
    dplyr::filter(!is.na(infection_time)) %>%
    dplyr::mutate(week = floor(infection_time/7)) %>%
    dplyr::count(role, week, name = "n_infections") %>%
    ggplot2::ggplot(ggplot2::aes(x = week, y = n_infections, fill = role)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::labs(x = "Week", y = "New infections", fill = "Age group",
                  title = "Weekly new infections by age group") +
    ggplot2::theme_classic(base_size = 14)
}


#' Timeline plot of household epidemics (RSV/VL simulator)
#'
#' Shows infection and detection events per person across selected households.
#'
#' @param households List of data frames from the RSV/VL simulator; each must
#'   contain \code{hh_id}, \code{person_id}, \code{role}, \code{infection_time},
#'   and \code{detection_time}.
#' @param max_hh Integer. Maximum number of households to display (default 15).
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.rsv_plot_timeline <- function(households, max_hh = 15) {
  df  <- do.call(rbind, households)
  sel <- unique(df$hh_id)[seq_len(min(length(unique(df$hh_id)), max_hh))]
  tl  <- df %>%
    dplyr::filter(hh_id %in% sel) %>%
    tidyr::pivot_longer(c(infection_time, detection_time),
                        names_to = "event_type", values_to = "day") %>%
    dplyr::filter(!is.na(day)) %>%
    dplyr::mutate(event_type = factor(event_type, levels = c("detection_time","infection_time")))

  ggplot2::ggplot(tl, ggplot2::aes(x = day, y = factor(person_id), color = event_type, shape = role)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~ hh_id, scales = "free_y") +
    ggplot2::scale_color_manual(
      name   = "Event",
      values = c(detection_time = "#1f77b4",  # blue
                 infection_time  = "#d62728") # red
    ) +
    ggplot2::labs(x = "Day", y = "Person ID", title = "Household Epidemic Timelines") +
    ggplot2::theme_bw(base_size = 12)
}


#' SAR by index viral load (RSV/VL simulator)
#'
#' Computes household secondary attack rate (SAR) and displays SAR by the
#' index case's viral load category on detection day.
#'
#' @param households List of data frames from the RSV/VL simulator; each must
#'   include \code{hh_id}, \code{infection_time}, \code{detection_time},
#'   and \code{viral_loads_test_days} (per-person list column).
#'
#' @return A \code{ggplot} object.
#' @keywords internal
.rsv_plot_sar_by_index_vl <- function(households) {
  df <- do.call(rbind, households)
  hh_sum <- df %>%
    dplyr::group_by(hh_id) %>%
    dplyr::summarize(n_total = dplyr::n(),
                     n_infected = sum(!is.na(infection_time)),
                     sar = ifelse(n_total > 1, (n_infected - 1)/(n_total - 1), NA_real_),
                     .groups = "drop") %>%
    dplyr::mutate(sar = ifelse(sar < 0, NA, sar))

  get_index_vl <- function(hh) {
    if (all(is.na(hh$infection_time))) return(NA_real_)
    tmin <- min(hh$infection_time, na.rm = TRUE)
    idx  <- which(hh$infection_time == tmin)[1]
    vl   <- hh$viral_loads_test_days[[idx]]
    tdet <- suppressWarnings(min(hh$detection_time, na.rm = TRUE))
    if (!is.finite(tdet) || length(vl) == 0) return(NA_real_)
    as.numeric(vl[as.character(tdet)])
  }

  idx_vl <- tibble::tibble(
    hh_id    = vapply(households, function(hh) unique(hh$hh_id)[1], character(1)),
    index_vl = vapply(households, get_index_vl, numeric(1))
  )

  hh_sum %>%
    dplyr::left_join(idx_vl, by = "hh_id") %>%
    dplyr::mutate(vl_category = cut(index_vl, breaks = c(-Inf,2,4,6,8,Inf),
                                    labels = c("0-2","2-4","4-6","6-8",">8"))) %>%
    dplyr::filter(!is.na(vl_category)) %>%
    ggplot2::ggplot(ggplot2::aes(x = vl_category, y = sar)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15), alpha = 0.5) +
    ggplot2::labs(x = "Index VL (log10 copies/mL, bins)", y = "Household SAR",
                  title = "Household SAR by index-case viral load category") +
    ggplot2::theme_classic(base_size = 14)
}


#' Make transmission plots (batch)
#'
#' Builds selected figures from RSV/VL simulation outputs.
#'
#' @param results Pipeline result list.
#' @param which Character vector of plot names among \code{"daily"},
#'   \code{"weekly"}, \code{"sar"}, \code{"timeline"}, or \code{"all"}.
#' @param print Logical; print each plot when \code{TRUE} (default).
#' @param index_vl_column Optional character; viral-load column (unused, kept for API compatibility).
#'
#' @return Named list of \code{ggplot} objects for the requested plots
#'   (a subset of \code{daily}, \code{weekly}, \code{timeline}, \code{sar}).
#' @keywords internal
.make_transmission_plots <- function(results,
                                     which = c("daily","weekly","sar","timeline"),
                                     print = TRUE,
                                     index_vl_column = NULL) {
  if (identical(which, "all")) which <- c("daily","weekly","sar","timeline")
  out <- list()

  # Check for RSV/VL households
  households <- NULL

  if (!is.null(results$raw_simulation) &&
      is.list(results$raw_simulation) &&
      length(results$raw_simulation) &&
      is.data.frame(results$raw_simulation[[1]]) &&
      all(c("hh_id","person_id","role","infection_time","detection_time")
          %in% names(results$raw_simulation[[1]]))) {
    households <- results$raw_simulation
  } else if (!is.null(results$households) &&
             is.list(results$households) &&
             length(results$households) &&
             is.data.frame(results$households[[1]])) {
    households <- results$households
  }

  if (!is.null(households)) {
    if ("daily"    %in% which) out$daily    <- .rsv_plot_daily(households)
    if ("weekly"   %in% which) out$weekly   <- .rsv_plot_weekly(households)
    if ("timeline" %in% which) out$timeline <- .rsv_plot_timeline(households)
    if ("sar"      %in% which) out$sar      <- .rsv_plot_sar_by_index_vl(households)
  }

  if (print && length(out)) for (nm in names(out)) print(out[[nm]])
  out
}
