#' @title Plotting Functions
#' @description Visualization functions for model results and transmission chains
#' @name plotting_functions
#' @keywords internal
NULL


#' Plot Posterior Distributions (Phi and Kappa)
#'
#' Creates violin plots of posterior distributions for susceptibility and
#' infectivity parameters by role.
#'
#' @param fit A \code{stanfit} object.
#'
#' @return A combined ggplot object.
#' @keywords internal
plot_posterior_distributions <- function(fit) {
  post <- rstan::extract(fit)

  process_param <- function(matrix_data, param_name, labels) {
    df <- as.data.frame(matrix_data)
    colnames(df) <- labels
    df %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = "group", values_to = "value") %>%
      dplyr::mutate(
        log_value = log10(value),
        group = factor(group, levels = labels),
        type = param_name
      )
  }

  role_labels <- c("Adult", "Infant", "Toddler", "Elderly")
  df_phi <- process_param(post$phi_by_role, "Susceptibility", paste0("phi_", role_labels))
  df_kappa <- process_param(post$kappa_by_role, "Infectivity", paste0("kappa_", role_labels))

  role_colors <- c(
    "phi_Adult" = "#00A1D5FF", "phi_Infant" = "#79AF97FF",
    "phi_Toddler" = "#DF8F44FF", "phi_Elderly" = "#B24745FF",
    "kappa_Adult" = "#00A1D5FF", "kappa_Infant" = "#79AF97FF",
    "kappa_Toddler" = "#DF8F44FF", "kappa_Elderly" = "#B24745FF"
  )

  create_violin <- function(data, title_y) {
    medians <- data %>%
      dplyr::group_by(group) %>%
      dplyr::summarize(med = stats::median(log_value), .groups = "drop")

    ggplot2::ggplot(data, ggplot2::aes(x = group, y = log_value, fill = group)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
      ggplot2::geom_violin(trim = FALSE, alpha = 0.7) +
      ggplot2::geom_boxplot(width = 0.1, outliers = FALSE, alpha = 0.7, fill = "white") +
      ggplot2::geom_text(data = medians, ggplot2::aes(x = group, y = med, label = round(med, 2)),
                         position = ggplot2::position_nudge(x = 0.2), hjust = 0, size = 4) +
      ggplot2::scale_fill_manual(values = role_colors) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::labs(x = "", y = title_y) +
      ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank()) +
      ggplot2::scale_x_discrete(labels = role_labels)
  }

  p1 <- create_violin(df_phi, "log10 Susceptibility")
  p2 <- create_violin(df_kappa, "log10 Infectivity")

  # Combine plots
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    ggpubr::ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  } else {
    p1
  }
}


#' Plot Covariate Coefficients (Forest Plot)
#'
#' Creates a forest plot showing the effects of covariates on susceptibility
#' and infectivity.
#'
#' @param fit A \code{stanfit} object.
#' @param stan_data The list data passed to Stan.
#'
#' @return A ggplot object.
#' @keywords internal
plot_covariate_effects <- function(fit, stan_data) {
  post <- rstan::extract(fit)

  process_beta <- function(beta_matrix, col_names, type_label) {
    if (is.null(beta_matrix)) return(NULL)
    df <- as.data.frame(beta_matrix)
    if (ncol(df) == 1 && length(col_names) == 1) colnames(df) <- col_names

    df %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = "covariate", values_to = "val") %>%
      dplyr::mutate(type = type_label)
  }

  df_susc <- NULL
  if (stan_data$K_susc > 0) {
    names_susc <- paste0("Susc_", 1:stan_data$K_susc)
    df_susc <- process_beta(post$beta_susc, names_susc, "Susceptibility")
  }

  df_inf <- NULL
  if (stan_data$K_inf > 0) {
    names_inf <- paste0("Inf_", 1:stan_data$K_inf)
    df_inf <- process_beta(post$beta_inf, names_inf, "Infectivity")
  }

  df_all <- rbind(df_susc, df_inf)

  if (is.null(df_all)) {
    message("No covariates found in model.")
    return(NULL)
  }

  ggplot2::ggplot(df_all, ggplot2::aes(x = val, y = covariate, color = type)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::stat_summary(
      fun.data = function(x) {
        data.frame(
          y = stats::median(x),
          ymin = stats::quantile(x, 0.025),
          ymax = stats::quantile(x, 0.975)
        )
      },
      geom = "pointrange", size = 1
    ) +
    ggplot2::labs(
      title = "Covariate Effects (Log Scale)",
      subtitle = "Negative = Protective",
      x = "Log Coefficient (Beta)", y = ""
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::scale_color_brewer(palette = "Set1")
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
#' @keywords internal
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


#' Plot Household Timeline with Transmission Chains
#'
#' Creates a timeline visualization of infections within a household,
#' showing transmission links with probabilities.
#'
#' @param trans_df Data frame from \code{reconstruct_transmission_chains}.
#' @param stan_data The list data passed to Stan.
#' @param target_hh_id Integer; household ID to plot.
#' @param start_date_str Character; start date string.
#' @param prob_cutoff Numeric; minimum probability for displaying links.
#' @param plot_width,plot_height Numeric; plot dimensions.
#'
#' @return A ggplot object.
#' @keywords internal
plot_household_timeline <- function(trans_df, stan_data, target_hh_id,
                                    start_date_str = "2024-10-08",
                                    prob_cutoff = 0,
                                    plot_width = 11,
                                    plot_height = 7) {

  target_hh_id <- as.integer(target_hh_id)
  start_date <- as.Date(start_date_str)

  # --- 1. PARAMETERS ---
  plot_aspect_ratio <- plot_width / plot_height
  curvature_val <- 0.3
  label_sag_factor <- 0.5

  # --- 2. DATA PREP ---
  hh_n_indices <- stan_data$hh_members[target_hh_id, ]
  hh_n_indices <- hh_n_indices[hh_n_indices != 0]

  if (length(hh_n_indices) == 0) return(.empty_plot("No members in household"))

  raw_p_ids <- if (!is.null(stan_data$p_id)) stan_data$p_id[hh_n_indices] else hh_n_indices

  map_n_to_p <- data.frame(
    n_id = hh_n_indices,
    p_id = raw_p_ids,
    role_id = stan_data$role_id[hh_n_indices],
    stringsAsFactors = FALSE
  )

  role_map_labels <- c("Adult", "Infant", "Toddler", "Elderly")

  roster_df <- map_n_to_p %>%
    dplyr::select(p_id, role_id) %>%
    dplyr::distinct() %>%
    dplyr::mutate(role_name = role_map_labels[role_id]) %>%
    dplyr::group_by(role_name) %>%
    dplyr::mutate(
      count_in_role = dplyr::n(),
      rank_in_role = dplyr::row_number(),
      label = ifelse(count_in_role > 1, paste(role_name, rank_in_role), role_name)
    ) %>%
    dplyr::ungroup()

  events_list <- list()
  for (i in 1:nrow(map_n_to_p)) {
    curr_n <- map_n_to_p$n_id[i]
    curr_p <- map_n_to_p$p_id[i]
    status_vec <- stan_data$I[curr_n, ]
    is_inf <- which(status_vec == 1)
    if (length(is_inf) > 0) {
      gaps <- c(0, diff(is_inf))
      starts <- is_inf[gaps != 1]
      if (length(starts) == 0) starts <- is_inf[1]
      for (k in seq_along(starts)) {
        events_list[[length(events_list) + 1]] <- data.frame(
          n_id = curr_n, p_id = curr_p, day_idx = starts[k], stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(events_list) == 0) return(.empty_plot("No infections in household"))
  events_df <- do.call(rbind, events_list)

  events_df <- events_df %>%
    dplyr::arrange(p_id, day_idx) %>%
    dplyr::group_by(p_id) %>%
    dplyr::mutate(type = ifelse(dplyr::row_number() == 1, "Primary", "Reinfection")) %>%
    dplyr::ungroup()

  layout_order <- events_df %>%
    dplyr::filter(type == "Primary") %>%
    dplyr::select(p_id, first_inf = day_idx) %>%
    dplyr::distinct() %>%
    dplyr::arrange(first_inf, p_id) %>%
    dplyr::mutate(y_rank = dplyr::row_number())

  max_rank <- max(layout_order$y_rank)
  layout_order$y_plot <- max_rank - layout_order$y_rank + 1

  plot_data <- events_df %>%
    dplyr::left_join(layout_order %>% dplyr::select(p_id, y_plot), by = "p_id") %>%
    dplyr::left_join(roster_df %>% dplyr::select(p_id, label, role_name), by = "p_id") %>%
    dplyr::mutate(date = start_date + day_idx) %>%
    dplyr::filter(!is.na(y_plot))

  y_labels <- plot_data %>%
    dplyr::select(y_plot, label) %>%
    dplyr::distinct() %>%
    dplyr::arrange(y_plot)

  # --- Filter Links ---
  hh_links <- trans_df %>%
    dplyr::filter(as.integer(hh_id) == target_hh_id, prob >= prob_cutoff)

  # Coordinate joins
  target_coords <- plot_data %>%
    dplyr::select(n_id, target_x = date, target_y = y_plot)
  edges_prep <- hh_links %>%
    dplyr::mutate(target = as.integer(target)) %>%
    dplyr::inner_join(target_coords, by = c("target" = "n_id"))

  source_coords <- plot_data %>%
    dplyr::select(n_id, source_x = date, source_y = y_plot, source_role = role_name)

  edges_hh <- edges_prep %>%
    dplyr::filter(source != "Community") %>%
    dplyr::mutate(source_int = as.integer(source)) %>%
    dplyr::inner_join(source_coords, by = c("source_int" = "n_id")) %>%
    dplyr::select(-source_int)

  comm_x <- min(plot_data$date) - 5
  comm_y <- max(plot_data$y_plot) + 0.5

  edges_comm <- edges_prep %>%
    dplyr::filter(source == "Community") %>%
    dplyr::mutate(source_x = comm_x, source_y = comm_y, source_role = "Community")

  all_edges <- dplyr::bind_rows(edges_hh, edges_comm)

  # --- 3. ROBUST GEOMETRY CALCULATION ---
  min_date_val <- min(c(plot_data$date, comm_x), na.rm = TRUE)
  max_date_val <- max(c(plot_data$date, comm_x), na.rm = TRUE)
  range_x <- as.numeric(max_date_val - min_date_val)
  if (range_x == 0) range_x <- 1

  min_y_val <- 0.5
  max_y_val <- max_rank + 1
  range_y <- max_y_val - min_y_val

  if (nrow(all_edges) > 0) {
    all_edges <- all_edges %>%
      dplyr::mutate(
        x1 = as.numeric(source_x), y1 = source_y,
        x2 = as.numeric(target_x), y2 = target_y,

        ny1 = (y1 - min_y_val) / range_y,
        ny2 = (y2 - min_y_val) / range_y,

        nx1 = ((x1 - as.numeric(min_date_val)) / range_x) * plot_aspect_ratio,
        nx2 = ((x2 - as.numeric(min_date_val)) / range_x) * plot_aspect_ratio,

        ndx = nx2 - nx1,
        ndy = ny2 - ny1,
        ndist = sqrt(ndx^2 + ndy^2),

        nmx = (nx1 + nx2) / 2,
        nmy = (ny1 + ny2) / 2,

        noffset = label_sag_factor * curvature_val * ndist,
        nshift_x = (ndy / ndist) * noffset,
        nshift_y = (-ndx / ndist) * noffset,
        nfinal_x = nmx + nshift_x,
        nfinal_y = nmy + nshift_y,

        lbl_x_num = (nfinal_x / plot_aspect_ratio) * range_x + as.numeric(min_date_val),
        lbl_y     = (nfinal_y * range_y) + min_y_val,
        lbl_x = as.Date(lbl_x_num, origin = "1970-01-01")
      )

    min_lbl_y <- min(all_edges$lbl_y, na.rm = TRUE)
    lower_limit <- min(0.5, min_lbl_y - 0.5)

    all_edges <- all_edges %>%
      dplyr::mutate(alpha_val = dplyr::case_when(
        prob <= 0.5 ~ 0.5,
        TRUE        ~ 1.0
      ))
  } else {
    lower_limit <- 0.5
  }

  role_colors <- c(
    "Adult"   = "#00A1D5FF",
    "Infant"  = "#79AF97FF",
    "Toddler" = "#DF8F44FF",
    "Elderly" = "#B24745FF",
    "Community" = "black"
  )

  # --- 4. PLOTTING ---
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = comm_x, y = comm_y + 0.25, label = "Community Risk",
                      fontface = "bold", color = "black", size = 8, hjust = 0) +
    ggplot2::annotate("point", x = comm_x, y = comm_y, shape = 18, size = 10, color = "black", alpha = 0.7)

  if (nrow(all_edges) > 0) {
    p <- p +
      ggplot2::geom_curve(
        data = all_edges,
        ggplot2::aes(x = source_x, y = source_y, xend = target_x, yend = target_y,
                     color = source_role, alpha = alpha_val),
        curvature = curvature_val,
        linetype = "dashed",
        linewidth = 0.8,
        arrow = ggplot2::arrow(length = ggplot2::unit(0.5, "cm"), type = "closed")
      ) +
      ggplot2::geom_label(
        data = all_edges,
        ggplot2::aes(x = lbl_x, y = lbl_y,
                     label = scales::percent(prob, accuracy = 1),
                     alpha = alpha_val),
        color = NA, fill = "white", size = 7.5,
        label.padding = ggplot2::unit(0.1, "lines")
      ) +
      ggplot2::geom_text(
        data = all_edges,
        ggplot2::aes(x = lbl_x, y = lbl_y,
                     label = scales::percent(prob, accuracy = 1),
                     color = source_role, alpha = alpha_val),
        size = 7.5
      )
  }

  p <- p +
    ggplot2::geom_point(
      data = plot_data %>% dplyr::filter(type == "Primary"),
      ggplot2::aes(x = date, y = y_plot, color = role_name),
      size = 10, alpha = 0.5
    ) +
    ggplot2::geom_point(
      data = plot_data %>% dplyr::filter(type == "Reinfection"),
      ggplot2::aes(x = date, y = y_plot, color = role_name),
      size = 10, alpha = 0.5
    ) +
    ggplot2::scale_x_date(date_labels = "%b %d", breaks = scales::breaks_pretty(n = 8)) +
    ggplot2::scale_y_continuous(breaks = y_labels$y_plot, labels = y_labels$label,
                                limits = c(lower_limit, max_rank + 1)) +
    ggplot2::scale_color_manual(values = role_colors) +
    ggplot2::scale_alpha_identity() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = paste("Household", target_hh_id), x = "", y = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      axis.text.y = ggplot2::element_text(face = "bold", size = 23, color = "black"),
      axis.text.x = ggplot2::element_text(face = "bold", size = 23, color = "black"),
      legend.position = "none"
    )

  return(p)
}


#' Create Empty Plot with Error Message
#'
#' @param message Error message to display.
#' @return A ggplot object.
#' @keywords internal
#' @noRd
.empty_plot <- function(message = "Plot not available") {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = message,
                      size = 5, color = "gray50") +
    ggplot2::theme_void() +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1)
}
