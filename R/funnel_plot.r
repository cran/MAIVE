#' Funnel plot (base graphics)
#'
#' Components provided:
#' - `get_funnel_plot()` as the public high-level entrypoint
#' - internal helpers to validate and map `(dat, result)` into plot inputs
#' - an internal base-graphics engine that draws onto the current device
#'
#' Device management (PNG/SVG/PDF), resolution, and any base64 encoding are the caller's
#' responsibility.
#'
#' @name maive_funnel_plot_helpers
#' @keywords internal
NULL

#' Minimum point size scale for WAIVE-adjusted points
#'
#' Used when downweighting adjusted funnel points (scaled between this value and 1.0).
#' @keywords internal
#' @noRd
MAIVE_ADJUSTED_POINT_MIN_SCALE <- 0.1

#' Default visual options for the funnel plot
#'
#' @return A list of plotting options.
#' @keywords internal
#' @noRd
maive_funnel_plot_opts <- function() {
  list(
    # Plot options
    xlab = "Effect Size",
    ylab = "Standard Error",
    yaxis = "sei", # seinv for precision
    digits = 3L,
    col = "black",
    # Legend options
    effect_shades = c("white", "black"),
    effect_pch = c(21, 19), # 21 is hollow circle, 19 is filled circle
    maive_shades = c("black", "black"),
    maive_pch = c(NA, NA), # NA for lines, will be overridden in legend
    ci_levels = c(90, 95, 99),
    ci_shades = c("white", "gray55", "gray75"),
    legend_texts = c(
      "Base effect",
      "Adjusted SE",
      "MAIVE fit",
      "95% CI bounds"
    ),
    pt_cex = 2,
    legend_inset = 0.01,
    text_color = "black",
    legend_bg = "white",
    legend_position = "bottomright",
    legend_bty = "o"
  )
}

#' Compute x-axis padding for funnel plot limits
#'
#' Dynamically sets padding so the right side is wider (to reduce legend overlap).
#' @param effect Numeric vector of effects.
#' @return List with `lower` and `upper` x-axis limits.
#' @keywords internal
#' @noRd
maive_funnel_padding <- function(effect) {
  # Set padding at either side - between 0 and 1, percentage of the effect range
  left_padding <- 0.05
  right_padding <- 0.25

  effect_all <- effect
  min_effect <- suppressWarnings(min(effect_all, na.rm = TRUE))
  max_effect <- suppressWarnings(max(effect_all, na.rm = TRUE))

  if (!is.finite(min_effect) || !is.finite(max_effect)) {
    return(list(lower = -1, upper = 1))
  }

  effect_range <- max_effect - min_effect
  if (!is.finite(effect_range) || effect_range <= 0) {
    effect_range <- max(abs(effect_all), na.rm = TRUE)
    if (!is.finite(effect_range) || effect_range == 0) {
      effect_range <- 1
    }
  }

  xlim_pad_left <- left_padding * effect_range
  xlim_pad_right <- right_padding * effect_range
  list(
    lower = min_effect - xlim_pad_left,
    upper = max_effect + xlim_pad_right
  )
}

#' Format axis labels with adaptive decimals
#'
#' @param ticks Numeric vector of tick positions.
#' @param digits Maximum decimals.
#' @return Character vector of tick labels.
#' @keywords internal
#' @noRd
maive_funnel_format_axis_labels <- function(ticks, digits = 3) {
  finite_ticks <- ticks[is.finite(ticks)]
  if (length(finite_ticks) == 0) {
    return(character(0))
  }

  if (all(abs(finite_ticks - round(finite_ticks)) < .Machine$double.eps^0.5)) {
    return(sprintf("%d", as.integer(round(ticks))))
  }

  unique_ticks <- sort(unique(finite_ticks))
  decimals <- digits
  if (length(unique_ticks) > 1) {
    diffs <- diff(unique_ticks)
    diffs <- diffs[diffs > 0]
    if (length(diffs) > 0) {
      min_diff <- min(diffs)
      if (is.finite(min_diff) && min_diff > 0) {
        decimals <- max(0, min(digits, ceiling(-log10(min_diff))))
      }
    }
  }

  formatC(ticks, format = "f", digits = decimals)
}

#' First non-NULL value helper
#'
#' @param ... Values to scan.
#' @return First non-NULL value, or NULL.
#' @keywords internal
#' @noRd
maive_funnel_first_non_null <- function(...) {
  vals <- list(...)
  for (val in vals) {
    if (!is.null(val)) {
      return(val)
    }
  }
  NULL
}

#' Normalize slope metadata structure for plotting
#'
#' Accepts either a structured slope list or legacy pieces (`slope_coef`, `is_quadratic_fit`).
#' @param primary Structured slope metadata list.
#' @param fallback_coef Legacy slope coefficient.
#' @param fallback_summary Legacy slope summary (logical or list).
#' @return Normalized slope metadata list.
#' @keywords internal
#' @noRd
maive_funnel_normalize_slope <- function(primary, fallback_coef, fallback_summary) {
  metadata <- list(
    type = NULL,
    quadratic = FALSE,
    coefficient = fallback_coef,
    detail = NULL
  )

  if (is.list(primary)) {
    metadata$type <- maive_funnel_first_non_null(primary$type, primary$slope_type)
    if (!is.null(primary$quadratic)) {
      metadata$quadratic <- isTRUE(primary$quadratic)
    }
    if (!is.null(primary$coefficient)) {
      metadata$coefficient <- primary$coefficient
    }
    metadata$detail <- maive_funnel_first_non_null(primary$detail, primary$slope_detail)
  } else if (is.logical(primary) && length(primary) == 1) {
    metadata$quadratic <- isTRUE(primary)
  }

  if (is.list(fallback_summary)) {
    if (!is.null(fallback_summary$quadratic)) {
      metadata$quadratic <- isTRUE(fallback_summary$quadratic)
    }
    if (is.null(metadata$type)) {
      metadata$type <- maive_funnel_first_non_null(fallback_summary$slope_type, fallback_summary$type)
    }
    if (is.null(metadata$detail)) {
      metadata$detail <- maive_funnel_first_non_null(fallback_summary$slope_detail, fallback_summary$detail)
    }
    if (is.null(metadata$coefficient) && !is.null(fallback_summary$coefficient)) {
      metadata$coefficient <- fallback_summary$coefficient
    }
  } else if (is.logical(fallback_summary) && length(fallback_summary) == 1) {
    metadata$quadratic <- isTRUE(fallback_summary)
  }

  if (is.null(metadata$type)) {
    metadata$type <- if (isTRUE(metadata$quadratic)) "quadratic" else "linear"
  }

  metadata$quadratic <- isTRUE(metadata$quadratic)
  metadata
}

# ---- High-level input prep helpers (internal) --------------------------------

#' Validate and extract core effect/SE vectors from `dat`
#'
#' @param dat Data frame containing numeric columns `bs` and `sebs`.
#' @return List with `effect` and `se` vectors.
#' @keywords internal
#' @noRd
maive_funnel_validate_dat <- function(dat) {
  if (is.null(dat)) {
    stop("dat must be provided.")
  }
  dat <- as.data.frame(dat)
  if (!all(c("bs", "sebs") %in% names(dat))) {
    stop("dat must contain numeric columns named 'bs' and 'sebs'.")
  }

  effect <- dat[["bs"]]
  se <- dat[["sebs"]]

  if (!is.numeric(effect) || !is.numeric(se)) {
    stop("dat$bs and dat$sebs must be numeric.")
  }
  if (length(effect) != length(se)) {
    stop("dat$bs and dat$sebs must have the same length.")
  }
  if (length(effect) == 0L) {
    stop("dat must contain at least one observation.")
  }
  if (any(!is.finite(effect)) || any(!is.finite(se))) {
    stop("dat$bs and dat$sebs must be finite.")
  }
  if (any(se <= 0)) {
    stop("dat$sebs must be strictly positive.")
  }

  list(effect = effect, se = se)
}

#' Validate MAIVE result object contains required fields
#'
#' @param result List returned by `maive()` or `waive()`.
#' @return Invisibly TRUE.
#' @keywords internal
#' @noRd
maive_funnel_validate_result <- function(result) {
  if (is.null(result)) {
    stop("result must be provided.")
  }
  if (!is.list(result)) {
    stop("result must be a list as returned by maive()/waive().")
  }
  if (is.null(result[["beta"]]) || is.null(result[["SE"]])) {
    stop("result must contain 'beta' and 'SE' fields.")
  }
  invisible(TRUE)
}

#' Normalize model type label
#'
#' @param model_type Character label.
#' @return Upper-cased model label (defaults to "MAIVE").
#' @keywords internal
#' @noRd
maive_funnel_normalize_model_type <- function(model_type) {
  if (is.null(model_type) || !nzchar(model_type)) {
    return("MAIVE")
  }
  toupper(as.character(model_type))
}

#' Normalize instrumentation flag
#'
#' @param instrument Optional indicator (0/1).
#' @param result MAIVE result (used for inference when `instrument` is NULL).
#' @return Integer 0 or 1.
#' @keywords internal
#' @noRd
maive_funnel_normalize_instrument <- function(instrument, result) {
  if (!is.null(instrument)) {
    inst <- suppressWarnings(as.integer(instrument))
    if (!is.na(inst) && inst %in% c(0L, 1L)) {
      return(inst)
    }
  }

  # Infer from the first-stage diagnostic when available:
  # - maive() sets F-test to "NA" (string) when instrument is disabled.
  ftest <- result[["F-test"]]
  if (is.null(ftest)) {
    return(0L)
  }
  if (is.character(ftest) && length(ftest) == 1L && identical(ftest, "NA")) {
    return(0L)
  }
  if (is.numeric(ftest) && (length(ftest) != 1L || is.na(ftest) || !is.finite(ftest))) {
    return(0L)
  }

  1L
}

#' Extract adjusted standard errors for plotting
#'
#' @param result MAIVE result.
#' @param instrument Integer 0/1.
#' @param n_points Expected vector length.
#' @return Numeric vector or NULL.
#' @keywords internal
#' @noRd
maive_funnel_extract_adjusted_se <- function(result, instrument, n_points) {
  if (instrument == 0L) {
    return(NULL)
  }
  se_adjusted <- result[["SE_instrumented"]]
  if (is.null(se_adjusted)) {
    return(NULL)
  }
  if (!is.numeric(se_adjusted) || length(se_adjusted) != n_points) {
    return(NULL)
  }
  if (!any(is.finite(se_adjusted))) {
    return(NULL)
  }
  se_adjusted
}

#' Extract WAIVE point weights for adjusted-point sizing
#'
#' @param result MAIVE/WAIVE result.
#' @param model_type Normalized model type.
#' @param instrument Integer 0/1.
#' @param n_points Expected vector length.
#' @return Numeric vector or NULL.
#' @keywords internal
#' @noRd
maive_funnel_extract_adjusted_weights <- function(result, model_type, instrument, n_points) {
  if (instrument == 0L) {
    return(NULL)
  }
  if (!identical(model_type, "WAIVE")) {
    return(NULL)
  }
  w <- result[["weights"]]
  if (is.null(w) || !is.numeric(w) || length(w) != n_points) {
    return(NULL)
  }
  if (!any(is.finite(w)) || !any(w > 0, na.rm = TRUE)) {
    return(NULL)
  }
  w
}

#' Extract slope metadata from MAIVE result
#'
#' @param result MAIVE result.
#' @return Slope metadata list compatible with the plot engine.
#' @keywords internal
#' @noRd
maive_funnel_extract_slope <- function(result) {
  slope_summary <- result[["is_quadratic_fit"]]
  list(
    type = if (is.list(slope_summary)) slope_summary[["slope_type"]] else NULL,
    quadratic = if (is.list(slope_summary)) slope_summary[["quadratic"]] else FALSE,
    coefficient = result[["slope_coef"]],
    detail = if (is.list(slope_summary)) slope_summary[["slope_detail"]] else NULL
  )
}

# Internal engine: draws the funnel plot on the current device.
#
# @keywords internal
maive_funnel_plot_engine <- function(
    effect,
    se,
    se_adjusted = NULL,
    adjusted_weights = NULL,
    intercept = NULL,
    intercept_se = NULL,
    slope_coef = NULL,
    is_quadratic_fit = NULL,
    instrument = 1,
    slope = NULL,
    model_type = "MAIVE") {
  funnel_opts <- maive_funnel_plot_opts()
  model_label <- toupper(model_type)

  slope_meta <- maive_funnel_normalize_slope(
    primary = slope,
    fallback_coef = slope_coef,
    fallback_summary = is_quadratic_fit
  )

  n_points <- length(effect)
  use_adjusted <- instrument != 0 &&
    !is.null(se_adjusted) &&
    length(se_adjusted) == n_points &&
    any(is.finite(se_adjusted))

  x_values <- effect
  y_values <- se

  plot_pch <- rep(funnel_opts$effect_pch[1], n_points)
  point_bg <- rep(funnel_opts$effect_shades[1], n_points)
  point_map <- seq_len(n_points)
  point_is_adjusted <- rep(FALSE, n_points)

  if (use_adjusted) {
    x_values <- c(x_values, effect)
    y_values <- c(y_values, se_adjusted)
    plot_pch <- c(plot_pch, rep(funnel_opts$effect_pch[2], n_points))
    point_bg <- c(point_bg, rep(funnel_opts$effect_shades[2], n_points))
    point_map <- c(point_map, seq_len(n_points))
    point_is_adjusted <- c(point_is_adjusted, rep(TRUE, n_points))
  }

  finite_points <- is.finite(x_values) & is.finite(y_values)
  x_values <- x_values[finite_points]
  y_values <- y_values[finite_points]
  plot_pch <- plot_pch[finite_points]
  point_bg <- point_bg[finite_points]
  point_map <- point_map[finite_points]
  point_is_adjusted <- point_is_adjusted[finite_points]

  base_point_cex <- max(0.1, funnel_opts$pt_cex * 0.6)
  point_cex <- rep(base_point_cex, length(x_values))

  adjusted_cex <- NULL
  if (
    use_adjusted &&
      !is.null(adjusted_weights) &&
      length(adjusted_weights) == n_points
  ) {
    weights_numeric <- suppressWarnings(as.numeric(adjusted_weights))
    weights_numeric[!is.finite(weights_numeric) | weights_numeric < 0] <- 0

    if (any(weights_numeric > 0)) {
      positive_max <- max(weights_numeric, na.rm = TRUE)
      if (is.finite(positive_max) && positive_max > 0) {
        adjusted_indices <- which(point_is_adjusted)
        if (length(adjusted_indices) > 0) {
          point_weights <- weights_numeric[point_map[adjusted_indices]]
          point_weights[!is.finite(point_weights) | point_weights <= 0] <- 0

          if (any(point_weights > 0)) {
            # Normalize weights to [0, 1] range
            normalized_weights <- point_weights / positive_max
            normalized_weights <- pmin(pmax(normalized_weights, 0), 1)

            # Scale point sizes from MAIVE_ADJUSTED_POINT_MIN_SCALE to 1.0 times the base size
            point_scale <- MAIVE_ADJUSTED_POINT_MIN_SCALE +
              (1 - MAIVE_ADJUSTED_POINT_MIN_SCALE) * normalized_weights

            adjusted_cex <- base_point_cex * point_scale
            point_cex[adjusted_indices] <- adjusted_cex
          }
        }
      }
    }
  }

  simple_mean <- mean(effect, na.rm = TRUE)

  padding <- maive_funnel_padding(effect)

  se_all <- y_values
  max_se <- max(se_all, na.rm = TRUE)
  if (!is.finite(max_se) || max_se <= 0) {
    max_se <- 1
  }

  se_pad <- max_se * 0.1
  ylim <- c(max_se + se_pad, 0)

  ci_levels <- funnel_opts$ci_levels
  ci_alpha <- 1 - (ci_levels / 100)
  valid_idx <- which(is.finite(ci_alpha) & ci_alpha > 0)
  ci_data <- NULL
  if (length(valid_idx) > 0) {
    ci_data <- data.frame(
      level = ci_levels[valid_idx],
      alpha = ci_alpha[valid_idx],
      z = stats::qnorm(1 - ci_alpha[valid_idx] / 2)
    )
  }
  max_z <- if (!is.null(ci_data) && nrow(ci_data) > 0) max(ci_data$z) else stats::qnorm(0.975)
  invisible(max_z)
  xlim <- c(padding$lower, padding$upper)
  if (!is.finite(xlim[1]) || !is.finite(xlim[2]) || xlim[1] == xlim[2]) {
    x_center <- ifelse(is.finite(simple_mean), simple_mean, 0)
    xlim <- c(x_center - 1, x_center + 1)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))

  graphics::plot(
    NA, NA,
    xlim = xlim,
    ylim = ylim,
    xlab = funnel_opts$xlab,
    ylab = funnel_opts$ylab,
    type = "n",
    axes = FALSE,
    xaxs = "i",
    yaxs = "i"
  )

  se_grid <- seq(0, max_se + se_pad, length.out = 400)

  shade_cols <- rep(funnel_opts$ci_shades, length.out = ifelse(is.null(ci_data), 0, nrow(ci_data)))
  contour_cols <- rep(c("gray90", "gray70", "gray50"), length.out = ifelse(is.null(ci_data), 0, nrow(ci_data)))
  outer_fill_col <- "gray90"
  outer_z <- NA_real_

  if (!is.null(ci_data) && nrow(ci_data) > 0) {
    names(shade_cols) <- ci_data$level
    names(contour_cols) <- ci_data$level

    outer_idx <- which.max(ci_data$z)
    outer_z <- ci_data$z[outer_idx]
    invisible(outer_z)

    p010_idx <- if (length(contour_cols) > 0) match(90, ci_data$level) else NA_integer_
    if (!is.na(p010_idx)) {
      outer_fill_col <- contour_cols[as.character(ci_data$level[p010_idx])]
    } else if (length(contour_cols) > 0) {
      outer_fill_col <- contour_cols[outer_idx]
    }

    graphics::rect(xlim[1], ylim[2], xlim[2], ylim[1], col = outer_fill_col, border = NA)

    draw_order <- order(ci_data$z, decreasing = TRUE)
    for (idx in draw_order) {
      level_value <- ci_data$level[idx]
      z_val <- ci_data$z[idx]
      left <- -z_val * se_grid
      right <- z_val * se_grid
      graphics::polygon(
        x = c(left, rev(right)),
        y = c(se_grid, rev(se_grid)),
        border = NA,
        col = shade_cols[as.character(level_value)]
      )
      graphics::lines(left, se_grid, col = contour_cols[as.character(level_value)], lwd = 1)
      graphics::lines(right, se_grid, col = contour_cols[as.character(level_value)], lwd = 1)
    }
  } else {
    graphics::rect(xlim[1], ylim[2], xlim[2], ylim[1], col = outer_fill_col, border = NA)
  }

  y_ticks <- pretty(c(0, max_se), n = 5)
  y_ticks <- unique(y_ticks[y_ticks >= 0 & y_ticks <= (max_se + se_pad)])
  if (length(y_ticks) > 0) {
    grid_col <- grDevices::adjustcolor("gray70", alpha.f = 0.4)
    graphics::abline(h = y_ticks, col = grid_col, lwd = 0.5)
  }

  draw_vertical_segment <- function(x_pos, lty, lwd, col) {
    if (!is.finite(x_pos)) {
      return()
    }
    graphics::segments(x_pos, ylim[1], x_pos, ylim[2], lty = lty, lwd = lwd, col = col)
  }

  if (!is.null(ci_data) && nrow(ci_data) > 0) {
    draw_vertical_segment(0, lty = 3, lwd = 1, col = "black")
  } else {
    graphics::abline(v = 0, lty = 3, col = "black")
  }

  draw_vertical_segment(simple_mean, lty = 4, lwd = 2, col = "black")

  se_curve_grid <- seq(0, max_se + se_pad, length.out = 200)

  if (!is.null(intercept) && !is.null(intercept_se) && !is.null(slope_meta$coefficient)) {
    slope_type <- tolower(as.character(maive_funnel_first_non_null(slope_meta$type, "linear")))

    draw_quadratic_or_linear <- function(power) {
      slope_value <- suppressWarnings(as.numeric(slope_meta$coefficient))
      if (length(slope_value) != 1 || !is.finite(slope_value)) {
        return()
      }
      se_term <- se_curve_grid^power
      x_pred <- intercept + slope_value * se_term
      graphics::lines(x_pred, se_curve_grid, lwd = 2, col = "black")

      vcov <- matrix(c(intercept_se^2, 0, 0, intercept_se^2), nrow = 2)
      se_fit <- sqrt(
        vcov[1, 1] +
          2 * se_term * vcov[1, 2] +
          se_term^2 * vcov[2, 2]
      )

      ci_lo <- x_pred - 1.96 * se_fit
      ci_hi <- x_pred + 1.96 * se_fit

      graphics::lines(ci_lo, se_curve_grid, lty = 2, col = "black")
      graphics::lines(ci_hi, se_curve_grid, lty = 2, col = "black")
    }

    draw_kinked <- function() {
      kink_effect <- NA_real_
      kink_location <- NA_real_

      if (is.list(slope_meta$coefficient)) {
        if (!is.null(slope_meta$coefficient$kink_effect)) {
          kink_effect <- as.numeric(slope_meta$coefficient$kink_effect)
        }
        if (!is.null(slope_meta$coefficient$kink_location)) {
          kink_location <- as.numeric(slope_meta$coefficient$kink_location)
        }
      }

      if (is.list(slope_meta$detail)) {
        if (is.na(kink_effect) && !is.null(slope_meta$detail$kink_effect)) {
          kink_effect <- as.numeric(slope_meta$detail$kink_effect)
        }
        if (is.na(kink_location) && !is.null(slope_meta$detail$kink_location)) {
          kink_location <- as.numeric(slope_meta$detail$kink_location)
        }
      }

      if (!is.finite(kink_effect) || !is.finite(kink_location)) {
        return()
      }

      se_term <- pmax(se_curve_grid - kink_location, 0)
      x_pred <- intercept + kink_effect * se_term

      graphics::lines(x_pred, se_curve_grid, lwd = 2, col = "black")

      ci_offset <- 1.96 * intercept_se
      graphics::lines(x_pred - ci_offset, se_curve_grid, lty = 2, col = "black")
      graphics::lines(x_pred + ci_offset, se_curve_grid, lty = 2, col = "black")

      draw_vertical_segment(
        kink_location,
        lty = 3,
        lwd = 1,
        col = grDevices::adjustcolor("black", alpha.f = 0.5)
      )
    }

    if (slope_type %in% c("quadratic") || (slope_meta$quadratic && slope_type == "linear")) {
      draw_quadratic_or_linear(power = 2)
    } else if (slope_type %in% c("linear")) {
      draw_quadratic_or_linear(power = 1)
    } else if (slope_type %in% c("kink", "kinked")) {
      draw_kinked()
    }
  }

  if (length(x_values) > 0) {
    graphics::points(
      x = x_values,
      y = y_values,
      pch = plot_pch,
      col = funnel_opts$col,
      bg = point_bg,
      cex = point_cex
    )
  }

  x_ticks <- pretty(xlim, n = 6)
  x_ticks_int <- sort(unique(round(x_ticks)))
  x_ticks_int <- x_ticks_int[x_ticks_int >= floor(xlim[1]) & x_ticks_int <= ceiling(xlim[2])]
  if (length(x_ticks_int) < 2) {
    x_ticks_int <- round(seq(from = xlim[1], to = xlim[2], length.out = 5))
    x_ticks_int <- sort(unique(x_ticks_int))
  }
  graphics::axis(
    1,
    at = x_ticks_int,
    labels = sprintf("%d", x_ticks_int)
  )

  if (length(y_ticks) > 0) {
    graphics::axis(
      2,
      at = y_ticks,
      labels = maive_funnel_format_axis_labels(y_ticks, digits = funnel_opts$digits),
      las = 1
    )
  }

  graphics::box()

  par_usr <- graphics::par("usr")
  y_span <- abs(par_usr[4] - par_usr[3])
  top_offset <- y_span * 0.04
  label_y_simple <- par_usr[4] - top_offset
  label_y_maive <- label_y_simple

  par_xpd_old <- graphics::par("xpd")
  on.exit(graphics::par(xpd = par_xpd_old), add = TRUE)
  graphics::par(xpd = NA)

  simple_mean_label <- paste0("Simple mean = ", round(simple_mean, 2))

  graphics::text(
    simple_mean,
    label_y_simple,
    labels = simple_mean_label,
    cex = 0.9,
    adj = c(0.5, 0)
  )

  if (!is.null(intercept) && !is.null(intercept_se)) {
    intercept_label <- if (instrument == 0) "Regression fit" else model_label
    maive_label <- paste0(
      intercept_label,
      " = ",
      round(intercept, 2),
      " (SE = ",
      round(intercept_se, 2),
      ")"
    )

    intercept_clamped <- max(xlim[1], min(xlim[2], intercept))

    simple_width <- graphics::strwidth(simple_mean_label, cex = 0.9)
    maive_width <- graphics::strwidth(maive_label, cex = 0.9)

    simple_left <- simple_mean - simple_width / 2
    simple_right <- simple_mean + simple_width / 2
    maive_left <- intercept_clamped - maive_width / 2
    maive_right <- intercept_clamped + maive_width / 2

    labels_overlap <- !(simple_right < maive_left || maive_right < simple_left)

    if (labels_overlap) {
      label_y_maive <- label_y_simple - y_span * 0.05
    }

    graphics::text(
      intercept_clamped,
      label_y_maive,
      labels = maive_label,
      cex = 0.9,
      adj = c(0.5, 0)
    )
  }

  graphics::par(xpd = par_xpd_old)

  p_value_labels <- expression()
  p_legend_fill <- character(0)
  if (!is.null(ci_data) && nrow(ci_data) > 0) {
    level_names <- as.character(ci_data$level)

    if ("90" %in% level_names) {
      p_value_labels <- c(p_value_labels, "p >= 0.10")
      p_legend_fill <- c(p_legend_fill, shade_cols["90"])
    }

    if ("95" %in% level_names) {
      p_value_labels <- c(p_value_labels, "0.10 > p >= 0.05")
      p_legend_fill <- c(p_legend_fill, shade_cols["95"])
    }

    if ("99" %in% level_names) {
      p_value_labels <- c(p_value_labels, "0.05 > p >= 0.01")
      p_legend_fill <- c(p_legend_fill, shade_cols["99"])
    }

    p_value_labels <- c(p_value_labels, "p < 0.01")
    p_legend_fill <- c(p_legend_fill, outer_fill_col)
  }

  fit_label <- if (instrument == 0) "Regression fit" else paste(model_label, "fit")

  legend_labels <- c(funnel_opts$legend_texts[1])
  legend_pch <- c(funnel_opts$effect_pch[1])
  legend_col <- c(funnel_opts$text_color)
  legend_pt_bg <- c(funnel_opts$effect_shades[1])
  legend_pt_cex <- c(base_point_cex)
  legend_lty <- c(NA)
  legend_lwd <- c(NA)

  if (use_adjusted) {
    legend_labels <- c(legend_labels, funnel_opts$legend_texts[2])
    legend_pch <- c(legend_pch, funnel_opts$effect_pch[2])
    legend_col <- c(legend_col, funnel_opts$text_color)
    legend_pt_bg <- c(legend_pt_bg, funnel_opts$effect_shades[2])
    adjusted_legend_cex <- if (!is.null(adjusted_cex) && any(is.finite(adjusted_cex))) {
      stats::median(adjusted_cex[is.finite(adjusted_cex)])
    } else {
      base_point_cex
    }
    legend_pt_cex <- c(legend_pt_cex, adjusted_legend_cex)
    legend_lty <- c(legend_lty, NA)
    legend_lwd <- c(legend_lwd, NA)
  }

  legend_labels <- c(legend_labels, "Simple mean")
  legend_pch <- c(legend_pch, NA)
  legend_col <- c(legend_col, "black")
  legend_pt_bg <- c(legend_pt_bg, NA)
  legend_pt_cex <- c(legend_pt_cex, NA)
  legend_lty <- c(legend_lty, 4)
  legend_lwd <- c(legend_lwd, 2)

  legend_labels <- c(legend_labels, fit_label)
  legend_pch <- c(legend_pch, NA)
  legend_col <- c(legend_col, "black")
  legend_pt_bg <- c(legend_pt_bg, NA)
  legend_pt_cex <- c(legend_pt_cex, NA)
  legend_lty <- c(legend_lty, 1)
  legend_lwd <- c(legend_lwd, 2)

  legend_labels <- c(legend_labels, funnel_opts$legend_texts[4])
  legend_pch <- c(legend_pch, NA)
  legend_col <- c(legend_col, "black")
  legend_pt_bg <- c(legend_pt_bg, NA)
  legend_pt_cex <- c(legend_pt_cex, NA)
  legend_lty <- c(legend_lty, 2)
  legend_lwd <- c(legend_lwd, 1)

  graphics::legend(
    funnel_opts$legend_position,
    legend = legend_labels,
    pch = legend_pch,
    col = legend_col,
    pt.bg = legend_pt_bg,
    pt.cex = legend_pt_cex,
    lty = legend_lty,
    lwd = legend_lwd,
    bg = funnel_opts$legend_bg,
    bty = funnel_opts$legend_bty,
    inset = funnel_opts$legend_inset
  )

  if (length(p_value_labels) > 0) {
    graphics::legend(
      "topright",
      legend = p_value_labels,
      fill = p_legend_fill,
      border = "gray40",
      bg = funnel_opts$legend_bg,
      bty = funnel_opts$legend_bty,
      inset = funnel_opts$legend_inset
    )
  }

  invisible(NULL)
}

#' Draw a funnel plot (base graphics)
#'
#' High-level wrapper that accepts the original input data and the result returned by
#' `maive()` / `waive()`, prepares the necessary vectors/metadata, and draws the funnel plot
#' on the currently active graphics device.
#'
#' Device management (PNG/SVG/PDF, resolution, base64 encoding) is intentionally left
#' to the caller.
#'
#' @param dat Data frame containing at least numeric columns `bs` (effect sizes)
#'   and `sebs` (standard errors).
#' @param result Result list returned by `maive()` or `waive()`.
#' @param instrument Optional indicator (0/1). If `NULL`, inferred from `result`.
#' @param model_type Label used for plot text/legend (e.g., `"MAIVE"` or `"WAIVE"`).
#' @return Invisibly returns `NULL`.
#' @export
get_funnel_plot <- function(dat, result, instrument = NULL, model_type = "MAIVE") {
  maive_funnel_validate_result(result)
  dat_extracted <- maive_funnel_validate_dat(dat)

  effect <- dat_extracted$effect
  se <- dat_extracted$se
  n_points <- length(effect)

  model_type_norm <- maive_funnel_normalize_model_type(model_type)
  instrument_norm <- maive_funnel_normalize_instrument(instrument, result)

  se_adjusted <- maive_funnel_extract_adjusted_se(result, instrument_norm, n_points)
  adjusted_weights <- maive_funnel_extract_adjusted_weights(result, model_type_norm, instrument_norm, n_points)
  slope <- maive_funnel_extract_slope(result)

  maive_funnel_plot_engine(
    effect = effect,
    se = se,
    se_adjusted = se_adjusted,
    adjusted_weights = adjusted_weights,
    intercept = result[["beta"]],
    intercept_se = result[["SE"]],
    slope = slope,
    instrument = instrument_norm,
    model_type = model_type_norm
  )
}
