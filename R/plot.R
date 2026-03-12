#' Plot Changepoint Detection Results
#'
#' Visualizes detected changepoints overlaid on the original time series.
#' Supports a series view with colored segments, a diagnostics view showing
#' segment means as a step function, or both panels stacked.
#'
#' @details
#' Three plot types are available:
#'
#' **`type = "series"` (default).** The original time series is drawn as a
#' line with segments colored by their segment ID. Vertical dashed lines
#' mark changepoint locations. This view immediately shows where the series
#' was split and how distinct the segments are.
#'
#' **`type = "diagnostics"`.** The segment means are drawn as a step
#' function (horizontal lines at each segment's mean), with vertical dashed
#' lines at changepoint locations and points marking the mean of each
#' segment. This view focuses on the magnitude and direction of each
#' detected shift.
#'
#' **`type = "both"`.** Stacks the series panel on top and the diagnostics
#' panel below using [patchwork::wrap_plots()], with shared time axes.
#'
#' @export
#' @param x \[`changepoint`]\cr
#'   An object of class `changepoint` as returned by
#'   [detect_cpts()].
#' @param type \[`character(1)`: `"series"`]\cr
#'   Plot type: `"series"` for colored segments with changepoint lines,
#'   `"diagnostics"` for segment mean step function, or `"both"` for
#'   both panels stacked.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object (or a [patchwork::wrap_plots()]
#'   composite when `type = "both"`).
#'
#' @seealso [detect_cpts()] for computing the changepoint analysis.
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- c(rnorm(200, 0, 1), rnorm(200, 5, 1))
#' cpts <- detect_cpts(x, method = "pelt")
#' plot(cpts)
#' plot(cpts, type = "diagnostics")
#' plot(cpts, type = "both")
#' }
plot.changepoint <- function(x, type = "series", ...) {
  check_missing(x)
  check_class(x, "changepoint")
  type <- check_match(type, c("series", "diagnostics", "both"))
  cpt_locs <- attr(x, "changepoint_locations")
  n_cpt <- attr(x, "n_changepoints")
  colors <- changepoint_colors_()
  observed_states <- levels(x$state)
  # Map state names to colors by name; fall back to cycling for unknowns
  seg_colors <- vapply(
    observed_states,
    function(s) {
      if (s %in% names(colors)) {
        colors[[s]]
      } else {
        colors[[(match(s, observed_states) - 1L) %% length(colors) + 1L]]
      }
    },
    character(1L)
  )
  names(seg_colors) <- observed_states
  if (type == "series" || type == "both") {
    p_series <- x |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = !!rlang::sym("time"),
          y = !!rlang::sym("value"),
          color = !!rlang::sym("state")
        )
      ) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::scale_color_manual(
        values = seg_colors,
        name = "State"
      )
    if (n_cpt > 0L) {
      cpt_df <- data.frame(xint = x$time[cpt_locs + 1L])
      p_series <- p_series +
        ggplot2::geom_vline(
          data = cpt_df,
          mapping = ggplot2::aes(xintercept = !!rlang::sym("xint")),
          linetype = "dashed",
          color = "red",
          linewidth = 0.6,
          inherit.aes = FALSE
        )
    }
    p_series <- p_series +
      ggplot2::labs(
        title = sprintf(
          "Changepoint Detection (%s): %d changepoint%s",
          attr(x, "method"),
          n_cpt,
          ifelse_(n_cpt == 1L, "", "s")
        ),
        x = "Time",
        y = "Value"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(color = "black", face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )
    if (type == "series") {
      return(p_series)
    }
  }
  if (type == "diagnostics" || type == "both") {
    seg_ids <- sort(unique(x$segment))
    step_df <- do.call(
      base::rbind,
      lapply(
        seg_ids,
        function(s) {
          seg_rows <- x[x$segment == s, ]
          data.frame(
            xmin = min(seg_rows$time),
            xmax = max(seg_rows$time),
            ymean = seg_rows$segment_mean[1L],
            state = seg_rows$state[1L],
            stringsAsFactors = FALSE
          )
        }
      )
    )
    step_df$xmid <- (step_df$xmin + step_df$xmax) / 2
    p_diag <- ggplot2::ggplot(step_df) +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = !!rlang::sym("xmin"),
          xend = !!rlang::sym("xmax"),
          y = !!rlang::sym("ymean"),
          yend = !!rlang::sym("ymean"),
          color = !!rlang::sym("state")
        ),
        linewidth = 1.2
      ) +
      ggplot2::geom_point(
        ggplot2::aes(
          x = !!rlang::sym("xmid"),
          y = !!rlang::sym("ymean"),
          color = !!rlang::sym("state")
        ),
        size = 3
      ) +
      ggplot2::scale_color_manual(
        values = seg_colors,
        name   = "State"
      )
    if (n_cpt > 0L) {
      cpt_df <- data.frame(xint = x$time[cpt_locs + 1L])
      p_diag <- p_diag +
        ggplot2::geom_vline(
          data = cpt_df,
          mapping = ggplot2::aes(xintercept = !!rlang::sym("xint")),
          linetype = "dashed",
          color = "red",
          linewidth = 0.6,
          inherit.aes = FALSE
        )
    }
    p_diag <- p_diag +
      ggplot2::labs(
        title = "Segment Means (Step Function)",
        x = "Time",
        y = "Segment Mean"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        axis.title = ggplot2::element_text(color = "black", face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )
    if (type == "diagnostics") {
      return(p_diag)
    }
  }
  p_series <- p_series +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
  patchwork::wrap_plots(p_series, p_diag, ncol = 1L, heights = c(1, 1))
}

#' Plot EWS Results
#'
#' @export
#' @param x \[`ews`]\cr
#'   Output of [detect_warnings()].
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @examples
#' set.seed(123)
#' ts_data <- stats::arima.sim(list(order = c(1, 1, 0), ar = 0.6), n = 200)
#' ews_roll <- detect_warnings(ts_data)
#' plot(ews_roll)
#'
plot.ews <- function(x, ...) {
  check_missing(x)
  check_class(x, "ews")
  d <- data.frame(
    value = attr(x, "orig_values"),
    time = attr(x, "orig_time")
  )
  p_ts <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("value"))
  ) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  if (attr(x, "method") == "rolling") {
    p_ts <- p_ts +
      ggplot2::labs(title = "Time Series", y = "Value", x = NULL)
    p_metrics <- plot_rolling_ews(x)
    patchwork::wrap_plots(p_ts, p_metrics, ncol = 1L)
  } else {
    warn <- data.frame(time = x$time[x$detected == 1])
    p_ts <- p_ts +
      ggplot2::geom_rug(
        data = warn,
        ggplot2::aes(
          x = !!rlang::sym("time"), color = "EWS"
        ),
        show.legend = TRUE,
        sides = "b",
        length = ggplot2::unit(0.05, "npc"),
        inherit.aes = FALSE
      ) +
      # Add point to control shape
      ggplot2::geom_point(
        data = data.frame(x = -Inf, y = -Inf),
        ggplot2::aes(
          x = !!rlang::sym("x"), y = !!rlang::sym("y"), color = "Void"
        ),
        show.legend = TRUE,
        na.rm = TRUE
      ) +
      ggplot2::labs(
        title = "Time Series with Detected Warnings",
        y = "Value",
        x = NULL
      ) +
      ggplot2::scale_color_manual(
        name = "EWS",
        values = c(EWS = "red", Void = "white"),
        labels = c("Detected", "")
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(
          override.aes = list(
            linetype = c(0, 0),
            shape = c("|", "."),
            size = 4
          )
        )
      )
    p_metrics <- plot_expanding_ews(x)
    p_cls <- plot_classification(attr(x, "classification"))
    patchwork::wrap_plots(
      p_ts, p_metrics, p_cls, ncol = 1L, heights = c(4, 8, 1)
    )
  }
}

plot_rolling_ews <- function(x) {
  cor <- attr(x, "cor")
  x$metric <- factor(
    x$metric,
    levels = names(cor),
    labels = lapply(
      paste0(
        names(cor), "~(tau == ", round(cor, 2), ")"
      ),
      str2lang
    )
  )
  ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("std"),
      color = !!rlang::sym("metric")
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::facet_wrap(
      stats::as.formula("~metric"),
      scales = "free_y",
      ncol = 3,
      drop = TRUE,
      labeller = ggplot2::label_parsed
    ) +
    ggplot2::labs(
      title = "Rolling EWS Indicators",
      y = "Metric Value",
      x = "Time"
    ) +
    ggplot2::scale_color_viridis_d() +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      strip.text = ggplot2::element_text(
        face = "bold", size = 10, margin = ggplot2::margin(b = 5)
      ),
      panel.spacing = ggplot2::unit(1, "lines"),
      axis.text = ggplot2::element_text(size = 8),
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(angle = 90, size = 10),
      plot.title = ggplot2::element_text(
        margin = ggplot2::margin(b = 5)
      ),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
    )
}

plot_expanding_ews <- function(x) {
  warn <- x[x$detected == 1, ]
  state_colors <- c(
    "Stable" = "#440154FF",
    "Vulnerable" = "#3B528BFF",
    "Weak Warning" = "#21908CFF",
    "Strong Warning" = "#5DC863FF",
    "Failing" = "#FDE725FF",
    "Warning" = "orange",
    "Critical" = "red"
  )
  ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("z_score"),
      color = !!rlang::sym("metric"),
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point(
      data = warn,
      ggplot2::aes(
        x = !!rlang::sym("time"),
        y = !!rlang::sym("z_score"),
        alpha = "EWS"
      ),
      size = 2.5
    ) +
    # Add line for alpha
    ggplot2::geom_line(
      data = warn,
      ggplot2::aes(
        x = !!rlang::sym("time"),
        y = !!rlang::sym("z_score"),
        color = !!rlang::sym("metric"),
        alpha = "Void"
      ),
      inherit.aes = FALSE,
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(-4, -2, 0, 2, 4)
    ) +
    ggplot2::labs(
      title = "Strength of Early Warning Signals",
      y = "Scaled Metric Value",
      x = NULL
    ) +
    ggplot2::scale_color_viridis_d(
      name = "EWS Indicator",
      guide = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          linewidth = 1,
          linetype = "solid",
          shape = NA
        )
      )
    ) +
    ggplot2::scale_alpha_manual(
      name = "EWS",
      values = c(EWS = 1, Void = 1),
      labels = c("Detected", "Not Detected"),
      guide = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          alpha = c(1, 1),
          linetype = c("blank", "solid")
        )
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::geom_hline(
      yintercept = c(1, -1) * attr(x, "threshold"),
      linetype = "dashed",
      color = "grey50"
    )
}

plot_classification <- function(x) {
  state_colors <- c(
    `Stable` = "#440154FF",
    `Vulnerable` = "#3B528BFF",
    `Warning` = "#5DC863FF",
    `Critical` = "#FDE725FF",
    `Failing` = "orange"
  )
  x$state <- factor(x$state, levels = names(state_colors))
  ggplot2::ggplot(x, ggplot2::aes(x = !!rlang::sym("time"), y = 1)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = !!rlang::sym("state"), height = 1),
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = state_colors,
      name = "System\n  State",
      drop = FALSE
    ) +
    ggplot2::labs(x = NULL, y = "") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

#' Plot Time Series Data with Detected Regime Stability
#'
#' @export
#' @param x \[`regimes`]\cr Output of [detect_regimes()].
#' @param points \[`logical(1)`]\cr Should a point be added for each
#'   observation?  The points are colored by regime stability
#'   (default: `FALSE`).
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @examples
#' set.seed(123)
#' ts_data <- stats::arima.sim(list(order = c(1, 1, 0), ar = 0.6), n = 200)
#' regimes <- detect_regimes(
#'   data = ts_data,
#'   method = "threshold",
#'   sensitivity = "medium"
#' )
#' plot(regimes)
#'
plot.regimes <- function(x, points = FALSE, ...) {
  check_missing(x)
  check_class(x, "regimes")
  check_flag(points)
  states <- factor(base::sort(dplyr::pull(x[, "stability", drop = FALSE], 1L)))
  data <- x |>
    dplyr::select(
      !!rlang::sym("value"),
      !!rlang::sym("time"),
      !!rlang::sym("stability")
    ) |>
    dplyr::rename(state = !!rlang::sym("stability"))
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("value"))
  )
  xmin <- rlang::sym(".min")
  xmax <- rlang::sym(".max")
  ymin <- rlang::sym(".neginf")
  ymax <- rlang::sym(".posinf")
  rects <- data |>
    dplyr::arrange(!!rlang::sym("time")) |>
    dplyr::mutate(
      .grouping_var = cumsum(
        !!rlang::sym("state") != dplyr::lag(
          !!rlang::sym("state"),
          default = dplyr::first(!!rlang::sym("state"))
        )
      ),
      .lag = dplyr::lag(
        !!rlang::sym("time"),
        default = dplyr::first(!!rlang::sym("time"))
      ),
      .lead = dplyr::lead(
        !!rlang::sym("time"),
        default = dplyr::last(!!rlang::sym("time"))
      )
    ) |>
    dplyr::group_by(
      !!rlang::sym(".grouping_var"), !!rlang::sym("state")
    ) |>
    dplyr::summarise(
      .neginf = -Inf,
      .posinf = Inf,
      .min = 0.5 * (min(!!rlang::sym("time")) + min(!!rlang::sym(".lag"))),
      .max = 0.5 * (max(!!rlang::sym("time")) + max(!!rlang::sym(".lead"))),
      .groups = "drop"
    )
  p <- p + ggplot2::geom_rect(
    data = rects,
    ggplot2::aes(
      xmin = !!xmin,
      xmax = !!xmax,
      ymin = !!ymin,
      ymax = !!ymax,
      fill = !!rlang::sym("state")
    ),
    alpha = 0.5,
    show.legend = TRUE,
    inherit.aes = FALSE
  ) +
    ggplot2::geom_line(linewidth = .5)
  if (points) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(fill = !!rlang::sym("state")),
      show.legend = FALSE,
      pch = 21
    )
  }
  p +
    ggplot2::scale_fill_brewer(
      palette = ifelse(
        n_unique(states) <= 8,
        "Accent",
        "Set3"
      ),
      limits = levels(states),
      name = "State",
      drop = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time", y = "Value") +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot Discovered Patterns
#'
#' @export
#' @param x \[`patterns`]\cr
#'   Output of [discover_patterns()].
#' @param n \[`integer(1)`: `10L`]\cr
#'   Maximum number of patterns to include in the plot.
#' @param prop \[`logical(1)`: `TRUE`]\cr
#'   Should outcome-specific count proportions be displayed in the plot?
#'   The default is `TRUE`. Ignored if `outcome` was not originally specified.
#' @param group \[`character(1)`]\cr
#'   Name of the outcome class to draw the plot for. If not provided, shows the
#'   counts and proportions by outcome class for each pattern. If provided,
#'   only the proportions of the specific class are drawn by pattern.
#'   Ignored if `outcome` was not originally specified.
#' @param global \[`logical(1)`: `TRUE`]\cr
#'   Should a line be added showing the global proportion when `group` is
#'   provided? Also colors the patterns according to whether the proportion is
#'   above or below the global value.
#' @param ... Ignored.
#' @return A `ggplot` object.
#' ngrams <- discover_patterns(engagement, type = "ngram")
#' plot(ngrams)
#'
plot.patterns <- function(x, n = 10L, prop = TRUE, group, global = TRUE, ...) {
  ifelse_(
    missing(group) || is.null(attr(x, "group")),
    plot_patterns_all(x = x, n = n, prop = prop),
    plot_patterns_group(x = x, n = n, group = group, global = global)
  )
}

plot_patterns_all <- function(x, n, prop) {
  p <- NULL
  x <- x |>
    dplyr::arrange(dplyr::desc(!!rlang::sym("count"))) |>
    dplyr::slice_head(n = n)
  if (is.null(attr(x, "group")) || !prop) {
    p <- x |>
      ggplot2::ggplot(
        ggplot2::aes(
          y = stats::reorder(!!rlang::sym("pattern"), !!rlang::sym("count")),
          x = !!rlang::sym("count")
        )
      ) +
      ggplot2::geom_col()
  } else {
    groups <- unique(attr(x, "group"))
    count_cols <- paste0("count_", groups)
    x <- x |>
      dplyr::select(tidyselect::all_of(c("pattern", count_cols))) |>
      tidyr::pivot_longer(
        cols = tidyselect::all_of(count_cols),
        names_to = "group",
        names_prefix = "count_",
        values_to = "count"
      ) |>
      dplyr::mutate(
        group = factor(!!rlang::sym("group"), levels = rev(groups))
      ) |>
      dplyr::group_by(!!rlang::sym("pattern")) |>
      dplyr::mutate(
        prop = !!rlang::sym("count") / sum(!!rlang::sym("count"))
      )
    p <- x |>
      ggplot2::ggplot(
        ggplot2::aes(
          y = stats::reorder(!!rlang::sym("pattern"), !!rlang::sym("count")),
          x = !!rlang::sym("count"),
          fill = !!rlang::sym("group")
        )
      ) +
      ggplot2::geom_col() +
      ggplot2::geom_text(
        ggplot2::aes(
          label = scales::label_percent(
            accuracy = 0.1,
            suffix = " %"
          )(!!rlang::sym("prop"))
        ),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "white",
        size = 4
      ) +
      ggplot2::scale_fill_discrete(
        name = NULL,
        guide = ggplot2::guide_legend(reverse = TRUE)
      )
  }
  p +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Count", y = "")
}

plot_patterns_group <- function(x, n, group, global) {
  count_col <- paste0("count_", group)
  global_prop <- mean(attr(x, "group") == group)
  x |>
    dplyr::select(tidyselect::all_of(c("pattern", "count", count_col))) |>
    dplyr::mutate(
      prop = !!rlang::sym(count_col) / !!rlang::sym("count"),
      group = factor(
        !!rlang::sym("prop") > global_prop,
        levels = c(TRUE, FALSE)
      )
    ) |>
    dplyr::arrange(dplyr::desc(!!rlang::sym("prop"))) |>
    dplyr::slice_head(n = n) |>
    ggplot2::ggplot(
      ggplot2::aes(
        y = stats::reorder(!!rlang::sym("pattern"), !!rlang::sym("prop")),
        x = !!rlang::sym("prop"),
        fill = onlyif(global, !!rlang::sym("group"))
      )
    ) +
    ggplot2::geom_col() +
    onlyif(
      global,
      ggplot2::geom_vline(
        xintercept = global_prop,
        color = "cyan",
        linetype = "dashed",
        lwd = 1
      )
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = scales::label_percent(
          accuracy = 0.1,
          suffix = " %"
        )(!!rlang::sym("prop"))
      ),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white",
      size = 4
    ) +
    ggplot2::scale_fill_manual(values = c("darkgreen", "gray25")) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Proportion", y = "")
}

#' Plot Hurst Analysis Results
#'
#' Visualizes rolling Hurst exponent results with state-colored backgrounds
#' and threshold bands. Supports displaying the original time series, the
#' Hurst trajectory, or both panels stacked.
#'
#' @details
#' Three plot types are available:
#'
#' **`"series"`.**
#' Draws the original time series as a line plot with the background
#' shaded according to the Hurst state classification at each time point.
#' The five states (`strong_antipersistent`, `antipersistent`,
#' `random_walk`, `persistent`, `strong_persistent`) are mapped to a
#' fixed color palette. This view answers the question: "What was the
#' system doing when it was in each long-range dependence regime?"
#'
#' **`"states"`.**
#' Plots the rolling Hurst exponent trajectory as a line, overlaid on
#' horizontal bands at the classification thresholds (0.4, 0.5, 0.6,
#' 0.7). A dashed reference line at \eqn{H = 0.5} marks the boundary
#' between antipersistent and persistent behavior. This view answers:
#' "How did the Hurst exponent evolve, and how close is it to regime
#' boundaries?"
#'
#' **`"both"`** (default).
#' Stacks the series panel on top and the states panel below using
#' [patchwork::wrap_plots()], with shared time axes. The top panel
#' suppresses its x-axis labels to avoid duplication.
#'
#' @export
#' @param x \[`hurst`]\cr
#'   An object of class `hurst`.
#' @param type \[`character(1)`: `"both"`]\cr
#'   Plot type. The available options are:
#'
#'   * `"series"`: Time series with state-colored background.
#'   * `"states"`: Hurst exponent trajectory with threshold bands.
#'   * `"both"`: Both panels stacked vertically.
#'
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object (or a
#'   [patchwork::wrap_plots()] composite when `type = "both"`).
#' @references
#' Peng, C.-K., Buldyrev, S.V., Havlin, S., Simons, M., Stanley, H.E.,
#' & Goldberger, A.L. (1994). Mosaic organization of DNA nucleotides.
#' \emph{Physical Review E}, 49(2), 1685--1689.
#' \doi{10.1103/PhysRevE.49.1685}
#'
#' @seealso [hurst()] for computing the rolling Hurst exponent.
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' h <- hurst(x, window = 50L)
#' plot(h, type = "series")
#' plot(h, type = "states")
#' plot(h, type = "both")
#' }
plot.hurst <- function(x, type = "both", ...) {
  check_missing(x)
  check_class(x, "hurst")
  type <- check_match(type, c("series", "states", "both"))
  colors <- hurst_state_colors()
  all_states <- c(
    "strong_antipersistent", "antipersistent", "random_walk",
    "persistent", "strong_persistent"
  )
  d <- dplyr::mutate(
    x,
    state_f = factor(
      !!rlang::sym("state"),
      levels = all_states
    )
  )
  # Build state-colored background rectangles
  rects <- build_state_rects_(d)
  if (type == "series" || type == "both") {
    val <- d[["value"]]
    val <- val[is.finite(val)]
    y_rng <- range(val, na.rm = TRUE)
    y_pad <- diff(y_rng) * 0.08
    y_lo  <- y_rng[1L] - y_pad
    y_hi  <- y_rng[2L] + y_pad
    rects[["ymin"]] <- y_lo
    rects[["ymax"]] <- y_hi
    p_series <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("value"))
    ) +
      ggplot2::geom_rect(
        data = rects,
        ggplot2::aes(
          xmin = !!rlang::sym("xmin"),
          xmax = !!rlang::sym("xmax"),
          ymin = !!rlang::sym("ymin"),
          ymax = !!rlang::sym("ymax"),
          fill = !!rlang::sym("state_f")
        ),
        alpha = 0.25,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_line(linewidth = 0.4) +
      ggplot2::scale_fill_manual(
        values = colors,
        limits = all_states,
        name = "State",
        drop = FALSE
      ) +
      ggplot2::coord_cartesian(
        ylim = c(y_lo, y_hi),
        clip = "off"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL, y = "Value") +
      ggplot2::theme(
        legend.position = "bottom",
        plot.margin = ggplot2::margin(t = 12, r = 5, b = 5, l = 5)
      )
    if (type == "series") {
      return(p_series)
    }
  }
  if (type == "states" || type == "both") {
    h_vals <- d[["hurst"]]
    h_vals <- h_vals[is.finite(h_vals)]
    h_max  <- max(h_vals, 1.0, na.rm = TRUE)
    h_top  <- h_max + 0.05
    p_states <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("hurst"))
    ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0, ymax = 0.4,
        fill = colors["strong_antipersistent"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.4, ymax = 0.5,
        fill = colors["antipersistent"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.5, ymax = 0.6,
        fill = colors["random_walk"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.6, ymax = 0.7,
        fill = colors["persistent"], alpha = 0.15
      ) +
      ggplot2::annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0.7, ymax = h_top,
        fill = colors["strong_persistent"], alpha = 0.15
      ) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::geom_hline(
        yintercept = 0.5,
        linetype = "dashed",
        color = "grey40"
      ) +
      ggplot2::scale_y_continuous(
        breaks = seq(0, ceiling(h_top * 4) / 4, 0.25)
      ) +
      ggplot2::coord_cartesian(
        ylim = c(-0.02, h_top),
        clip = "off"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Time", y = "Hurst Exponent") +
      ggplot2::theme(
        plot.margin = ggplot2::margin(t = 8, r = 5, b = 5, l = 5)
      )
    if (type == "states") {
      return(p_states)
    }
  }
  p_series <- p_series +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  patchwork::wrap_plots(p_series, p_states, ncol = 1L, heights = c(1, 1))
}

#' Plot Hurst Early Warning Signals
#'
#' Creates a three-panel stacked visualization of early warning signals
#' derived from rolling Hurst exponent analysis.
#'
#' @details
#' The plot consists of three vertically stacked panels:
#'
#' **Panel 1 -- Hurst trajectory with warning backgrounds.**
#' The rolling Hurst exponent is drawn as a line. Behind it,
#' colored rectangles shade each time region according to its warning
#' level (green = none, yellow = moderate, orange = high, red =
#' critical). A dashed line at \eqn{H = 0.5} marks the random-walk
#' boundary. This panel links changes in long-range dependence to the
#' composite warning score.
#'
#' **Panel 2 -- Warning level step plot.**
#' A step function shows the discrete warning level (0--4) over time.
#' Sudden jumps reveal the onset and duration of elevated warning
#' periods. Reading this jointly with Panel 1 shows whether warning
#' escalation coincides with Hurst crossing a regime boundary.
#'
#' **Panel 3 -- Indicator heatmap.**
#' Each of the selected binary indicators is displayed as a row, with
#' time on the x-axis. Active indicators (value 1) are shaded with a
#' warm color gradient; inactive indicators are light grey. This
#' reveals which specific mechanisms (e.g., flickering, trend, spectral
#' shift) contribute to each warning, enabling targeted interpretation
#' rather than relying solely on the composite score.
#'
#' The `indicators` argument controls which rows appear in Panel 3.
#' Subsetting to a few indicators of interest can simplify the display
#' for presentation.
#'
#' @export
#' @param x \[`hurst_ews`]\cr
#'   An object of class `hurst_ews` as returned by
#'   [detect_hurst_warnings()].
#' @param indicators \[`character()`]\cr
#'   Which indicators to display in the heatmap panel. Defaults to all
#'   10 indicators.
#' @param ... Additional arguments (currently unused).
#' @return A [patchwork::wrap_plots()] composite of three
#'   [ggplot2::ggplot()] panels.
#'
#' @references
#' Scheffer, M., Bascompte, J., Brock, W.A., Brovkin, V., Carpenter,
#' S.R., Dakos, V., Held, H., van Nes, E.H., Rietkerk, M., &
#' Sugihara, G. (2009). Early-warning signals for critical transitions.
#' \emph{Nature}, 461, 53--59. \doi{10.1038/nature08227}
#'
#' @seealso [detect_hurst_warnings()] for computing the warning signals
#'   displayed by this plot.
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' h <- hurst(x, window = 50L)
#' ews <- detect_hurst_warnings(h)
#'
#' # Default: all indicators
#' plot(ews)
#'
#' # Subset to specific indicators
#' plot(ews, indicators = c("flickering", "trend_up", "spectral_shift"))
#' }
plot.hurst_ews <- function(x, indicators = NULL, ...) {
  check_missing(x)
  check_class(x, "hurst_ews")
  all_indicators <- c(
    "extreme_low", "extreme_high", "trend_up", "trend_down",
    "high_volatility", "flickering", "variance_ratio",
    "spectral_shift", "autocorr_increase", "state_persistence"
  )
  if (is.null(indicators)) indicators <- all_indicators
  indicators <- check_match(
    indicators, all_indicators, several.ok = TRUE
  )
  warning_colors <- c(
    "0" = "#2ca02c",
    "1" = "#98df8a",
    "2" = "#ffcc00",
    "3" = "#ff8c00",
    "4" = "#d62728"
  )
  # Panel 1: Hurst line with warning background
  warn_rects <- build_warning_rects_(x)
  p_hurst <- ggplot2::ggplot(
    x,
    ggplot2::aes(x = !!rlang::sym("time"), y = !!rlang::sym("hurst"))
  )
  if (nrow(warn_rects) > 0L) {
    p_hurst <- p_hurst +
      ggplot2::geom_rect(
        data = warn_rects,
        ggplot2::aes(
          xmin = !!rlang::sym("xmin"),
          xmax = !!rlang::sym("xmax"),
          ymin = -Inf, ymax = Inf,
          fill = !!rlang::sym("level")
        ),
        alpha = 0.3,
        inherit.aes = FALSE
      )
  }
  p_hurst <- p_hurst +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_hline(
      yintercept = 0.5, linetype = "dashed", color = "grey40"
    ) +
    ggplot2::scale_fill_manual(
      values = warning_colors,
      name = "Warning",
      labels = c("None", "Low", "Moderate", "High", "Critical"),
      drop = FALSE
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "Hurst") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  # Panel 2: Warning level step plot
  p_level <- ggplot2::ggplot(
    x,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("warning_level")
    )
  ) +
    ggplot2::geom_step(linewidth = 0.5, color = "#d62728") +
    ggplot2::scale_y_continuous(
      breaks = 0:4,
      labels = c("None", "Low", "Mod", "High", "Crit"),
      limits = c(-0.2, 4.2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = NULL, y = "Warning Level") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  # Panel 3: Indicator heatmap
  indicator_display <- c(
    extreme_low       = "Extreme Low",
    extreme_high      = "Extreme High",
    trend_up          = "Trend Up",
    trend_down        = "Trend Down",
    high_volatility   = "Volatility",
    flickering        = "Flickering",
    variance_ratio    = "Variance Ratio",
    spectral_shift    = "Spectral Shift",
    autocorr_increase = "Autocorrelation",
    state_persistence = "Persistence"
  )
  long <- tidyr::pivot_longer(
    dplyr::select(x, !!rlang::sym("time"),
                  tidyselect::all_of(indicators)),
    cols = tidyselect::all_of(indicators),
    names_to = "indicator",
    values_to = "intensity"
  )
  long <- dplyr::mutate(
    long,
    indicator_label = factor(
      indicator_display[!!rlang::sym("indicator")],
      levels = rev(indicator_display[indicators])
    )
  )
  p_heat <- ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("indicator_label"),
      fill = !!rlang::sym("intensity")
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colours = c("#f0f0f0", "#ffcc00", "#ff8c00", "#d62728"),
      limits = c(0, 1),
      na.value = "#f0f0f0",
      name = "Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time", y = NULL) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid = ggplot2::element_blank()
    )
  patchwork::wrap_plots(
    p_hurst, p_level, p_heat,
    ncol = 1L,
    heights = c(2, 1, 3)
  )
}

#' @noRd
build_state_rects_ <- function(d) {
  n <- nrow(d)
  if (n == 0L) {
    return(
      data.frame(
        xmin = numeric(0),
        xmax = numeric(0),
        state_f = factor(character(0))
      )
    )
  }
  states <- d$state
  time <- d$time
  # Find run-length encoding of states
  rle_states <- rle(states)
  n_runs <- length(rle_states$lengths)
  rects <- data.frame(
    xmin = numeric(n_runs),
    xmax = numeric(n_runs),
    state_f = factor(
      rle_states$values,
      levels = c(
        "strong_antipersistent", "antipersistent", "random_walk",
        "persistent", "strong_persistent"
      )
    )
  )
  pos <- 1L
  for (i in seq_len(n_runs)) {
    run_start <- pos
    run_end <- pos + rle_states$lengths[i] - 1L
    rects$xmin[i] <- time[run_start]
    rects$xmax[i] <- time[run_end]
    pos <- run_end + 1L
  }
  rects
}

#' @noRd
build_warning_rects_ <- function(x) {
  n <- nrow(x)
  if (n == 0L) {
    return(
      data.frame(
        xmin = numeric(0),
        xmax = numeric(0),
        level = factor(character(0))
      )
    )
  }
  levels_vec <- x$warning_level
  time <- x$time
  rle_levels <- rle(levels_vec)
  n_runs <- length(rle_levels$lengths)
  rects <- data.frame(
    xmin = numeric(n_runs),
    xmax = numeric(n_runs),
    level = factor(
      as.character(rle_levels$values),
      levels = as.character(0:4)
    )
  )
  pos <- 1L
  for (i in seq_len(n_runs)) {
    run_start <- pos
    run_end <- pos + rle_levels$lengths[i] - 1L
    rects$xmin[i] <- time[run_start]
    rects$xmax[i] <- time[run_end]
    pos <- run_end + 1L
  }
  rects
}

#' Plot Multivariate EWS Results
#'
#' Creates a multi-panel visualization of multivariate early warning signal
#' analysis results. For expanding windows, shows the standardized metric
#' strengths with threshold lines and a system-state classification ribbon.
#' For rolling windows, shows faceted metric trends with Kendall's tau
#' annotations. Both methods optionally include a dimension-reduction panel
#' showing MAF1 and PC1 trajectories.
#'
#' @details
#' **Expanding window layout** (top to bottom):
#'
#'   1. Dimension reduction panel (MAF1 + PC1) with warning points.
#'   2. Standardized metric strengths with threshold line(s); points
#'      mark individual metric threshold crossings.
#'   3. System-state classification ribbon (Stable through Failing).
#'
#' **Rolling window layout**:
#'
#' 1. Dimension reduction panel (MAF1 + PC1) with warning points.
#' 2. Faceted panels per metric showing standardized EWS values;
#'    panel titles include Kendall's tau trend statistic.
#'
#' @export
#' @param x \[`multi_ews`]\cr
#'   An object of class `"multi_ews"` produced by
#'   [detect_multivariate_warnings()].
#' @param include_dr \[`logical(1)`: `TRUE`]\cr
#'   Whether to include the dimension reduction (MAF1/PC1) panel.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] or [patchwork::wrap_plots()] object.
#' @seealso [detect_multivariate_warnings()] for computing the EWS.
#' @family multivariate EWS
#' @concept early warning signals
#' @examples
#' \donttest{
#' tip <- generate_tipping_data(n_time = 100, n_vars = 3, tipping_point = 60)
#' ews <- detect_multivariate_warnings(tip, method = "expanding", window = 50)
#' plot(ews)
#' plot(ews, include_dr = FALSE)
#' }
plot.multi_ews <- function(x, include_dr = TRUE, ...) {
  check_missing(x)
  check_class(x, "multi_ews")
  check_flag(include_dr)
  method <- attr(x, "method")
  ifelse_(
    method == "expanding",
    mews_plot_expanding_(x, include_dr),
    mews_plot_rolling_(x, include_dr)
  )
}

#' @noRd
mews_theme <- function(base_size = 11L) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        colour = "grey90",
        linetype = "dashed"
      ),
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black"),
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        colour = "grey95",
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(
        5, 0, 5, 0)
      ),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = ggplot2::rel(1.2),
        hjust = 0
      ),
      plot.subtitle = ggplot2::element_text(
        color = "grey40",
        size = ggplot2::rel(0.9)
      ),
      legend.position = "left",
      legend.key = ggplot2::element_rect(fill = "white", colour = NA)
    )
}


mews_state_colors <- function() {
  c(
    "Stable"     = "#440154FF",
    "Vulnerable" = "#3B528BFF",
    "Warning"    = "#21908CFF",
    "Critical"   = "#5DC863FF",
    "Failing"    = "#FDE725FF"
  )
}

mews_plot_expanding_ <- function(ews_object, include_dr) {
  classification_data <- attr(ews_object, "classification")
  state_colors <- mews_state_colors()
  if (!is.null(classification_data)) {
    classification_data$state <- factor(
      classification_data$state, levels = names(state_colors)
    )
  }
  # Main metric strength plot
  p <- ggplot2::ggplot(
    ews_object,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("z_score"),
      color = !!rlang::sym("metric")
    )
  ) +
    ggplot2::geom_line(alpha = 0.8) +
    ggplot2::geom_point(
      data  = function(d) d[d$detected == 1L, ],
      ggplot2::aes(color = !!rlang::sym("metric")),
      size  = 3.5,
      shape = 19,
      alpha = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = attr(ews_object, "threshold"),
      linetype = "dashed",
      color = "grey50"
    )
  if (attr(ews_object, "tail_direction") == "two.tailed") {
    p <- p + ggplot2::geom_hline(
      yintercept = -attr(ews_object, "threshold"),
      linetype = "dashed",
      color = "grey50"
    )
  }
  p <- p +
    ggplot2::labs(
      title = "Multivariate Early Warning Signals (Expanding Window)",
      subtitle = "Points indicate metric strength crossing the threshold.",
      y = "Strength of EWS",
      x = NULL,
      color = "EWS Indicator"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(shape = 19, size = 3)
      )
    ) +
    mews_theme() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  # Dimension reduction panel
  if (include_dr && !is.null(attr(ews_object, "dimension_reduction"))) {
    p_dr <- mews_plot_dr_(ews_object)
    if (!is.null(classification_data)) {
      p_ribbon <- mews_plot_ribbon_(classification_data, state_colors)
      patchwork::wrap_plots(
        p_dr, p, p_ribbon,
        ncol = 1L,
        heights = c(2, 4, 1)
      )
    } else {
      patchwork::wrap_plots(p_dr, p, ncol = 1L, heights = c(2, 4))
    }
  } else {
    if (!is.null(classification_data)) {
      p_ribbon <- mews_plot_ribbon_(classification_data, state_colors)
      patchwork::wrap_plots(p, p_ribbon, ncol = 1L, heights = c(4, 1))
    } else {
      p + ggplot2::labs(x = "Time Point")
    }
  }
}

mews_plot_rolling_ <- function(ews_object, include_dr) {
  # Build correlation labels
  cor_vals <- attr(ews_object, "cor")
  cor_names <- sub("\\.tau$", "", names(cor_vals))
  cor_data <- data.frame(
    metric = cor_names,
    tau = round(cor_vals, 2L),
    stringsAsFactors = FALSE
  )
  cor_data$label <- paste0(cor_data$metric, ": tau=", cor_data$tau)
  plot_data <- merge(as.data.frame(ews_object), cor_data, by = "metric")
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("std"),
      color = !!rlang::sym("metric")
    )
  ) +
    ggplot2::geom_line(alpha = 0.8, linewidth = 0.5) +
    ggplot2::facet_wrap(
      ggplot2::vars(!!rlang::sym("label")),
      scales = "free_y"
    ) +
    ggplot2::labs(
      title  = "Multivariate Early Warning Signals (Rolling Window)",
      subtitle = paste0(
        "Each panel shows the trend of an EWS metric.",
        " Kendall's tau indicates trend strength."
      ),
      y = "Standardized EWS",
      x = onlyif(!include_dr, "Time Point")
    ) +
    mews_theme() +
    ggplot2::theme(
      legend.position  = "none",
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = 8)
    )

  if (include_dr) {
    p <- p + ggplot2::theme(
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }
  if (include_dr && !is.null(attr(ews_object, "dimension_reduction"))) {
    p_dr <- mews_plot_dr_(ews_object) +
      ggplot2::theme(legend.position = "right")
    patchwork::wrap_plots(p_dr, p, ncol = 1L, heights = c(2, 3.5))
  } else {
    p
  }
}

mews_plot_dr_ <- function(ews_object) {
  dr_data <- attr(ews_object, "dimension_reduction")
  stopifnot_(
    !is.null(dr_data),
    "No dimension reduction data found in {.arg ews_object}."
  )
  plot_df <- tidyr::pivot_longer(
    dr_data,
    cols      = c(
      !!rlang::sym("maf1_scaled"),
      !!rlang::sym("pc1_scaled")
    ),
    names_to  = "dimension",
    values_to = "scaled_value"
  )
  plot_df$dimension <- ifelse(
    plot_df$dimension == "maf1_scaled", "MAF1", "PC1"
  )
  plot_df$warning <- ifelse(
    plot_df$dimension == "MAF1",
    plot_df$maf1_warning,
    plot_df$pc1_warning
  )
  y_vals <- plot_df$scaled_value[!is.na(plot_df$scaled_value)]
  if (length(y_vals) > 0L) {
    y_range  <- range(y_vals, na.rm = TRUE)
    y_buffer <- diff(y_range) * 0.05
    y_min    <- max(y_range[1L] - y_buffer, -6)
    y_max    <- min(y_range[2L] + y_buffer,  6)
  } else {
    y_min <- -4
    y_max <-  4
  }
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x     = !!rlang::sym("time"),
      y     = !!rlang::sym("scaled_value"),
      color = !!rlang::sym("dimension")
    )
  ) +
    ggplot2::geom_line(alpha = 0.8, linewidth = 1) +
    ggplot2::geom_point(
      data  = function(d) d[d$warning == TRUE, ],
      ggplot2::aes(
        x     = !!rlang::sym("time"),
        y     = !!rlang::sym("scaled_value"),
        color = !!rlang::sym("dimension")
      ),
      size  = 3,
      alpha = 0.8
    ) +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
    ggplot2::labs(
      title    = "Dimension Reduction Time Series",
      subtitle = "MAF1 and PC1 components with warning points highlighted",
      y        = "Scaled Component Value",
      x        = "Time Point",
      color    = "Component"
    ) +
    ggplot2::scale_color_manual(
      values = c("MAF1" = "#3B82F6", "PC1" = "#F59E0B")
    ) +
    mews_theme() +
    ggplot2::theme(
      legend.position    = "left",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    )
}

mews_plot_ribbon_ <- function(classification_data, state_colors) {
  ggplot2::ggplot(
    classification_data,
    ggplot2::aes(x = !!rlang::sym("time"), y = 1)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = !!rlang::sym("state")),
      height = 1
    ) +
    ggplot2::scale_fill_manual(
      values = state_colors,
      name   = "System State",
      drop   = FALSE
    ) +
    ggplot2::labs(x = "Time Point", y = "") +
    mews_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )
}


#' Plot Potential Analysis Results
#'
#' Visualizes the estimated potential landscape, probability density, or
#' both. Wells are marked as points and barriers as vertical dashed lines.
#'
#' @details
#' Three plot types are available:
#'
#' **`type = "landscape"` (default).** Plots the quasi-potential
#' \eqn{U(x) = -\log P(x)} as a function of state-space position
#' \eqn{x}. Wells (local minima) are shown as blue points and barriers
#' (local maxima) as red vertical dashed lines. This view directly shows
#' the stability landscape: valleys are stable attractors, peaks are
#' unstable equilibria, and the depth of each valley relative to the
#' nearest peak quantifies how resilient each attractor is to
#' perturbations.
#'
#' **`type = "density"`.** Plots the kernel density estimate \eqn{P(x)}.
#' Peaks in density correspond to wells in potential and vice versa. This
#' view is useful for checking the quality of the density estimate and
#' confirming that the identified attractors correspond to genuine modes
#' in the data distribution.
#'
#' **`type = "both"`.** Stacks the density plot on top and the potential
#' landscape below using [patchwork::wrap_plots()]. The shared x-axis
#' makes it easy to see how features in density map to features in
#' potential.
#'
#' When the object contains rolling-window results (i.e., `window` was
#' specified), an additional panel showing the number of wells over time
#' is appended below the landscape panel(s).
#'
#' @export
#' @param x \[`potential`\]\cr
#'   An object of class `potential` as returned by [potential_analysis()].
#' @param type \[`character(1)`: `"landscape"`\]\cr
#'   Plot type: `"landscape"`, `"density"`, or `"both"`.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object (or a
#'   [patchwork::wrap_plots()] composite when `type = "both"` or when
#'   rolling results are present).
#'
#' @seealso [potential_analysis()] for computing the potential landscape.
#' @examples
#' \donttest{
#' set.seed(42)
#' bimodal <- c(rnorm(1000, -2, 0.5), rnorm(1000, 2, 0.5))
#' pa <- potential_analysis(bimodal)
#' plot(pa, type = "landscape")
#' plot(pa, type = "density")
#' plot(pa, type = "both")
#' }
plot.potential <- function(x, type = "landscape", ...) {
  check_missing(x)
  check_class(x, "potential")
  type <- check_match(type, c("landscape", "density", "both"))
  colors <- potential_colors_()
  landscape <- x$landscape
  wells <- x$wells
  barriers <- x$barriers
  panels <- list()
  if (type == "density" || type == "both") {
    p_density <- ggplot2::ggplot(
      landscape,
      ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("density"))
    ) +
      ggplot2::geom_line(
        linewidth = 0.7,
        color = colors["density"]
      ) +
      ggplot2::geom_area(
        alpha = 0.15,
        fill = colors["density"]
      )
    if (nrow(barriers) > 0L) {
      p_density <- p_density +
        ggplot2::geom_vline(
          xintercept = barriers$location,
          linetype = "dashed",
          color = colors["barrier"],
          linewidth = 0.5
        )
    }
    p_density <- p_density +
      ggplot2::labs(
        title = "Probability Density Estimate",
        x = ifelse_(type == "both", NULL, "State (x)"),
        y = "Density P(x)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.title = ggplot2::element_text(color = "black", face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )
    if (type == "both") {
      p_density <- p_density +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
    if (type == "density") {
      return(p_density)
    }
    panels <- c(panels, list(p_density))
  }
  if (type == "landscape" || type == "both") {
    p_potential <- ggplot2::ggplot(
      landscape,
      ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("potential"))
    ) +
      ggplot2::geom_line(
        linewidth = 0.7,
        color = colors["potential"]
      )
    # Mark wells
    if (nrow(wells) > 0L) {
      # Find potential values at well locations via interpolation
      well_potential <- stats::approx(
        landscape$x, landscape$potential,
        xout = wells$location, rule = 2
      )$y
      well_df <- data.frame(
        x = wells$location,
        y = well_potential
      )
      p_potential <- p_potential +
        ggplot2::geom_point(
          data = well_df,
          ggplot2::aes(
            x = !!rlang::sym("x"),
            y = !!rlang::sym("y")
          ),
          color = colors["well"],
          size = 3,
          shape = 19,
          inherit.aes = FALSE
        )
    }
    if (nrow(barriers) > 0L) {
      p_potential <- p_potential +
        ggplot2::geom_vline(
          xintercept = barriers$location,
          linetype = "dashed",
          color = colors["barrier"],
          linewidth = 0.5
        )
    }
    p_potential <- p_potential +
      ggplot2::labs(
        title = ifelse_(
          type == "both",
          "Potential Landscape",
          "Stability Landscape: Potential U(x) = -log P(x)"
        ),
        x = "State (x)",
        y = "Potential U(x)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.title = ggplot2::element_text(color = "black", face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )

    if (type == "landscape" && is.null(x$rolling)) {
      return(p_potential)
    }
    panels <- c(panels, list(p_potential))
  }
  if (!is.null(x$rolling)) {
    rolling <- x$rolling
    p_rolling <- ggplot2::ggplot(
      rolling,
      ggplot2::aes(
        x = !!rlang::sym("time"),
        y = !!rlang::sym("n_wells")
      )
    ) +
      ggplot2::geom_step(linewidth = 0.6, color = "#2196F3") +
      ggplot2::scale_y_continuous(
        breaks = function(lims) {
          seq(floor(lims[1L]), ceiling(lims[2L]), by = 1L)
        }
      ) +
      ggplot2::labs(
        title = "Number of Wells Over Time",
        x = "Time",
        y = "Number of Wells"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 13, face = "bold"),
        axis.title = ggplot2::element_text(color = "black", face = "bold"),
        axis.text = ggplot2::element_text(color = "black")
      )
    panels <- c(panels, list(p_rolling))
  }
  if (length(panels) == 1L) {
    return(panels[[1L]])
  }
  patchwork::wrap_plots(panels, ncol = 1L)
}

#' Colour palette for potential landscape plots
#' @noRd
potential_colors_ <- function() {
  c(
    well = "#2196F3",
    barrier = "#F44336",
    potential = "#333333",
    density = "#4CAF50"
  )
}

#' Plot Resilience Results
#'
#' Visualizes resilience metric trajectories alongside the original time
#' series, either as a combined ribbon plot or as faceted line plots.
#'
#' @details
#' Two complementary visualizations are available:
#'
#' **`type = "ribbons"` (default).** Produces a single-panel plot combining
#' the original time series (black line, upper portion) with a set of
#' gradient-colored horizontal bands below it. Each band corresponds to one
#' resilience metric (plus the composite score). Color at each time point
#' encodes the directional score from [classify_resilience()], ranging from
#' green ("Excellent") through yellow ("Fair") to red ("Troubled"). This
#' layout answers the question *"when and where does the system shift between
#' resilience states?"* at a glance, because all metrics share the same time
#' axis. Requires classified data; call [classify_resilience()] first.
#'
#' **`type = "lines"`.** Produces a faceted plot with one panel per raw
#' metric, each showing the metric's rolling-window value over time. No
#' classification is needed. This view answers the question *"what is the
#' raw trajectory of each metric?"* and is useful for diagnosing which
#' metrics drive changes visible in the ribbon plot. Panels use free y-axes
#' so that metrics on different scales are individually legible.
#'
#' @export
#' @param x \[`resilience`\]\cr
#'   An object of class `resilience`, optionally augmented by
#'   [classify_resilience()] (required for `type = "ribbons"`).
#' @param type \[`character(1)`: `"ribbons"`\]\cr
#'   Plot type. The available options are:
#'   * `"ribbons"`: Time series with gradient-colored risk ribbon below.
#'     Requires classified data from [classify_resilience()].
#'   * `"lines"`: Faceted line plots of raw metric values.
#' @param ... Additional arguments (currently unused).
#' @return A [ggplot2::ggplot()] object that can be further customized with
#'   standard ggplot2 functions (e.g., `+ theme()`, `+ labs()`).
#'
#' @seealso [resilience()] for computing the metrics;
#'   [classify_resilience()] for directional normalization required by the
#'   ribbon plot.
#' @examples
#' \donttest{
#' set.seed(20)
#' x <- cumsum(rnorm(500))
#' res <- resilience(x, window = 50L)
#'
#' # Faceted raw metric trajectories (no classification needed)
#' plot(res, type = "lines")
#'
#' # Ribbon plot (requires classification first)
#' cls <- classify_resilience(res)
#' plot(cls, type = "ribbons")
#' }
plot.resilience <- function(x, type = "ribbons", ...) {
  check_missing(x)
  check_class(x, "resilience")
  type <- check_match(type, c("ribbons", "lines"))
  ifelse_(
    type == "ribbons",
    plot_resilience_ribbons_(x, ...),
    plot_resilience_lines_(x, ...)
  )
}

resilience_display_names <- function() {
  c(
    vsi            = "Stability",
    arch_lm        = "Volatility",
    cv             = "Variability",
    recovery_time  = "Recovery time",
    recovery_slope = "Recovery rate",
    sample_entropy = "Complexity",
    dfa_alpha      = "Memory",
    ac_ratio       = "Correlation",
    composite      = "OVERALL"
  )
}

plot_resilience_ribbons_ <- function(x, show_metrics = NULL, ...) {
  has_scores <- any(grepl("_score$", names(x)))
  stopifnot_(
    has_scores,
    c(
      "Resilience ribbon plot requires classified data.",
      `i` = "Run {.fn classify_resilience} first."
    )
  )
  display <- resilience_display_names()
  all_metrics <- attr(x, "metrics")
  # Default: exclude recovery_time and recovery_slope (matches original)
  if (is.null(show_metrics)) {
    plot_metrics <- all_metrics[
      !all_metrics %in% c("recovery_time", "recovery_slope")
    ]
  } else if (length(show_metrics) == 1L && show_metrics == "all") {
    plot_metrics <- all_metrics
  } else {
    plot_metrics <- show_metrics
  }
  # Order: composite first, then individual metrics
  ordered_metrics <- c("composite", plot_metrics)
  time_vals <- x[["time"]]
  ts_vals   <- x[["value"]]
  ts_range <- diff(range(ts_vals, na.rm = TRUE))
  ts_min   <- min(ts_vals, na.rm = TRUE)
  ts_max   <- max(ts_vals, na.rm = TRUE)
  # Compress time series into upper 60% of its range
  ts_compressed_range <- ts_range * 0.6
  ts_new_max <- ts_max
  ts_new_min <- ts_new_max - ts_compressed_range
  # Ribbon sizing
  ribbon_height  <- ts_range * 0.08
  ribbon_spacing <- ribbon_height * 1.3
  ribbon_start_y <- ts_new_min - ts_range * 0.15
  ts_df <- data.frame(time = time_vals, value = ts_vals)
  p <- ggplot2::ggplot(ts_df, ggplot2::aes(
    x = !!rlang::sym("time"), y = !!rlang::sym("value")
  )) +
    ggplot2::geom_line(color = "black", linewidth = 0.5) +
    ggplot2::labs(
      title = "System Resilience Analysis",
      x = "Time",
      y = "Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.text = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(color = "black", face = "bold"),
      axis.text = ggplot2::element_text(color = "black")
    )
  time_diff <- ifelse_(
    length(time_vals) > 1L,
    stats::median(diff(time_vals), na.rm = TRUE),
    1
  )
  t_min <- min(time_vals, na.rm = TRUE)
  t_max <- max(time_vals, na.rm = TRUE)
  for (i in seq_along(ordered_metrics)) {
    metric <- ordered_metrics[i]
    score_col <- ifelse_(
      metric == "composite",
      "composite_score",
      paste0(metric, "_score")
    )
    if (!score_col %in% names(x)) {
      next
    }
    ribbon_y_bottom <- ribbon_start_y - (i * ribbon_spacing)
    ribbon_y_center <- ribbon_y_bottom + ribbon_height / 2
    scores <- x[[score_col]]
    scores[is.na(scores)] <- -1
    tile_df <- data.frame(
      x = time_vals,
      y = ribbon_y_center,
      severity = scores
    )
    p <- p + ggplot2::geom_tile(
      data   = tile_df,
      ggplot2::aes(
        x = !!rlang::sym("x"),
        y = !!rlang::sym("y"),
        fill = !!rlang::sym("severity")
      ),
      width  = time_diff,
      height = ribbon_height,
      alpha  = 0.8
    )
    label_color <- ifelse_(metric == "composite", "darkblue", "black")
    label_size  <- ifelse_(metric == "composite", 3.5, 3)
    disp_name  <- ifelse_(metric == "composite", "OVERALL", display[metric])
    p <- p + ggplot2::annotate(
      "text",
      x = t_max + (t_max - t_min) / 70,
      y = ribbon_y_center,
      label = disp_name,
      hjust = 0,
      size = label_size,
      fontface = "bold",
      color = label_color
    )
  }
  # Extend x-axis to make room for labels; set y limits
  n_metrics <- length(ordered_metrics)
  y_bottom  <- ribbon_start_y - (n_metrics * ribbon_spacing) - ribbon_height
  p <- p +
    ggplot2::coord_cartesian(
      xlim = c(t_min, t_max + (t_max - t_min) / 6),
      ylim = c(y_bottom, ts_new_max + ts_range * 0.1),
      clip = "off"
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c(
        "#cccccc",
        "#2d7d32",
        "#66bb6a",
        "#ffeb99",
        "#ffcc66",
        "#ff9999",
        "#cc0000"
      ),
      values = (c(-1, 0, 0.2, 0.4, 0.6, 0.8, 1.0) + 1) / 2,
      name = "Resilience Level",
      breaks = c(-1, 0, 0.2, 0.4, 0.6, 0.8, 1.0),
      labels = c(
        "Insufficient", "Excellent", "Solid", "Fair",
        "Vulnerable", "Failing", "Troubled"
      ),
      na.value = "#cccccc",
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barwidth       = 24,
        barheight      = 1
      )
    )
  p
}

plot_resilience_lines_ <- function(x, ...) {
  metrics <- attr(x, "metrics")
  display <- resilience_display_names()
  long <- tidyr::pivot_longer(
    dplyr::select(x, !!rlang::sym("time"), tidyselect::all_of(metrics)),
    cols = tidyselect::all_of(metrics),
    names_to = "metric",
    values_to = "score"
  )
  long <- dplyr::mutate(
    long,
    metric_label = factor(
      display[!!rlang::sym("metric")],
      levels = display[metrics]
    )
  )
  ggplot2::ggplot(
    long,
    ggplot2::aes(
      x = !!rlang::sym("time"),
      y = !!rlang::sym("score"),
      color = !!rlang::sym("metric_label")
    )
  ) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::facet_wrap(
      ggplot2::vars(!!rlang::sym("metric_label")),
      scales = "free_y",
      ncol = 2L
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Time", y = "Value", color = "Metric") +
    ggplot2::theme(legend.position = "none")
}
