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
#' @param x \[`changepoint`\]\cr
#'   An object of class `changepoint` as returned by
#'   [detect_cpts()].
#' @param type \[`character(1)`: `"series"`\]\cr
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
  # Map state names to colours by name; fall back to cycling for unknowns
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
