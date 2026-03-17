#' Print Changepoint Detection Results
#'
#' @export
#' @param x \[`changepoint`]\cr A `changepoint` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.changepoint <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a Summary of Changepoint Detection Results
#'
#' @export
#' @param x \[`summary.changepoint`]\cr A `summary.changepoint` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
print.summary.changepoint <- function(x, ...) {
  check_missing(x)
  check_class(x, "summary.changepoint")
  cpts <- x$cpts
  changes <- x$changes
  seg_stats <- x$segments
  n_cpt <- attr(cpts, "n_changepoints")
  cpt_locs <- attr(cpts, "changepoint_locations")
  cat("Changepoint Detection Summary\n\n")
  cat("Method           :", attr(cpts, "method"), "\n")
  cat("Change type      :", attr(cpts, "type"), "\n")
  cat("Penalty          :", attr(cpts, "penalty"), "\n")
  cat("Min segment      :", attr(cpts, "min_segment"), "\n")
  cat("N observations   :", nrow(cpts), "\n")
  cat("N changepoints   :", n_cpt, "\n")
  if (n_cpt > 0L) {
    cat("Locations        :", paste(cpt_locs, collapse = ", "), "\n")
  }
  cat("\nSegment statistics:\n")
  for (i in seq_len(nrow(seg_stats))) {
    seg_regime <- cpts$regime[cpts$segment == seg_stats$segment[i]][1L]
    seg_state  <- as.character(
      cpts$state[cpts$segment == seg_stats$segment[i]][1L]
    )
    cat(
      sprintf(
        paste0(
          "  Segment %d (regime %d) [%s]: t = [%s, %s], ",
          "n = %d, mean = %.4f, var = %.4f\n"
        ),
        seg_stats$segment[i],
        seg_regime,
        seg_state,
        seg_stats$start[i],
        seg_stats$end[i],
        seg_stats$length[i],
        seg_stats$mean[i],
        seg_stats$variance[i]
      )
    )
  }
  if (n_cpt > 0L) {
    cat("\nChange details:\n")
    for (i in seq_len(n_cpt)) {
      from_seg <- seg_stats[seg_stats$segment == i, ]
      to_seg <- seg_stats[seg_stats$segment == i + 1L, ]
      shift <- to_seg$mean - from_seg$mean
      direction <- ifelse_(shift > 0, "increase", "decrease")
      cpt_type_val <- cpts$changepoint_type[cpts$changepoint][i]
      cat(
        sprintf(
          paste0(
            "  CP %d at t = %d [%s]: mean %.4f -> %.4f (%s of %.4f) ",
            "var %.4f -> %.4f\n"
          ),
          i, cpt_locs[i], cpt_type_val %||% "change",
          from_seg$mean, to_seg$mean, direction, abs(shift),
          from_seg$variance, to_seg$variance
        )
      )
    }
  }
  invisible(x)
}

#' Print EWS Detection Results
#'
#' @export
#' @param x \[`ews`]\cr
#'   EWS detection result from [detect_warnings()].
#' @param ... Additional arguments passed to the generic print method.
#' @return `x` (invisibly).
#' @examples
#' set.seed(123)
#' ts_data <- stats::arima.sim(list(order = c(1, 1, 0), ar = 0.6), n = 200)
#' ews <- detect_warnings(ts_data)
#' print(ews)
#'
print.ews <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print Regime Detection Results
#'
#' @export
#' @param x \[`regimes`]\cr
#'   Regime detection result from [detect_regimes()].
#' @param ... Additional arguments passed to the generic print method.
#' @return `x` (invisibly).
#' @examples
#' set.seed(123)
#' ts_data <- stats::arima.sim(list(order = c(1, 1, 0), ar = 0.6), n = 200)
#' regimes <- detect_regimes(
#'   data = ts_data,
#'   method = "threshold",
#'   sensitivity = "medium"
#' )
#' print(regimes)
#'
print.regimes <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print Discovered Patterns
#'
#' @export
#' @param x \[`patterns`]\cr
#'   Pattern discovery result from [discover_patterns()].
#' @param ... Additional arguments passed to the generic print method.
#' @return `x` (invisibly).
#' @examples
#' ngrams <- discover_patterns(engagement, type = "ngram")
#' print(ngrams)
#'
print.patterns <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a Hurst Rolling Analysis Object
#'
#' Prints the underlying tibble of rolling Hurst exponent results. Delegates
#' to the default tibble print method.
#'
#' @export
#' @param x \[`hurst`]\cr
#'   A `hurst` object as returned by [hurst()] with `states = TRUE`.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.hurst <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print Hurst Early Warning Signal Results
#'
#' Prints the underlying tibble of early warning signal results. Delegates
#' to the default tibble print method.
#'
#' @export
#' @param x \[`hurst_ews`]\cr
#'   A `hurst_ews` object as returned by [detect_hurst_warnings()].
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.hurst_ews <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a Global Hurst Exponent Result
#'
#' Prints a concise summary of a single global Hurst exponent estimate,
#' including the method, series length, exponent value, \eqn{R^2}, and
#' state classification.
#'
#' @export
#' @param x \[`hurst_global`]\cr
#'   A `hurst_global` object as returned by [hurst()] with
#'   `states = FALSE`.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
print.hurst_global <- function(x, ...) {
  check_missing(x)
  check_class(x, "hurst_global")
  cat("Hurst Exponent Analysis\n\n")
  cat("Method :", x$method, "\n")
  cat("N      :", x$n, "\n")
  cat("Hurst  :", round(x$hurst, 4L), "\n")
  cat("R\u00b2     :", round(x$r_squared, 4L), "\n")
  cat("State  :", hurst_classify(x$hurst), "\n")
  invisible(x)
}

#' Print Hurst Early Warning Signal Analysis Summary
#'
#' @export
#' @param x \[`summary.hurst_ews`]\cr A `summary.hurst_ews` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
print.summary.hurst_ews <- function(x, ...) {
  check_missing(x)
  check_class(x, "summary.hurst_ews")
  n <- x$n
  counts <- x$counts
  labels <- x$labels
  max_level <- x$max_level
  cat("Hurst Early Warning Signal Summary\n\n")
  cat("Number of time points:", n, "\n")
  cat("Max warning level:", labels[max_level + 1L], "\n")
  cat("Warning distribution:\n")
  for (i in seq_along(labels)) {
    pct <- round(100 * counts[i] / n, 1L)
    cat("  ", labels[i], ": ", counts[i], " (", pct, "%)\n", sep = "")
  }
  invisible(x)
}

#' Print a Multivariate EWS Object
#'
#' Print method for `"multi_ews"` objects.
#' Delegates to tibble's default print.
#'
#' @export
#' @param x \[`multi_ews`]\cr A `multi_ews` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.multi_ews <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a Summary of a  Multivariate EWS Object
#'
#' @export
#' @param x \[`summary.multi_ews`]\cr A `summary.multi_ews` object.
#' @param ... Not used.
#' @return `x`, invisibly.
print.summary.multi_ews <- function(x, ...) {
  check_missing(x)
  check_class(x, "summary.multi_ews")
  method <- attr(x, "method")
  cat("Multivariate EWS Summary\n\n")
  cat("Method:", method, "\n")
  cat("Total time points:", length(unique(x$time)), "\n")
  cat("Metrics computed:", paste(unique(x$metric), collapse = ", "), "\n")
  cat("\n")
  if (method == "expanding") {
    cat("Per-metric warnings:\n")
    warned <- x[x$detected == 1L, ]
    if (nrow(warned) > 0L) {
      tbl <- table(warned$metric)
      for (nm in names(tbl)) {
        cat("  ", nm, ":", tbl[[nm]], "time points\n")
      }
    } else {
      cat("  No warnings detected.\n")
    }
    cls <- attr(x, "classification")
    if (!is.null(cls)) {
      cat("\nSystem state distribution:\n")
      tbl <- table(cls$state)
      total <- sum(tbl)
      for (nm in names(tbl)) {
        pct <- round(100 * tbl[[nm]] / total, 1L)
        cat("  ", nm, ":", tbl[[nm]], "(", pct, "%)\n")
      }
    }
  }
  if (method == "rolling") {
    cor_vals <- attr(x, "cor")
    if (!is.null(cor_vals)) {
      cat("Kendall's tau trend statistics:\n")
      tau_clean <- cor_vals
      names(tau_clean) <- sub("\\.tau$", "", names(tau_clean))
      strong <- names(tau_clean)[tau_clean > 0.7 & !is.na(tau_clean)]
      if (length(strong) > 0L) {
        cat(
          "  Strong upward trends (tau > 0.7):",
          paste(strong, collapse = ", "),
          "\n"
        )
      } else {
        cat("  No metrics show strong upward trends.\n")
      }
    }
  }
}


#' Print a Potential Analysis Object
#'
#' Prints a concise summary of the potential analysis: the number of
#' wells and barriers, well locations, and analysis settings.
#'
#' @export
#' @param x \[`potential`\]\cr
#'   A `potential` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
print.potential <- function(x, ...) {
  cat("Potential Analysis\n\n")
  cat("Detrend  :", attr(x, "detrend"), "\n")
  cat("N points :", length(x$values), "\n")
  cat("N bins   :", attr(x, "n_bins"), "\n")
  bw_str <- ifelse_(
    is.null(attr(x, "bandwidth")),
    "auto (Silverman)",
    as.character(round(attr(x, "bandwidth"), 4))
  )
  cat("Bandwidth:", bw_str, "\n")
  win_str <- ifelse_(
    is.null(attr(x, "window")),
    "global",
    as.character(attr(x, "window"))
  )
  cat("Window   :", win_str, "\n")
  cat("Wells    :", x$n_wells, "\n")
  cat("Barriers :", nrow(x$barriers), "\n")
  if (x$n_wells > 0L) {
    cat("\nWell locations:\n")
    for (i in seq_len(nrow(x$wells))) {
      cat(
        sprintf(
          "Well %d: x = %.4f  (depth = %.4f, width = %.4f)\n",
          i, x$wells$location[i], x$wells$depth[i], x$wells$width[i]
        )
      )
    }
  }
  if (!is.null(x$rolling)) {
    n_wells_range <- range(x$rolling$n_wells, na.rm = TRUE)
    cat(sprintf(
      "\n  Rolling: n_wells range = [%d, %d]\n",
      n_wells_range[1L], n_wells_range[2L]
    ))
  }
  invisible(x)
}

#' Print a Resilience Object
#'
#' Print method for `"resilience"` objects. Dispatches
#' to the tibble print method, displaying the time series data and computed
#' resilience metrics.
#'
#' @export
#' @param x \[`resilience`\]\cr
#'   A `resilience` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.resilience <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a Sensitivity Object
#'
#' Print method for `"sensitivity_ews"` objects.
#'   Dispatches to the tibble print method.
#'
#' @export
#' @param x \[`sensitivity_ews`\]\cr
#'   A `sensitivity_ews` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.sensitivity_ews <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a EWS Sensitivity Analysis Summary
#'
#' @export
#' @param x \[`summary.sensitivity_ews`]\cr A `summary.sensitivity_ews` object.
#' @param ... Not used.
#' @return `x`, invisibly.
print.summary.sensitivity_ews <- function(x, ...) {
  cat("Sensitivity Analysis Summary\n\n")
  cat("Metric        :", x$metric, "\n")
  cat("Combinations  :", x$n_total, "\n")
  cat(
    "Tau range     : [", sprintf("%.4f", x$tau_range[1L]),
    ", ", sprintf("%.4f", x$tau_range[2L]), "]\n", sep = ""
  )
  cat("Tau mean      :", sprintf("%.4f", x$tau_mean), "\n")
  cat("Positive tau  : ", x$n_positive, "/", x$n_total, "\n", sep = "")
  cat("Negative tau  : ", x$n_negative, "/", x$n_total, "\n\n", sep = "")
  cat(
    "Most robust   : window =", x$most_robust$window,
    ", detrend =", x$most_robust$detrend,
    ", tau =", sprintf("%.4f", x$most_robust$tau), "\n")
  cat(
    "Least robust  : window =", x$least_robust$window,
    ", detrend =", x$least_robust$detrend,
    ", tau =", sprintf("%.4f", x$least_robust$tau), "\n\n")
  if (x$consistent) {
    direction <- ifelse_(x$tau_mean > 0, "positive", "negative")
    cat("Consistency   : ALL tau values are", direction, "\n")
    cat("Assessment    : Signal is ROBUST across parameter choices.\n")
  } else {
    cat("Consistency   : Mixed signs detected\n")
    cat(
      "Assessment    : Signal DEPENDS on parameter choices;",
      "interpret with caution.\n"
    )
  }
}

#' Print a Spectral EWS Analysis Results
#'
#' Print method for `"spectral"` objects. Dispatches to the tibble print method.
#'
#' @export
#' @param x \[`spectral`\]\cr
#'   A `spectral` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.spectral <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a Summary of Spectral EWS Analysis Results
#'
#' @export
#' @param x \[`summary.spectral`]\cr A `summary.spectral` object.
#' @param ... Not used.
#' @return `x`, invisibly.
print.summary.spectral <- function(x, ...) {
  cat("Spectral Early Warning Signal Summary\n\n")
  cat("Method  :", x$method, "\n")
  cat("Detrend :", x$detrend, "\n")
  cat("Window  :", x$window, "\n")
  cat("Align   :", x$align, "\n")
  cat("N       :", x$n, "\n\n")
  cat("Mean spectral exponent (beta):", round(x$mean_beta, 4L), "\n")
  cat("Mean spectral ratio          :", round(x$mean_ratio, 4L), "\n")
  cat("Mean R-squared               :", round(x$mean_r_squared, 4L), "\n")
  if (!is.null(x$state_counts)) {
    props <- prop.table(x$state_counts)
    cat("\nState distribution:\n\n")
    for (i in seq_along(x$state_counts)) {
      if (x$state_counts[i] > 0L) {
        cat(
          sprintf(
            "    %-12s %5d  (%5.1f%%)\n",
            names(x$state_counts)[i],
            x$state_counts[i],
            props[i] * 100
          )
        )
      }
    }
  }
  invisible(x)
}

#' Print a Surrogate Test Object
#'
#' Print method for `"surrogate_test"` objects. Displays the observed tau,
#' p-value, significance, and method settings.
#'
#' @export
#' @param x \[`surrogate_test`\]\cr
#'   A `surrogate_test` object.
#' @param ... Additional arguments (currently unused).
#' @return `x`, invisibly.
print.surrogate_test <- function(x, ...) {
  check_missing(x)
  check_class(x, "surrogate_test")
  cat("Surrogate Significance Test\n\n")
  cat("Method      :", x$method, "\n")
  cat("Metric      :", x$metric, "\n")
  cat("Window      :", x$window, "\n")
  cat("Surrogates  :", x$n_surrogates, "\n")
  cat("Observed tau:", round(x$observed_tau, 4L), "\n")
  cat("p-value     :", surrogate_format_p_(x$p_value), "\n")
  sig_label <- ifelse_(
    is.na(x$significant), "NA",
    ifelse_(x$significant, "YES (p < 0.05)", "NO (p >= 0.05)")
  )
  cat("Significant :", sig_label, "\n")
  invisible(x)
}

#' Print a Summary of Surrogate Test Results
#'
#' @export
#' @param x \[`summary.surrotate_test`]\cr A `summary.surrogate_test` object.
#' @param ... Not used.
#' @return `x`, invisibly.
print.summary.surrogate_test <- function(x, ...) {
  cat("Surrogate Test Summary\n\n")
  cat("Method         :", x$method, "\n")
  cat("Metric         :", x$metric, "\n")
  cat("Window         :", x$window, "\n")
  cat("Test           :", x$test, "\n")
  cat("N surrogates   :", x$n_surrogates, "\n")
  cat("Valid surrogates:", length(x$valid_taus), "\n\n")
  cat("Observed:\n")
  cat("  Kendall's tau  :", round(x$observed_tau, 4L), "\n")
  cat("  p-value        :", surrogate_format_p_(x$p_value), "\n")
  sig_label <- ifelse_(
    is.na(x$significant), "NA",
    ifelse_(x$significant, "YES (p < 0.05)", "NO (p >= 0.05)")
  )
  cat("Significant    :", sig_label, "\n\n")
  cat("Surrogate distribution:\n\n")
  cat("Mean           :", round(x$surr_mean, 4L), "\n")
  cat("SD             :", round(x$surr_sd, 4L), "\n")
  cat("Quantiles:\n")
  q_names <- c("Min", "2.5%", "25%", "50%", "75%", "97.5%", "Max")
  for (i in seq_along(surr_quantiles)) {
    cat(sprintf("    %-6s : %8.4f\n", q_names[i], x$surr_quantiles[i]))
  }
  # Rank of observed tau within surrogates
  if (length(x$valid_taus) > 0L && !is.na(x$observed_tau)) {
    rank_pct <- mean(x$valid_taus <= x$observed_tau) * 100
    cat(
      sprintf(
        "\n  Observed tau ranks at the %.1f%% percentile of surrogates.\n",
        rank_pct
      )
    )
  }
}

#' Print a Trend Object
#'
#' Print method for `"trend"` objects. Dispatches to the tibble print method.
#'
#' @export
#' @param x \[`trend`\]\cr
#'   A `trend` object.
#' @param ... Additional arguments passed to the tibble print method.
#' @return `x`, invisibly.
print.trend <- function(x, ...) {
  NextMethod(generic = "print", object = x, ...)
}

#' Print a Trend Classification Summary
#'
#' @export
#' @param x \[`summary.trend`]\cr A `summary.trend` object.
#' @param ... Not used.
#' @return `x`, invisibly.
print.summary.trend <- function(x, ...) {
  cat("Trend Classification Summary\n\n")
  cat("Method :", x$method, "\n")
  cat("Window :", x$window, "\n")
  cat("Align  :", x$align, "\n")
  cat("N      :", x$n, "\n\n")
  cat("State distribution:\n")
  for (i in seq_along(x$counts)) {
    if (x$counts[i] > 0L) {
      cat(
        sprintf(
          "  %-12s %5d  (%5.1f%%)\n",
          names(x$counts)[i],
          x$counts[i],
          x$props[i] * 100
        )
      )
    }
  }
}
