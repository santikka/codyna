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
#' @family hurst
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
