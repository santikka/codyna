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
