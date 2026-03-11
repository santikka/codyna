#' Summarize Changepoint Detection Results
#'
#' Prints a structured summary of the changepoint analysis including the
#' number and location of detected changepoints, the direction and
#' magnitude of each change, and per-segment statistics.
#'
#' @describeIn detect_cpts Summary method for `"changepoint"`
#'   objects.
#'
#' @export
#' @param object \[`changepoint`]\cr
#'   A `changepoint` object.
#' @param ... Additional arguments (currently unused).
#' @return A `summary.changepoint` object, which is a list with the following
#'   elements:
#'
#'   * `cpts`: The original `changepoint` object.
#'   * `changes` (a `data.frame` with columns `location`, `from_mean`,
#'     `to_mean`, `mean_shift`, `from_var`, `to_var`)
#'   * `segments` (a `data.frame` with columns `segment`, `start`, `end`,
#'     `length`, `mean`, `variance`).
#'
summary.changepoint <- function(object, ...) {
  check_missing(object)
  check_class(object, "changepoint")
  n_cpt <- attr(object, "n_changepoints")
  cpt_locs <- attr(object, "changepoint_locations")
  seg_ids <- sort(unique(object$segment))
  seg_stats <- do.call(
    base::rbind,
    lapply(
      seg_ids,
      function(s) {
        rows <- object[object$segment == s, ]
        data.frame(
          segment = s,
          start = min(rows$time),
          end = max(rows$time),
          length = nrow(rows),
          mean = round(rows$segment_mean[1L], 4),
          variance = round(rows$segment_var[1L], 4),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  changes <- NULL
  if (n_cpt > 0L) {
    changes <- data.frame(
      location = integer(n_cpt),
      from_mean = numeric(n_cpt),
      to_mean = numeric(n_cpt),
      mean_shift = numeric(n_cpt),
      from_var = numeric(n_cpt),
      to_var = numeric(n_cpt),
      stringsAsFactors = FALSE
    )
    for (i in seq_len(n_cpt)) {
      from_seg <- seg_stats[seg_stats$segment == i, ]
      to_seg <- seg_stats[seg_stats$segment == i + 1L, ]
      shift <- to_seg$mean - from_seg$mean
      direction <- ifelse_(shift > 0, "increase", "decrease")
      cpt_type_val <- object$changepoint_type[object$changepoint][i]
      changes$location[i] <- cpt_locs[i]
      changes$from_mean[i] <- from_seg$mean
      changes$to_mean[i] <- to_seg$mean
      changes$mean_shift[i] <- round(shift, 4)
      changes$from_var[i] <- from_seg$variance
      changes$to_var[i] <- to_seg$variance
    }
  }
  structure(
    list(
      cpts = object,
      changes = changes,
      segments = seg_stats
    ),
    class = "summary.changepoint"
  )
}

#' Summarize Hurst Early Warning Signal Results
#'
#' Prints a compact summary of the warning level distribution across all
#' time points, including counts and percentages for each level and the
#' maximum warning level observed.
#'
#' @export
#' @param object \[`hurst_ews`]\cr
#'   A `hurst_ews` object as returned by [detect_hurst_warnings()].
#' @param ... Additional arguments (currently unused).
#' @return `summary.hurst_ews` object. A `list` containing the warning
#'   labels and the warning levels.
summary.hurst_ews <- function(object, ...) {
  check_missing(object)
  check_class(object, "hurst_ews")
  labels <- c("none", "low", "moderate", "high", "critical")
  counts <- table(factor(object$warning_label, levels = labels))
  max_level <- max(object$warning_level, na.rm = TRUE)
  structure(
    list(
      n = nrow(object),
      labels = labels,
      counts = counts,
      max_level = max(object$warning_level, na.rm = TRUE)
    ),
    class = c("summary.hurst_ews", "list")
  )
}

#' Summarize a Multivariate EWS Object
#'
#' @describeIn detect_multivariate_warnings Summary method for `"multi_ews"`
#'   objects. Provides a condensed overview including per-metric warning
#'   counts and system state distribution.
#'
#' @export
#' @param object \[`multi_ews`]\cr
#'   A `multi_ews` object.
#' @param ... Additional arguments (currently unused).
#' @return A `summary.multi_ews` object. Identical to `object`.
summary.multi_ews <- function(object, ...) {
  structure(
    object,
    class = c("summary.multi_ews", "data.frame")
  )
}
