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

#' Summarize EWS Sensitivity Analysis Results
#'
#' Summary method for `"sensitivity_ews"` objects. Computes the range of
#' Kendall tau values, identifies the most and least robust parameter
#' combinations, and assesses overall consistency.
#'
#' @export
#' @param object \[`sensitivity_ews`\]\cr
#'   A `sensitivity_ews` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `metric`, `tau_range`, `tau_mean`,
#'   `n_positive`, `n_negative`, `n_total`, `most_robust` (combination
#'   with highest |tau|), `least_robust` (combination with lowest |tau|),
#'   and `consistent` (logical), returned invisibly.
summary.sensitivity_ews <- function(object, ...) {
  check_missing(object)
  check_class(object, "sensitivity_ews")
  metric_name <- attr(object, "metric") %||% "unknown"
  tau_df <- unique(object[, c("window", "detrend", "tau")])
  tau_vals <- tau_df$tau[!is.na(tau_df$tau)]
  n_total <- length(tau_vals)
  n_positive <- sum(tau_vals > 0)
  n_negative <- sum(tau_vals < 0)
  tau_mean <- mean(tau_vals, na.rm = TRUE)
  tau_range <- range(tau_vals, na.rm = TRUE)
  # Most robust: highest absolute tau
  abs_tau <- abs(tau_df$tau)
  most_idx <- which.max(abs_tau)
  least_idx <- which.min(abs_tau)
  most_robust <- tau_df[most_idx, ]
  least_robust <- tau_df[least_idx, ]
  # Consistency: all tau values have the same sign
  consistent <- n_total > 0L && (n_positive == n_total || n_negative == n_total)
  structure(
    list(
      metric = metric_name,
      tau_range = tau_range,
      tau_mean = tau_mean,
      n_positive = n_positive,
      n_negative = n_negative,
      n_total = n_total,
      most_robust = most_robust,
      least_robust = least_robust,
      consistent = consistent
    ),
    class = "summary.sensitivity_ews"
  )
}

#' Summarize Spectral EWS Analysis Results
#'
#' Summary method for `"spectral"` objects. Shows the state distribution,
#' mean spectral exponent, and analysis settings.
#'
#' @export
#' @param object \[`spectral`\]\cr
#'   A `spectral` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `mean_beta`, `mean_ratio`, `mean_r_squared`,
#'   `state_counts` (if states present), `window`, `method`, `detrend`, and
#'   `align`, returned invisibly.
summary.spectral <- function(object, ...) {
  check_missing(object)
  check_class(object, "spectral")
  mean_beta <- mean(object$spectral_exponent, na.rm = TRUE)
  mean_ratio <- mean(object$spectral_ratio, na.rm = TRUE)
  mean_r2 <- mean(object$r_squared, na.rm = TRUE)
  state_counts <- NULL
  if ("state" %in% names(object)) {
    state_counts <- table(object$state)
  }
  structure(
    list(
      n = nrow(object),
      mean_beta = mean_beta,
      mean_ratio = mean_ratio,
      mean_r_squared = mean_r2,
      state_counts = state_counts,
      window = attr(object, "window"),
      method = attr(object, "method"),
      detrend = attr(object, "detrend"),
      align = attr(object, "align")
    ),
    class = "summary.spectral"
  )
}

#' Summarize Surrogate Test Results
#'
#' Summary method for `"surrogate_test"` objects. Provides detailed statistics
#' of the surrogate distribution including quantiles and the position of the
#' observed value.
#'
#' @export
#' @param object \[`surrogate_test`\]\cr
#'   A `surrogate_test` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `observed_tau`, `p_value`, `significant`,
#'   `surrogate_quantiles` (named numeric vector), `surrogate_mean`,
#'   `surrogate_sd`, `method`, `metric`, `window`, and `n_surrogates`.
summary.surrogate_test <- function(object, ...) {
  check_missing(object)
  check_class(object, "surrogate_test")
  valid_taus <- object$surrogate_taus[!is.na(object$surrogate_taus)]
  surr_quantiles <- stats::quantile(
    valid_taus, probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1),
    na.rm = TRUE
  )
  surr_mean <- mean(valid_taus, na.rm = TRUE)
  surr_sd <- stats::sd(valid_taus, na.rm = TRUE)
  structure(
    list(
      observed_tau = object$observed_tau,
      valid_taus = valid_taus,
      p_value = object$p_value,
      significant = object$significant,
      surrogate_quantiles = surr_quantiles,
      surrogate_mean = surr_mean,
      surrogate_sd = surr_sd,
      method = object$method,
      metric = object$metric,
      window = object$window,
      n_surrogates = object$n_surrogates
    ),
    class = "summary.surrogate_test"
  )
}

#' Summarize Trend Classification Results
#'
#' Summary method for `"trend"` objects. Shows the distribution of trend
#' classifications and the analysis settings.
#'
#' @export
#' @param object \[`trend`\]\cr
#'   A `trend` object.
#' @param ... Additional arguments (currently unused).
#' @return A list with elements `counts` (frequency table of states),
#'   `proportions` (relative frequencies), `window`, `method`, and `align`.
summary.trend <- function(object, ...) {
  check_missing(object)
  check_class(object, "trend")
  counts <- table(object$state)
  props  <- prop.table(counts)
  structure(
    list(
      n = nrow(object),
      counts = counts,
      proportions = props,
      window = attr(object, "window"),
      method = attr(object, "method"),
      align = attr(object, "align")
    ),
    class = "summary.trend"
  )
}
