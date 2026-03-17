#' Compute and Classify Rolling Time Series Trends
#'
#' Performs rolling-window trend analysis on a time series, classifying each
#' time point as **ascending**, **descending**, **flat**, or **turbulent**.
#' Turbulence is detected from the volatility of the trend metric itself,
#' not from the raw series variance.
#'
#' @details
#' The analysis proceeds in three stages:
#'
#' **1. Rolling metric calculation.** A trend metric is computed for each
#' window position. The available methods are:
#'
#' * `"slope"` (default): Rate of change estimated by one of four
#'     sub-methods:
#'
#'       * `"ols"`: Ordinary Least Squares slope. Efficient but
#'         sensitive to outliers.
#'       * `"robust"`: Theil-Sen estimator (median of all pairwise
#'         slopes). Highly robust to outliers.
#'       * `"spearman"` / `"kendall"`: Rank-based correlation scaled
#'         by the ratio of standard deviations to produce a slope-like
#'         metric. Robust to non-linear relationships.
#'
#'   * `"ar1_phi1"`: First-order autoregressive coefficient. Values
#'     near 1 indicate a strong trend; near 0 indicate no trend.
#'   * `"growth_factor"`: Ratio of the last to the first value in each
#'     window (\eqn{V_{end} / V_{start}}). Suitable for multiplicative or
#'     exponential processes.
#'
#' **2. Initial classification.** Each point is classified based on the
#' metric value relative to a neutral value (0 for slope, 1 for ar1_phi1
#' and growth_factor) plus or minus `epsilon`:
#'
#' - **ascending**: metric > neutral + epsilon
#' - **descending**: metric < neutral - epsilon
#' - **flat**: metric within the epsilon band
#'
#'
#' **3. Turbulence re-classification.** A combined volatility score is
#' computed on a rolling sub-window of the metric itself:
#' \deqn{V_{combined} = CV(metric) + 0.5 \times \frac{Range(metric)}{|mean(metric)|}}
#' If this score exceeds `turbulence_threshold`, the point is re-classified
#' as **turbulent**. The `flat_to_turbulent_factor` introduces hysteresis:
#' flat segments require a higher volatility score to be reclassified,
#' preventing minor noise from triggering false turbulence detections.
#'
#' @export
#' @param data \[`numeric` or `ts`\]\cr
#'   The time series to analyze. Accepts a numeric vector or a `ts` object.
#' @param window \[`integer(1)`: `NULL`\]\cr
#'   Rolling window size. If `NULL`, an adaptive size is used:
#'   `max(3, min(n, round(n / 10)))`.
#' @param method \[`character(1)`: `"slope"`\]\cr
#'   Trend quantification method. One of `"slope"`, `"ar1_phi1"`, or
#'   `"growth_factor"`.
#' @param slope_method \[`character(1)`: `"robust"`\]\cr
#'   Slope estimation method when `method = "slope"`. One of `"ols"`,
#'   `"robust"`, `"spearman"`, or `"kendall"`.
#' @param epsilon \[`numeric(1)`: `0.05`\]\cr
#'   Threshold for the flat classification band around the neutral value.
#' @param turbulence_threshold \[`numeric(1)`: `5`\]\cr
#'   Baseline threshold for the combined volatility score.
#' @param flat_to_turbulent_factor \[`numeric(1)`: `1.5`\]\cr
#'   Multiplier applied to `turbulence_threshold` for points already
#'   classified as flat. Must be >= 1.
#' @param min_points \[`integer(1)`: `3L`\]\cr
#'   Minimum non-NA observations required in a window.
#' @param align \[`character(1)`: `"center"`\]\cr
#'   Window alignment: `"center"`, `"right"`, or `"left"`.
#' @return An object of class `"trend"` (a tibble) with columns:
#'
#' * `time`: the time index.
#' * `value`: the input time series values.
#' * `metric`: the rolling trend metric.
#' * `state`: factor with levels `"ascending"`, `"descending"`,
#'   `"flat"`, `"turbulent"`, `"Missing_Data"`, `"Initial"`.
#'
#' @seealso [compare_ts()] for visualizing trend classifications against
#'   the original series; [resilience()] for rolling resilience metrics.
#' @examples
#' \donttest{
#' set.seed(123)
#' x <- cumsum(rnorm(200, sd = 2))
#' tr <- compute_trend(x, window = 20)
#' tr
#' table(tr$state)
#' plot(tr)
#' }
compute_trend <- function(data, window = NULL, method = "slope",
                          slope_method = "robust", epsilon = 0.05,
                          turbulence_threshold = 5,
                          flat_to_turbulent_factor = 1.5,
                          min_points = 3L, align = "center") {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n  <- length(values)
  method <- check_match(method, c("slope", "ar1_phi1", "growth_factor"))
  slope_method <- check_match(
    slope_method,
    c("ols", "robust", "spearman", "kendall")
  )
  align <- check_match(align, c("center", "right", "left"))
  check_range(epsilon, type = "numeric", scalar = TRUE, min = 0, max = Inf)
  check_range(
    turbulence_threshold, type = "numeric", scalar = TRUE, min = 0, max = Inf
  )
  check_range(
    flat_to_turbulent_factor, type = "numeric", scalar = TRUE, min = 1, max = Inf
  )
  if (is.null(window)) {
    window <- max(3L, min(n, round(n / 10)))
  }
  check_range(window, type = "integer", min = 2L, max = n)
  if (min_points > window) {
    min_points <- window
  }
  metric <- rep(NA_real_, n)
  half_left  <- (window - 1L) %/% 2L
  half_right <- window - 1L - half_left
  iter_start <- 1L
  iter_end   <- n
  if (align == "center") {
    iter_start <- 1L + half_left
    iter_end   <- n - half_right
  } else if (align == "right") {
    iter_start <- window
  } else {
    iter_end <- n - window + 1L
  }
  for (k in iter_start:iter_end) {
    w_idx <- if (align == "center") {
      (k - half_left):(k + half_right)
    } else if (align == "right") {
      (k - window + 1L):k
    } else {
      k:(k + window - 1L)
    }
    wd <- values[w_idx]
    if (sum(!is.na(wd)) < min_points) {
      next
    }
    metric[k] <- switch(
      method,
      slope = trend_slope_(seq_along(wd), wd, slope_method, min_points),
      ar1_phi1 = trend_ar1_(wd, max(2L, min_points)),
      growth_factor = {
        clean <- wd[!is.na(wd)]
        if (length(clean) >= 2L && abs(clean[1L]) > 1e-10) {
          clean[length(clean)] / clean[1L]
        } else {
          NA_real_
        }
      }
    )
  }
  if (align == "center" && any(!is.na(metric))) {
    first_valid <- min(which(!is.na(metric)))
    if (first_valid > 1L) {
      metric[seq_len(first_valid - 1L)] <- metric[first_valid]
    }
    last_valid <- max(which(!is.na(metric)))
    if (last_valid < n) {
      metric[(last_valid + 1L):n] <- metric[last_valid]
    }
  }
  state <- rep("Initial", n)
  na_data <- is.na(values)
  state[na_data] <- "Missing_Data"
  neutral <- ifelse_(method %in% c("ar1_phi1", "growth_factor"), 1, 0)
  upper <- neutral + epsilon
  lower <- neutral - epsilon
  valid_idx <- !is.na(metric) & !na_data
  state[valid_idx] <- ifelse(
    metric[valid_idx] > upper, "ascending",
    ifelse(metric[valid_idx] < lower, "descending", "flat")
  )
  if (sum(valid_idx) >= min_points) {
    vol_win <- min(max(3L, window %/% 2L), sum(valid_idx))
    metrics_valid <- metric[valid_idx]
    orig_indices  <- which(valid_idx)
    if (length(metrics_valid) >= vol_win) {
      for (j in vol_win:length(metrics_valid)) {
        seg <- metrics_valid[(j - vol_win + 1L):j]
        if (sum(!is.na(seg)) < max(2L, min_points - 1L)) {
          next
        }
        combined <- trend_volatility_(seg)
        if (is.na(combined)) {
          next
        }
        idx_orig   <- orig_indices[j]
        base_trend <- state[idx_orig]
        eff_thresh <- ifelse_(
          base_trend == "flat",
          turbulence_threshold * flat_to_turbulent_factor,
          turbulence_threshold
        )
        if (base_trend != "Missing_Data" && combined > eff_thresh) {
          state[idx_orig] <- "turbulent"
        }
      }
    }
  }
  state <- factor(
    state,
    levels = c(
      "ascending", "descending", "flat",
      "turbulent", "Missing_Data", "Initial"
    )
  )
  out <- tibble::tibble(
    time = time,
    value = values,
    metric = metric,
    state = state
  )
  structure(
    out,
    window = window,
    align = align,
    method = method,
    slope_method = if (method == "slope") slope_method else NULL,
    epsilon = epsilon,
    class = c("trend", "tbl_df", "tbl", "data.frame")
  )
}

#' AR(1) coefficient for a single window
#'
#' Computes the first-order autoregressive coefficient via OLS, bounded to
#' (-1, 1). Returns `NA_real_` when data is insufficient or variance is
#' near zero.
#'
#' @param window_data Numeric vector.
#' @param min_points Minimum non-NA points required.
#' @return Single numeric, or `NA_real_`.
#' @noRd
trend_ar1_ <- function(window_data, min_points = 3L) {
  clean <- stats::na.omit(as.numeric(window_data))
  if (length(clean) < min_points) {
    return(NA_real_)
  }
  if (stats::var(clean, na.rm = TRUE) < 1e-10) {
    return(NA_real_)
  }
  model <- try_(
    suppressWarnings(
      stats::ar.ols(
        clean, aic = FALSE, order.max = 1L,
        demean = FALSE, intercept = TRUE
      )
    )
  )
  if (inherits(model, "try-error")) {
    model <- NULL
  }
  if (is.null(model) || length(model$ar) == 0L) {
    return(NA_real_)
  }
  max(-0.999, min(0.999, model$ar[1L]))
}

trend_slope_ <- function(x_vals, y_vals, slope_method, min_points) {
  valid <- !is.na(x_vals) & !is.na(y_vals)
  if (sum(valid) < min_points) return(NA_real_)
  xc <- x_vals[valid]
  yc <- y_vals[valid]
  if (length(unique(xc)) < 2L && slope_method != "robust") {
    return(NA_real_)
  }
  if (slope_method == "robust" && length(xc) < 2L) {
    return(NA_real_)
  }
  if (stats::var(xc, na.rm = TRUE) < 1e-10) {
    return(ifelse_(stats::var(yc, na.rm = TRUE) < 1e-10, 0, NA_real_))
  }
  if (stats::var(yc, na.rm = TRUE) < 1e-10) {
    return(0)
  }
  switch(
    slope_method,
    ols = stats::cov(xc, yc) / stats::var(xc),
    robust = {
      np <- length(xc)
      idx <- which(upper.tri(matrix(0, np, np)), arr.ind = TRUE)
      dx <- xc[idx[, 2L]] - xc[idx[, 1L]]
      dy <- yc[idx[, 2L]] - yc[idx[, 1L]]
      keep <- dx != 0
      if (!any(keep)) {
        return(NA_real_)
      }
      stats::median(dy[keep] / dx[keep], na.rm = TRUE)
    },
    spearman = ,
    kendall = {
      cor_val <- try_(
        stats::cor(xc, yc, method = slope_method, use = "complete.obs")
      )
      if (inherits(cor_val, "try-error")) {
        cor_val <- NA_real_
      }
      sd_x <- stats::sd(xc, na.rm = TRUE)
      if (is.na(cor_val) || is.na(sd_x) || sd_x < 1e-10) {
        NA_real_
      } else {
        cor_val * (stats::sd(yc, na.rm = TRUE) / sd_x)
      }
    }
  )
}

trend_volatility_ <- function(seg) {
  m_sd   <- stats::sd(seg, na.rm = TRUE)
  m_mean <- abs(mean(seg, na.rm = TRUE))
  m_rng  <- suppressWarnings(diff(range(seg, na.rm = TRUE)))
  if (is.na(m_sd) || is.na(m_mean) || is.na(m_rng)) {
    return(NA_real_)
  }
  m_sd / (m_mean + 1e-8) + 0.5 * m_rng / (m_mean + 1e-8)
}
