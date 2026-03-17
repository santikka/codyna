#' Parameter Sensitivity Analysis for Early Warning Signals
#'
#' Evaluates how the Kendall tau trend statistic of an EWS metric varies
#' across different rolling window sizes and detrending methods. This
#' reveals whether a detected signal is robust to parameter choices or an
#' artifact of a specific configuration.
#'
#' @details
#' The analysis proceeds in four steps:
#'
#' **1. Data preparation.** The input time series is parsed via
#' `prepare_timeseries_data()`. If `windows` is `NULL`, a sequence of 8
#' window sizes is generated spanning from `max(20, n/20)` to
#' `min(n/2, 200)` on a linear grid, rounded to integers. This range
#' covers small windows (high temporal resolution but noisy estimates)
#' through large windows (smooth estimates but low resolution).
#'
#' **2. Detrending.** For each detrend method, the full series is
#' preprocessed: `"none"` uses the raw series; `"linear"` removes the
#' best-fit linear trend via OLS residuals. Detrending before computing
#' EWS metrics eliminates confounding from deterministic trends that
#' could inflate autocorrelation or variance artificially.
#'
#' **3. Rolling metric computation.** For each (window, detrend)
#' combination, the specified metric is computed on right-aligned rolling
#' windows. The available metrics are:
#'
#'   * `"ar1"`: Lag-1 autocorrelation estimated via OLS AR(1).
#'     Rising AR(1) is a classic signature of critical slowing down
#'     (Scheffer et al., 2009).
#'   * `"sd"`: Standard deviation. Increasing variance is the
#'     second canonical EWS indicator alongside AR(1).
#'   * `"variance"`: Variance (SD squared). Equivalent information
#'     to SD but on the original scale, which may be preferable for
#'     some comparisons.
#'   * `"skewness"`: Third standardized moment. Non-zero skewness
#'     indicates asymmetric fluctuations, which can precede transitions
#'     in systems with asymmetric potential wells.
#'   * `"kurtosis"`: Excess kurtosis (fourth standardized moment
#'     minus 3). Heavy-tailed fluctuations (positive kurtosis) can
#'     signal proximity to a transition.
#'   * `"cv"`: Coefficient of variation (SD / |mean|). Useful when
#'     the mean is non-zero and relative dispersion matters more than
#'     absolute dispersion.
#'
#' **4. Kendall tau trend test.** For each (window, detrend)
#' combination, Kendall's rank correlation between the metric values and
#' their time indices is computed. A consistently positive tau across
#' many parameter combinations provides robust evidence for an upward
#' trend in the EWS metric (and thus for approaching a critical
#' transition). Sensitivity is assessed by examining whether tau remains
#' consistently positive (or negative) regardless of window size and
#' detrending.
#'
#' @export
#' @param data \[`numeric` or `ts`]\cr
#'   The time series to analyze. Accepts a numeric vector or a `ts` object.
#' @param metric \[`character(1)`: `"ar1"`]\cr
#'   EWS metric to test. One of `"ar1"`, `"sd"`, `"variance"`,
#'   `"skewness"`, `"kurtosis"`, `"cv"`.
#' @param windows \[`integer()`: `NULL`]\cr
#'   Vector of window sizes to evaluate. If `NULL`, a default sequence
#'   of 8 sizes is generated spanning `max(20, n/20)` to
#'   `min(n/2, 200)`.
#' @param detrend_methods \[`character()`: `c("none", "linear")`]\cr
#'   Detrending methods to compare. `"none"` uses the raw series;
#'   `"linear"` removes a linear trend via OLS residuals.
#' @param method \[`character(1)`: `"rolling"`]\cr
#'   Analysis method. Currently only `"rolling"` is supported.
#' @return An object of class `"sensitivity_ews"` (inheriting from
#'   [tibble::tbl_df]) with the following columns:
#'
#'     * `window`: integer window size used.
#'     * `detrend`: character detrend method applied.
#'     * `time`: time index for each metric value.
#'     * `score`: computed metric value at that time point.
#'     * `tau`: Kendall's tau for this (window, detrend)
#'       combination.
#'
#'   Attributes stored on the object:
#'
#'   * `metric`: the EWS metric name.
#'   * `windows`: integer vector of window sizes evaluated.
#'   * `detrend_methods`: character vector of detrend methods.
#'
#' @references
#' Scheffer, M., Bascompte, J., Brock, W.A., Brovkin, V., Carpenter,
#' S.R., Dakos, V., Held, H., van Nes, E.H., Rietkerk, M., &
#' Sugihara, G. (2009). Early-warning signals for critical transitions.
#' \emph{Nature}, 461, 53--59. \doi{10.1038/nature08227}
#'
#' Dakos, V., Carpenter, S.R., Brock, W.A., Ellison, A.M., Guttal, V.,
#' Ives, A.R., Kefi, S., Livina, V., Seekell, D.A., van Nes, E.H., &
#' Scheffer, M. (2012). Methods for detecting early warnings of
#' critical transitions in time series illustrated using simulated
#' ecological data. \emph{PLoS ONE}, 7(7), e41010.
#' \doi{10.1371/journal.pone.0041010}
#'
#' Lenton, T.M. (2011). Early warning of climate tipping points.
#' \emph{Nature Climate Change}, 1(4), 201--209.
#' \doi{10.1038/nclimate1143}
#'
#' @seealso [resilience()] for rolling resilience metric computation;
#'   [compute_trend()] for trend classification; [hurst()] for
#'   Hurst exponent analysis.
#' @examples
#' \donttest{
#' set.seed(42)
#' # Generate a series with increasing variance (approaching transition)
#' x <- cumsum(rnorm(300, sd = seq(0.5, 2, length.out = 300)))
#' sa <- sensitivity_ews(x, metric = "ar1")
#' sa
#' summary(sa)
#' plot(sa, type = "heatmap")
#' plot(sa, type = "lines")
#' }
sensitivity_ews <- function(data, metric = "ar1", windows = NULL,
                            detrend_methods = c("none", "linear"),
                            method = "rolling") {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)
  metric <- check_match(
    metric, c("ar1", "sd", "variance", "skewness", "kurtosis", "cv")
  )
  detrend_methods <- check_match(
    detrend_methods, c("none", "linear"), several.ok = TRUE
  )
  method <- check_match(method, c("rolling"))
  if (is.null(windows)) {
    w_min <- max(20L, round(n / 20))
    w_max <- min(round(n / 2), 200L)
    if (w_min >= w_max) {
      windows <- unique(as.integer(round(
        seq(w_min, w_max, length.out = 2L)
      )))
    } else {
      windows <- unique(
        as.integer(
          round(
            seq(w_min, w_max, length.out = 8L)
          )
        )
      )
    }
  } else {
    windows <- as.integer(round(windows))
  }
  # Validate each window size
  windows <- windows[windows >= 2L & windows <= n]
  stopifnot_(
    length(windows) > 0L,
    "No valid window sizes. All values must be between 2 and {n}."
  )
  # Compute metrics for all (window, detrend) combinations
  results <- vector("list", length(windows) * length(detrend_methods))
  idx <- 0L
  for (dm in detrend_methods) {
    detrended <- sensitivity_detrend_(values, dm)
    for (w in windows) {
      idx <- idx + 1L
      rolling <- sensitivity_rolling_(detrended, w, metric)
      if (length(rolling$score) == 0L) {
        results[[idx]] <- tibble::tibble(
          window = integer(0),
          detrend = character(0),
          time = numeric(0),
          score = numeric(0),
          tau = numeric(0)
        )
        next
      }
      tau_val <- sensitivity_tau_(rolling$score)
      results[[idx]] <- tibble::tibble(
        window = rep(as.integer(w), length(rolling$score)),
        detrend = rep(dm, length(rolling$score)),
        time = time[rolling$time],
        score = rolling$score,
        tau = rep(tau_val, length(rolling$score))
      )
    }
  }
  out <- do.call(rbind, results)
  out <- tibble::as_tibble(out)
  structure(
    out,
    metric = metric,
    windows = windows,
    detrend_methods = detrend_methods,
    class = c("sensitivity_ews", "tbl_df", "tbl", "data.frame")
  )
}

#' Detrend a time series
#'
#' Removes a linear trend via OLS residuals. Returns the original series
#' when `method = "none"`.
#'
#' @param values Numeric vector.
#' @param method One of `"none"` or `"linear"`.
#' @return Numeric vector of same length.
#' @noRd
sensitivity_detrend_ <- function(values, method) {
  if (method == "none") {
    return(values)
  }
  n <- length(values)
  time_idx <- seq_len(n)
  valid <- !is.na(values)
  if (sum(valid) < 3L) return(values)
  fit <- try_(stats::lm(values[valid] ~ time_idx[valid]))
  if (inherits(fit, "try-error")) {
    return(values)
  }
  out <- values
  out[valid] <- stats::residuals(fit)
  out
}

#' Compute a single EWS metric on a numeric window
#'
#' @param x Numeric vector (window data, no NAs expected).
#' @param metric One of "ar1", "sd", "variance", "skewness", "kurtosis", "cv".
#' @return Single numeric value.
#' @noRd
sensitivity_metric_ <- function(x, metric) {
  n <- length(x)
  if (n < 3L) {
    return(NA_real_)
  }
  switch(
    metric,
    ar1 = {
      if (stats::var(x, na.rm = TRUE) < 1e-10) {
        return(NA_real_)
      }
      model <- try_(
        stats::ar.ols(
          x,
          aic = FALSE,
          order.max = 1L,
          demean = FALSE,
          intercept = TRUE
        )
      )
      if (inherits(model, "try-error") || is.null(model) ||
          length(model$ar) == 0L) {
        return(NA_real_)
      }
      max(-0.999, min(0.999, model$ar[1L]))
    },
    sd = stats::sd(x, na.rm = TRUE),
    variance = stats::var(x, na.rm = TRUE),
    skewness = {
      mu <- mean(x, na.rm = TRUE)
      s <- stats::sd(x, na.rm = TRUE)
      if (is.na(s) || s < 1e-10) {
        return(NA_real_)
      }
      mean(((x - mu) / s)^3, na.rm = TRUE)
    },
    kurtosis = {
      mu <- mean(x, na.rm = TRUE)
      s <- stats::sd(x, na.rm = TRUE)
      if (is.na(s) || s < 1e-10) {
        return(NA_real_)
      }
      mean(((x - mu) / s)^4, na.rm = TRUE) - 3
    },
    cv = {
      mu <- mean(x, na.rm = TRUE)
      if (abs(mu) < .Machine$double.eps) {
        return(NA_real_)
      }
      stats::sd(x, na.rm = TRUE) / abs(mu)
    }
  )
}

#' Compute rolling metric values for a detrended series
#'
#' @param values Numeric vector (already detrended).
#' @param window Integer window size.
#' @param metric Character metric name.
#' @return List with `time` (integer indices) and `score` (metric values).
#' @noRd
sensitivity_rolling_ <- function(values, window, metric) {
  n <- length(values)
  n_windows <- n - window + 1L
  if (n_windows < 1L) {
    return(
      list(
        time = integer(0),
        score = numeric(0)
      )
    )
  }
  scores <- vapply(seq_len(n_windows), function(i) {
    wd <- values[i:(i + window - 1L)]
    clean <- wd[!is.na(wd)]
    if (length(clean) < 3L) {
      return(NA_real_)
    }
    sensitivity_metric_(clean, metric)
  }, numeric(1L))
  # Right-aligned: metric at position window, window+1, ..., n
  list(
    time = seq(window, n),
    score = scores
  )
}

#' Compute Kendall tau of a metric series against time
#'
#' @param scores Numeric vector of metric values.
#' @return Single numeric Kendall tau value, or `NA_real_`.
#' @noRd
sensitivity_tau_ <- function(scores) {
  valid <- !is.na(scores)
  if (sum(valid) < 4L) {
    return(NA_real_)
  }
  time_idx <- seq_along(scores)[valid]
  vals <- scores[valid]
  result <- try_(
    suppressWarnings(
      stats::cor.test(time_idx, vals, method = "kendall")
    )
  )
  if (inherits(result, "try-error")) {
    return(NA_real_)
  }
  result$estimate
}

#' Color palette for sensitivity heatmap
#' @noRd
sensitivity_colors_ <- function() {
  c(
    "#2166AC",
    "#4393C3",
    "#92C5DE",
    "#D1E5F0",
    "#F7F7F7",
    "#FDDBC7",
    "#F4A582",
    "#D6604D",
    "#B2182B"
  )
}
