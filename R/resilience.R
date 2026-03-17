#' Resilience Metrics for Time Series Data
#'
#' Computes rolling-window resilience metrics for univariate time series data,
#' quantifying stability, volatility, variability, memory, complexity,
#' and recovery characteristics. Metric computations are identical to
#' `calculate_rolling_resilience_metrics()` from the Resilience4 toolkit.
#'
#' @details
#' The function slides a fixed-size window across the time series. Within each
#' window position, eight resilience-relevant metrics are computed on the local
#' segment. The resulting metric trajectories reveal *when* and *how* resilience
#' changes over time, not just whether the system is resilient on average.
#'
#' The eight metrics and their mathematical definitions are:
#'
#' 1. `vsi` (Variance of Sub-window Variances) The variance of rolling
#'     variances computed on sub-windows within each main window
#'     (\eqn{\mathrm{Var}(\sigma^2_{w_1}, \ldots, \sigma^2_{w_k})}). Captures
#'     instability in volatility. Lower values indicate more stable variance
#'     and thus better resilience.
#' 2. `arch_lm` (ARCH-LM Statistic) The Engle (1982) test statistic
#'     \eqn{n \cdot R^2} from regressing squared residuals on their own lags.
#'     Detects autoregressive conditional heteroscedasticity (volatility
#'     clustering). Lower values indicate more homogeneous variance and better
#'     resilience.
#' 3. `cv` (Coefficient of Variation) \eqn{\mathrm{CV} = \sigma / |\mu|},
#'     the ratio of standard deviation to the absolute mean. Measures relative
#'     dispersion. Lower values indicate tighter fluctuations around the mean
#'     and better resilience.
#'  4. `recovery_time` Median number of time steps required to return within
#'     a threshold of the baseline after a shock event. Shocks are detected as
#'     deviations exceeding `shock_threshold` baseline standard deviations.
#'     Lower values indicate faster recovery and better resilience.
#'  5. `recovery_slope` Median OLS slope of the series over a short window
#'     immediately following each detected shock. Positive slopes indicate
#'     upward recovery toward baseline; the sign and magnitude together convey
#'     recovery direction and speed. Higher (more positive) values indicate
#'     better resilience.
#'  6. `sample_entropy` \eqn{\mathrm{SampEn} = -\ln(A / B)}, where \eqn{A}
#'     is the number of template matches at embedding dimension \eqn{m + 1} and
#'     \eqn{B} is the number of matches at dimension \eqn{m} (Richman &
#'     Moorman, 2000). Quantifies signal complexity and regularity. Higher
#'     values indicate greater unpredictability and better resilience (overly
#'     regular dynamics may signal rigidity).
#'  7. `dfa_alpha` (DFA Scaling Exponent) The slope of the log-log
#'     fluctuation plot from Detrended Fluctuation Analysis (Peng et al.,
#'     1994). Values near 0.5 indicate uncorrelated noise; values near 1.0
#'     indicate \eqn{1/f} noise; values near 1.5 indicate Brownian motion.
#'     Lower values (closer to white noise) indicate better resilience.
#'  8. `ac_ratio` (Autocorrelation Ratio) \eqn{|\mathrm{acf}[1]| /
#'     \mathrm{mean}(|\mathrm{acf}[2{:}k]|)}, the ratio of lag-1
#'     autocorrelation magnitude to the mean magnitude of higher-order lags.
#'     High values suggest strong short-term memory with faster decorrelation.
#'     Higher values indicate better resilience.
#'
#' Note that metrics are directional: for some (vsi, arch_lm, cv,
#' recovery_time, dfa_alpha) lower values signal better resilience, while for
#' others (recovery_slope, sample_entropy, ac_ratio) higher values are
#' preferable. The [classify_resilience()] function handles this directionality
#' automatically.
#'
#' @export
#' @param data \[`ts`, `numeric()`, `data.frame`]\cr
#'   Univariate time series data.
#' @param window \[`integer(1)`: `50L`]\cr
#'   A positive `integer` specifying the rolling window size.
#'   Must be at least `20`.
#' @param align \[`character(1)`: `"right"`]\cr
#'   Alignment of the window. The available options are:
#'   `"right"`, `"center"`, and `"left"`.
#' @param metrics \[`character()`]\cr
#'   A vector of metrics to calculate, or `"all"` for all eight metrics.
#' @param shock_threshold \[`numeric(1)`: `2`]\cr
#'   Number of baseline standard deviations for shock detection.
#' @param recovery_threshold \[`numeric(1)`: `0.1`]\cr
#'   Fraction of baseline SD within which recovery is declared.
#' @param recovery_window \[`integer(1)`: `5L`]\cr
#'   Number of points after a shock to fit recovery slope.
#' @param baseline_window \[`integer(1)`: `10L`]\cr
#'   Number of initial points used to establish baseline.
#' @param edim \[`integer(1)`: `2L`]\cr
#'   Embedding dimension for sample entropy.
#' @param r \[`numeric(1)`: `NULL`]\cr
#'   Tolerance for sample entropy template matching.
#'   If `NULL`, uses `0.2 * sd(x)`.
#' @param tau \[`integer(1)`: `1L`]\cr
#'   Time delay for sample entropy embedding.
#' @param max_lag \[`integer(1)`: `5L`]\cr
#'   Maximum lag for autocorrelation ratio.
#' @param lags \[`integer(1)`: `1L`]\cr
#'   Number of lags for ARCH-LM test.
#' @param demean \[`logical(1)`: `FALSE`]\cr
#'   Whether to demean the series before ARCH-LM computation.
#' @return An object of class `"resilience"` (inherits from
#'   [tibble::tbl_df]) with the following columns:

#'   * `time`: integer or Date index of each observation.
#'   * `value`: the original time series values.
#'   * `vsi`: Variance of Sub-window Variances (if requested).
#'   * `arch_lm`: ARCH-LM statistic (if requested).
#'   * `cv`: Coefficient of Variation (if requested).
#'   * `recovery_time`: median recovery time in steps (if requested).
#'   * `recovery_slope`: median post-shock OLS slope (if requested).
#'   * `sample_entropy`: Sample Entropy (if requested).
#'   * `dfa_alpha`: DFA scaling exponent (if requested).
#'   * `ac_ratio`: Autocorrelation Ratio (if requested).
#'
#'   The returned object carries three attributes:
#'
#'   * `window`: the window size used.
#'   * `align`: the alignment method (`"right"`, `"center"`, or
#'       `"left"`).
#'   * `metrics`: character vector of computed metric names.
#'
#'   Leading or trailing rows contain `NA` values for the metrics depending on
#'   the alignment setting, because the rolling window cannot be computed at
#'   those positions.
#'
#' @references
#' Richman, J. S. & Moorman, J. R. (2000). Physiological time-series analysis
#' using approximate entropy and sample entropy. *American Journal of
#' Physiology-Heart and Circulatory Physiology*, 278(6), H2039--H2049.
#' \doi{10.1152/ajpheart.2000.278.6.H2039}
#'
#' Peng, C.-K., Buldyrev, S. V., Havlin, S., Simons, M., Stanley, H. E., &
#' Goldberger, A. L. (1994). Mosaic organization of DNA nucleotides.
#' *Physical Review E*, 49(2), 1685--1689.
#' \doi{10.1103/PhysRevE.49.1685}
#'
#' Engle, R. F. (1982). Autoregressive Conditional Heteroscedasticity with
#' Estimates of the Variance of United Kingdom Inflation. *Econometrica*,
#' 50(4), 987--1007. \doi{10.2307/1912773}
#'
#' @seealso [classify_resilience()] for directional normalization and
#'   categorical scoring; [plot.resilience()] for visualization.
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' res <- resilience(x, window = 50L)
#' res
#'
#' # Visualize raw metric trajectories
#' plot(res, type = "lines")
#' }
resilience <- function(data, window = 50L, align = "right", metrics = "all",
                       shock_threshold = 2, recovery_threshold = 0.1,
                       recovery_window = 5L, baseline_window = 10L,
                       edim = 2L, r = NULL, tau = 1L, max_lag = 5L,
                       lags = 1L, demean = FALSE) {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)
  align <- check_match(align, c("right", "center", "left"))
  available_metrics <- c(
    "vsi", "arch_lm", "cv", "recovery_time", "recovery_slope",
    "sample_entropy", "dfa_alpha", "ac_ratio"
  )
  metrics <- check_match(
    metrics,
    c(available_metrics, "all"),
    several.ok = TRUE
  )
  metrics <- ifelse_("all" %in% metrics, available_metrics, metrics)
  check_range(window, type = "integer", min = 20L, max = n)
  total_windows <- n - window + 1L
  results_list <- vector("list", total_windows)
  for (i in seq_len(total_windows)) {
    window_data <- values[i:(i + window - 1L)]
    window_data <- window_data[!is.na(window_data)]
    wd_len <- length(window_data)
    if (wd_len < (window / 2)) {
      row <- stats::setNames(
        as.list(rep(NA_real_, length(metrics))),
        metrics
      )
    } else {
      adaptive_min_scale <- max(3L, wd_len %/% 20L)
      adaptive_max_scale <- max(adaptive_min_scale + 2L, wd_len %/% 4L)
      row <- list()
      for (m in metrics) {
        row[[m]] <- switch(
          m,
          vsi = resilience_vsi(window_data),
          arch_lm = resilience_arch_lm(window_data, lags, demean),
          cv = resilience_cv(window_data),
          recovery_time = resilience_recovery_time(
            window_data, shock_threshold, recovery_threshold, baseline_window
          ),
          recovery_slope = resilience_recovery_slope(
            window_data, shock_threshold, recovery_window, baseline_window
          ),
          sample_entropy = resilience_sample_entropy(
            window_data, edim, r, tau
          ),
          dfa_alpha = resilience_dfa_alpha(
            window_data,
            min_scale = adaptive_min_scale,
            max_scale = adaptive_max_scale,
            min_obs = max(15L, wd_len %/% 3L),
            min_r_squared = 0.3
          ),
          ac_ratio = resilience_ac_ratio(window_data, max_lag)
        )
      }
    }
    results_list[[i]] <- row
  }
  # Convert to data frame
  metrics_df <- do.call(
    rbind,
    lapply(results_list, function(r) as.data.frame(r, stringsAsFactors = FALSE))
  )
  pad_start <- 0L
  pad_end <- 0L
  if (align == "right") {
    pad_start <- window - 1L
  } else if (align == "center") {
    pad_start <- floor((window - 1L) / 2)
    pad_end <- ceiling((window - 1L) / 2)
  } else if (align == "left") {
    pad_end <- window - 1L
  }
  if (pad_start > 0L) {
    na_pad <- data.frame(
      matrix(NA_real_, nrow = pad_start, ncol = ncol(metrics_df))
    )
    names(na_pad) <- names(metrics_df)
    metrics_df <- rbind(na_pad, metrics_df)
  }
  if (pad_end > 0L) {
    na_pad <- data.frame(
      matrix(NA_real_, nrow = pad_end, ncol = ncol(metrics_df))
    )
    names(na_pad) <- names(metrics_df)
    metrics_df <- rbind(metrics_df, na_pad)
  }
  out <- tibble::as_tibble(
    cbind(data.frame(time = time, value = values), metrics_df)
  )
  structure(
    out,
    window = window,
    align = align,
    metrics = metrics,
    class = c("resilience", "tbl_df", "tbl", "data.frame")
  )
}

#' @noRd
resilience_vsi <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 10L) return(NA_real_)
  window_size <- max(5L, floor(n / 4))
  if (n < window_size * 2L) {
    return(stats::var(x, na.rm = TRUE))
  }
  n_windows <- n - window_size + 1L
  rolling_vars <- numeric(n_windows)
  for (i in seq_len(n_windows)) {
    rolling_vars[i] <- stats::var(x[i:(i + window_size - 1L)], na.rm = TRUE)
  }
  stats::var(rolling_vars, na.rm = TRUE)
}

resilience_arch_lm <- function(x, lags = 1L, demean = FALSE) {
  x <- as.numeric(x[!is.na(x)])
  n <- length(x)
  if (n < max(10L, lags + 3L)) {
    return(NA_real_)
  }
  if (demean) {
    x <- x - mean(x)
  }
  x2 <- x^2
  X <- matrix(NA_real_, nrow = n, ncol = lags)
  for (i in seq_len(lags)) {
    X[(i + 1L):n, i] <- x2[1L:(n - i)]
  }
  y <- x2[(lags + 1L):n]
  X_clean <- X[(lags + 1L):n, , drop = FALSE]
  if (length(y) < 3L) {
    return(NA_real_)
  }
  fit <- try_(stats::lm(y ~ X_clean))
  if (inherits(fit, "try-error")) {
    return(NA_real_)
  }
  r_sq <- suppressWarnings(summary(fit)$r.squared)
  length(y) * r_sq
}

resilience_cv <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3L) {
    return(NA_real_)
  }
  mu <- mean(x)
  if (abs(mu) < .Machine$double.eps) {
    return(NA_real_)
  }
  stats::sd(x) / abs(mu)
}

#' @noRd
resilience_recovery_time <- function(x, shock_threshold = 2,
                                     recovery_threshold = 0.1,
                                     baseline_window = 10L) {
  x <- x[!is.na(x) & is.finite(x)]
  n <- length(x)
  if (n < 20L) return(NA_real_)
  baseline_data <- x[seq_len(min(baseline_window, n))]
  if (length(baseline_data) < 3L) {
    return(NA_real_)
  }
  baseline <- mean(baseline_data)
  baseline_sd <- stats::sd(baseline_data)
  if (is.na(baseline_sd) || baseline_sd == 0) {
    return(NA_real_)
  }
  deviations <- abs(x - baseline)
  shock_points <- which(deviations > shock_threshold * baseline_sd)
  if (length(shock_points) == 0L) {
    return(0)
  }
  recovery_times <- numeric(length(shock_points))
  for (i in seq_along(shock_points)) {
    shock_idx <- shock_points[i]
    recovery_threshold_value <- recovery_threshold * baseline_sd
    if (shock_idx < n) {
      remaining <- abs(x[(shock_idx + 1L):n] - baseline)
      recovered_indices <- which(remaining <= recovery_threshold_value)
      if (length(recovered_indices) > 0L) {
        recovery_times[i] <- recovered_indices[1L]
      } else {
        recovery_times[i] <- n - shock_idx
      }
    } else {
      recovery_times[i] <- 0
    }
  }
  stats::median(recovery_times, na.rm = TRUE)
}

resilience_recovery_slope <- function(x, shock_threshold = 2,
                                      recovery_window = 5L,
                                      baseline_window = 10L) {
  x <- x[!is.na(x) & is.finite(x)]
  n <- length(x)
  if (n < 20L) {
    return(NA_real_)
  }
  baseline_data <- x[seq_len(min(baseline_window, n))]
  if (length(baseline_data) < 3L) {
    return(NA_real_)
  }
  baseline <- mean(baseline_data)
  baseline_sd <- stats::sd(baseline_data)
  if (is.na(baseline_sd) || baseline_sd == 0) return(NA_real_)
  deviations <- abs(x - baseline)
  shock_points <- which(deviations > shock_threshold * baseline_sd)
  if (length(shock_points) == 0L) return(0)
  recovery_slopes <- numeric(length(shock_points))
  for (i in seq_along(shock_points)) {
    shock_idx <- shock_points[i]
    end_idx <- min(shock_idx + recovery_window, n)
    if (end_idx - shock_idx < 2L) {
      recovery_slopes[i] <- NA_real_
      next
    }
    time_points <- shock_idx:end_idx
    values <- x[time_points]
    valid_idx <- !is.na(values) & is.finite(values)
    if (sum(valid_idx) < 2L) {
      recovery_slopes[i] <- NA_real_
      next
    }
    fit <- try_(stats::lm(values[valid_idx] ~ time_points[valid_idx]))
    if (inherits(fit, "try-error")) {
      recovery_slopes[i] <- NA_real_
    } else {
      recovery_slopes[i] <- stats::coef(fit)[2L]
    }
  }
  stats::median(recovery_slopes, na.rm = TRUE)
}

resilience_sample_entropy <- function(x, edim = 2L, r = NULL, tau = 1L) {
  x <- as.numeric(x[!is.na(x)])
  n <- length(x)
  if (n < 10L) {
    return(NA_real_)
  }
  if (is.null(r)) {
    r <- 0.2 * stats::sd(x)
  }
  if (is.na(r) || r < 0 || !is.finite(r) || n <= edim * tau) {
    return(NA_real_)
  }
  # Create embedded vectors
  create_embedded <- function(data, m, delay) {
    nn <- length(data)
    num_vec <- nn - (m - 1L) * delay
    if (num_vec <= 0L) {
      return(matrix(nrow = 0L, ncol = m))
    }
    emb <- matrix(NA_real_, nrow = num_vec, ncol = m)
    for (i in seq_len(num_vec)) {
      indices <- i + (0L:(m - 1L)) * delay
      emb[i, ] <- data[indices]
    }
    emb
  }
  # Count template matches (Chebyshev distance <= r)
  count_matches <- function(vectors, tol) {
    nv <- nrow(vectors)
    if (nv < 2L) {
      return(0L)
    }
    matches <- 0L
    for (i in seq_len(nv - 1L)) {
      for (j in (i + 1L):nv) {
        if (max(abs(vectors[i, ] - vectors[j, ])) <= tol) {
          matches <- matches + 1L
        }
      }
    }
    matches
  }
  vectors_m <- create_embedded(x, edim, tau)
  vectors_m1 <- create_embedded(x, edim + 1L, tau)
  n_m1 <- nrow(vectors_m1)
  if (n_m1 < 2L) {
    return(NA_real_)
  }
  vectors_m_trunc <- vectors_m[seq_len(n_m1), , drop = FALSE]
  A <- count_matches(vectors_m1, r)
  B <- count_matches(vectors_m_trunc, r)
  -log(A / B)
}

resilience_dfa_alpha <- function(x, min_scale = 4L, max_scale = NULL,
                                 scales_per_octave = 8L, poly_order = 1L,
                                 min_obs = 30L, min_r_squared = 0.6) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < min_obs) return(NA_real_)
  if (is.null(max_scale)) max_scale <- floor(n / 6)
  if (min_scale >= max_scale) return(NA_real_)
  # Integrated profile
  profile <- cumsum(x - mean(x, na.rm = TRUE))
  # Log-spaced scales
  scales <- unique(round(
    exp(seq(log(min_scale), log(max_scale), length.out = 20L))
  ))
  scales <- scales[scales >= min_scale & scales <= max_scale]
  if (length(scales) < 4L) return(NA_real_)
  fluct <- numeric(length(scales))
  for (si in seq_along(scales)) {
    s <- scales[si]
    n_seg <- floor(n / s)
    if (n_seg < 3L) {
      fluct[si] <- NA_real_
      next
    }
    all_residuals <- c()
    for (seg in seq_len(n_seg)) {
      start_idx <- (seg - 1L) * s + 1L
      end_idx <- seg * s
      segment <- profile[start_idx:end_idx]
      tp <- seq_along(segment)
      fit <- try_(
        stats::lm(
          segment ~ stats::poly(tp, degree = poly_order, raw = TRUE)
        )
      )
      if (!inherits(fit, "try-error")) {
        all_residuals <- c(all_residuals, stats::residuals(fit))
      }
    }
    if (length(all_residuals) > 0L) {
      fluct[si] <- sqrt(mean(all_residuals^2))
    } else {
      fluct[si] <- NA_real_
    }
  }
  valid <- !is.na(fluct) & fluct > 0
  if (sum(valid) < 4L) return(NA_real_)
  log_s <- log10(scales[valid])
  log_f <- log10(fluct[valid])
  fit <- try_(stats::lm(log_f ~ log_s))
  if (inherits(fit, "try-error")) {
    return(NA_real_)
  }
  alpha <- stats::coef(fit)[2L]
  r_sq <- summary(fit)$r.squared
  if (r_sq < min_r_squared) {
    warning_(
      "Poor linear fit in DFA (R^2 = {round(r_sq, 3)})."
    )
  }
  alpha
}

resilience_ac_ratio <- function(x, max_lag = 5L) {
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 10L) return(NA_real_)
  if (n <= max_lag) return(NA_real_)
  ac_result <- try_(
    stats::acf(x, lag.max = max_lag, plot = FALSE, na.action = stats::na.pass)
  )
  if (inherits(ac_result, "try-error")) {
    return(NA_real_)
  }
  ac_vals <- as.numeric(ac_result$acf)[-1L]
  if (length(ac_vals) < 2L) {
    return(NA_real_)
  }
  ac_lag1 <- ac_vals[1L]
  ac_higher <- mean(abs(ac_vals[-1L]), na.rm = TRUE)
  if (is.na(ac_higher) || ac_higher == 0) {
    return(Inf)
  }
  abs(ac_lag1) / ac_higher
}

#' Classify Resilience Scores
#'
#' Applies directional normalization to resilience metrics and computes
#' a weighted composite score with category labels.
#'
#' @details
#' Raw metric values from [resilience()] are not directly comparable because
#' they live on different scales and some are "lower is better" while others
#' are "higher is better." This function resolves both problems in three
#' steps:
#'
#' **1. Directional normalization.** Each metric is rescaled to a \[0, 1\]
#' score where 0 means "Excellent resilience" and 1 means "Troubled
#' resilience," regardless of the original direction. Metrics where lower raw
#' values indicate better resilience (vsi, arch_lm, cv, recovery_time,
#' dfa_alpha) are scaled so that low raw values map to scores near 0.
#' Metrics where higher raw values indicate better resilience
#' (recovery_slope, sample_entropy, ac_ratio) are inverted accordingly.
#'
#' The `"empirical_state_aware"` method (the default) derives scaling
#' thresholds from empirically observed means for known resilience phases
#' (Stable, Recovery, Shock, Turbulent, Volatile). "Excellent" and
#' "Troubled" anchor points are set from quantiles of the good-state and
#' bad-state phase means, respectively. Alternative methods
#' (`"empirical_percentile"`, `"percentile"`, `"minmax"`, `"zscore"`,
#' `"divide_max"`) use purely data-driven scaling.
#'
#' **2. Smoothing.** After scaling, each metric's score trajectory is
#' smoothed with a rolling median (or mean) of width `smooth_window` to
#' reduce transient noise in the classification. Set `smooth_window = 1L`
#' to disable smoothing.
#'
#' **3. Weighted composite.** The per-metric scores are combined into a
#' single `composite_score` using empirically derived weights that reflect
#' each metric's importance for distinguishing resilience phases.
#' The default weights are: vsi = 0.256, cv = 0.226, arch_lm = 0.188,
#' ac_ratio = 0.070, recovery_slope = 0.058, dfa_alpha = 0.036,
#' sample_entropy = 0.028, recovery_time = 0.016.
#'
#' Both per-metric and composite scores are binned into six ordered
#' categories:
#'
#' * \strong{Excellent} (0.00--0.15)
#' * \strong{Solid} (0.15--0.30)
#' * \strong{Fair} (0.30--0.50)
#' * \strong{Vulnerable} (0.50--0.70)
#' * \strong{Failing} (0.70--0.85)
#' * \strong{Troubled} (0.85--1.00)
#'
#' @export
#' @param data \[`resilience`\]\cr
#'   An object of class `resilience` as returned by [resilience()].
#' @param method \[`character(1)`: `"empirical_state_aware"`\]\cr
#'   Normalization method for scaling raw metric values. The available
#'   options are:
#'   * `"empirical_state_aware"`: State-aware scaling using empirical
#'     thresholds derived from known resilience phases.
#'   * `"empirical_percentile"`: Hardcoded percentile-bin scaling.
#'   * `"percentile"`: Empirical CDF (rank-based) scaling.
#'   * `"minmax"`: Min-max scaling to \[0, 1\].
#'   * `"zscore"`: Z-score standardization.
#'   * `"divide_max"`: Division by absolute maximum.
#' @param weights \[`numeric()`: `NULL`\]\cr
#'   Named numeric vector of metric weights for the composite score.
#'   If `NULL`, empirically derived default weights are used.
#' @param smooth_window \[`integer(1)`: `10L`\]\cr
#'   Window size for post-scaling median smoothing (applied per-metric).
#'   Set to `1L` to disable smoothing.
#' @param smooth_method \[`character(1)`: `"median"`\]\cr
#'   Smoothing method: `"median"` or `"mean"`.
#' @return The input `"resilience"` object augmented with additional columns:
#'
#'     * `{metric}_score`: directionally normalized score in \[0, 1\]
#'       for each computed metric (e.g., `vsi_score`, `cv_score`).
#'     * `{metric}_category`: ordered factor label
#'       (`"Excellent"` through `"Troubled"`) for each metric.
#'     * `composite_score`: weighted average of all per-metric scores.
#'     * `composite_category`: ordered factor label for the composite
#'       score.
#'
#'   All original columns (`time`, `value`, and the raw metric columns) are
#'   preserved. The returned object retains the `"resilience"` class and
#'   its attributes.
#'
#' @seealso [resilience()] for computing the raw metrics;
#'   [plot.resilience()] for visualization of classified results.
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' res <- resilience(x, window = 50L)
#' cls <- classify_resilience(res)
#' cls
#'
#' # Ribbon plot with classified scores
#' plot(cls)
#' }
classify_resilience <- function(data, method = "empirical_state_aware",
                                weights = NULL, smooth_window = 10L,
                                smooth_method = "median") {
  check_missing(data)
  check_class(data, "resilience")
  method <- check_match(
    method,
    c(
      "empirical_state_aware",
      "empirical_percentile",
      "percentile",
      "minmax",
      "zscore",
      "divide_max"
    )
  )
  smooth_method <- check_match(smooth_method, c("median", "mean"))
  metrics <- attr(data, "metrics")
  default_w <- default_resilience_weights()
  if (is.null(weights)) {
    weights <- default_w[metrics]
  } else {
    stopifnot_(
      is.numeric(weights) && length(weights) == length(metrics),
      "Argument {.arg weights} must be a named numeric vector of
      length {length(metrics)}."
    )
    if (is.null(names(weights))) names(weights) <- metrics
  }
  # Category breaks and labels
  breaks <- c(0, 0.15, 0.30, 0.50, 0.70, 0.85, 1.0)
  labels <- c(
    "Excellent", "Solid", "Fair", "Vulnerable", "Failing", "Troubled"
  )
  # Scale each metric directionally
  n <- nrow(data)
  weighted_mat <- matrix(0, nrow = n, ncol = length(metrics))
  total_weight <- 0
  for (mi in seq_along(metrics)) {
    m <- metrics[mi]
    raw <- data[[m]]
    scored <- scale_directional(raw, m, method)
    # Clamp to [0, 1]
    scored <- pmax(0, pmin(1, scored))
    # Apply smoothing
    if (smooth_window > 1L && length(scored) >= smooth_window) {
      scored <- resilience_smooth_(scored, smooth_window, smooth_method)
    }
    data[[paste0(m, "_score")]] <- scored
    data[[paste0(m, "_category")]] <- cut(
      scored,
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE
    )
    w <- ifelse_(m %in% names(weights), weights[[m]], 0.01)
    weighted_mat[, mi] <- ifelse(is.na(scored), 0, scored * w)
    total_weight <- total_weight + w
  }
  data$composite_score <- rowSums(weighted_mat, na.rm = TRUE) / total_weight
  data$composite_category <- cut(
    data$composite_score,
    breaks = breaks,
    labels = labels,
    include.lowest = TRUE
  )
  structure(
    data,
    class = c("resilience", "tbl_df", "tbl", "data.frame")
  )
}

is_lower_better <- function(metric) {
  metric %in% c("vsi", "arch_lm", "cv", "recovery_time", "dfa_alpha")
}

default_resilience_weights <- function() {
  c(
    vsi = 0.256,
    cv = 0.226,
    arch_lm = 0.188,
    ac_ratio = 0.070,
    recovery_slope = 0.058,
    dfa_alpha = 0.036,
    sample_entropy = 0.028,
    recovery_time = 0.016
  )
}

resilience_empirical_means <- function() {
  data.frame(
    true_phase = c("Recovery", "Shock", "Stable", "Turbulent", "Volatile"),
    vsi = c(24.65, 75.61, 1455.42, 2857.7, 15764.89),
    arch_lm = c(18.05, 17.97, 9.47, 16.9, 2.03),
    cv = c(0.13, 0.16, 0.07, 0.16, 0.16),
    recovery_time = c(6.57, 6.35, 6.19, 6.17, 3.92),
    recovery_slope = c(0.61, -0.83, -0.58, 0.51, 0.28),
    sample_entropy = c(1.75, 0.93, 1.53, 1.57, 1.6),
    dfa_alpha = c(0.97, 0.86, 1.3, 1.45, 1.04),
    ac_ratio = c(1.39, 1.33, 2.8, 2.29, 1.79),
    stringsAsFactors = FALSE
  )
}

scale_directional <- function(values, metric, method) {
  if (metric == "sample_entropy") {
    inf_mask <- is.infinite(values)
    if (any(inf_mask)) {
      finite_vals <- values[!inf_mask & !is.na(values)]
      if (length(finite_vals) > 0L) {
        values[inf_mask] <- max(finite_vals) * 1.1
      }
    }
  }
  valid <- values[!is.na(values) & is.finite(values)]
  if (length(valid) < 2L) {
    return(rep(0.5, length(values)))
  }
  lower <- is_lower_better(metric)
  scaled <- switch(
    method,
    empirical_state_aware = {
      empirical <- resilience_empirical_means()
      good_states <- c("Stable", "Recovery")
      bad_states <- c("Shock")
      good_vals <- empirical[[metric]][empirical$true_phase %in% good_states]
      bad_vals  <- empirical[[metric]][empirical$true_phase %in% bad_states]
      good_vals <- good_vals[!is.na(good_vals) & is.finite(good_vals)]
      bad_vals  <- bad_vals[!is.na(bad_vals) & is.finite(bad_vals)]
      if (length(good_vals) < 1L || length(bad_vals) < 1L) {
        # Fall back to percentile
        sc <- scale_directional(values, metric, "empirical_percentile")
        return(sc)
      }
      if (lower) {
        excellent_th <- stats::quantile(good_vals, 0.75, na.rm = TRUE)
        troubled_th  <- stats::quantile(bad_vals, 0.25, na.rm = TRUE)
      } else {
        excellent_th <- stats::quantile(good_vals, 0.25, na.rm = TRUE)
        troubled_th  <- stats::quantile(bad_vals, 0.75, na.rm = TRUE)
      }
      sc <- numeric(length(values))
      for (i in seq_along(values)) {
        v <- values[i]
        if (is.na(v) || !is.finite(v)) { sc[i] <- NA_real_; next }
        if (lower) {
          if (v <= excellent_th) {
            sc[i] <- 0
          } else if (v >= troubled_th) {
            sc[i] <- 1
          } else {
            rng <- troubled_th - excellent_th
            sc[i] <- if (rng > 0) (v - excellent_th) / rng else 0.5
          }
        } else {
          if (v >= excellent_th) {
            sc[i] <- 0
          } else if (v <= troubled_th) {
            sc[i] <- 1
          } else {
            rng <- excellent_th - troubled_th
            sc[i] <- if (rng > 0) 1 - (v - troubled_th) / rng else 0.5
          }
        }
      }
      pmax(0, pmin(1, sc))
    },
    empirical_percentile = {
      pctiles <- stats::quantile(
        valid, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE
      )
      sc <- numeric(length(values))
      for (i in seq_along(values)) {
        v <- values[i]
        if (is.na(v) || !is.finite(v)) { sc[i] <- NA_real_; next }
        if (lower) {
          if      (v <= pctiles[1L]) sc[i] <- 0.0
          else if (v <= pctiles[2L]) sc[i] <- 0.2
          else if (v <= pctiles[3L]) sc[i] <- 0.4
          else if (v <= pctiles[4L]) sc[i] <- 0.6
          else if (v <= pctiles[5L]) sc[i] <- 0.8
          else                       sc[i] <- 1.0
        } else {
          if      (v >= pctiles[5L]) sc[i] <- 0.0
          else if (v >= pctiles[4L]) sc[i] <- 0.2
          else if (v >= pctiles[3L]) sc[i] <- 0.4
          else if (v >= pctiles[2L]) sc[i] <- 0.6
          else if (v >= pctiles[1L]) sc[i] <- 0.8
          else                       sc[i] <- 1.0
        }
      }
      sc
    },
    percentile = {
      sc <- stats::ecdf(valid)(values)
      if (!lower) sc <- 1 - sc
      sc
    },
    minmax = {
      mn <- min(valid); mx <- max(valid); rng <- mx - mn
      sc <- if (rng == 0) rep(0.5, length(values)) else (values - mn) / rng
      if (lower) sc else 1 - sc
    },
    zscore = {
      mu <- mean(valid); s <- stats::sd(valid)
      sc <- if (s == 0) rep(0, length(values)) else (values - mu) / s
      if (lower) sc else -sc
    },
    divide_max = {
      mx <- max(abs(valid))
      sc <- if (mx == 0) rep(0, length(values)) else values / mx
      if (lower) sc else 1 - sc
    }
  )
  scaled
}

resilience_smooth_ <- function(values, window_size, method) {
  n <- length(values)
  if (window_size <= 1L || n < window_size) return(values)
  smoothed <- numeric(n)
  half <- floor(window_size / 2)
  smooth_fn <- if (method == "median") stats::median else mean
  for (i in seq_len(n)) {
    s <- max(1L, i - half)
    e <- min(n, i + half)
    smoothed[i] <- smooth_fn(values[s:e], na.rm = TRUE)
  }
  smoothed
}

