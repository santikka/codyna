#' Rolling Spectral Early Warning Signals
#'
#' Computes rolling-window spectral metrics for a univariate time series,
#' estimating the spectral exponent (beta) and spectral ratio to detect
#' shifts in noise color that precede critical transitions. Reddening of
#' the power spectrum (increasing beta) is a hallmark of critical slowing
#' down.
#'
#' @details
#' The analysis proceeds in four stages:
#'
#' **1. Window extraction.** A rolling window of size `window` slides
#' across the series with the specified `align`ment. Each window
#' position produces one set of spectral metrics.
#'
#' **2. Detrending.** Before computing the spectrum, each window segment
#' is detrended according to `detrend`:
#'
#' * `"none"`: No detrending. Suitable when the data are already
#'   stationary within each window.
#' * `"linear"` (default): Removes a least-squares linear trend
#'   from each window. Appropriate for slowly drifting series.
#' * `"diff"`: First-differences the window segment. Effective
#'   for integrated (unit-root) processes but reduces the segment
#'   length by one observation.
#'
#' **3. Spectral estimation.** The power spectral density is estimated
#' by one of two methods:
#'
#' * `"periodogram"` (default): Computes the raw periodogram via
#'   [stats::spectrum()] with `method = "pgram"`. Fast and
#'   nonparametric. A 10 percent cosine taper is applied by default.
#' * `"ar"`: Fits an autoregressive model and derives the spectrum
#'   via [stats::spectrum()] with `method = "ar"`. Smoother estimates,
#'   but assumes the AR model is appropriate.
#'
#' **4. Metrics and classification.** For each window:
#'
#' * The **spectral exponent** (beta) is the negative slope of the
#'   log(power) vs. log(frequency) regression. Positive beta indicates
#'   reddened spectra (more low-frequency power). Rising beta over time
#'   signals critical slowing down.
#' * The **spectral ratio** is the ratio of total power at
#'   frequencies below the median to total power above the median.
#'   Values > 1 indicate a red-shifted spectrum.
#' * The **R-squared** of the log-log regression quantifies how
#'   well the power-law model fits the spectrum.
#'
#' **State classification** based on the spectral exponent:
#'
#' * beta < 0.5: `white_noise` -- uncorrelated fluctuations
#'   with flat spectrum.
#' * 0.5 <= beta < 1.5: `pink_noise` -- 1/f-like
#'   dynamics, moderate long-range correlations.
#' * 1.5 <= beta < 2.5: `red_noise` -- strong low-frequency
#'   dominance, consistent with critical slowing down.
#' * beta >= 2.5: `brownian` -- Brownian-motion-like
#'   dynamics with very strong low-frequency power.
#'
#' **Edge padding.** Window positions at the series boundaries that lack
#' sufficient data are filled by a three-step strategy: (1) forward-fill
#' from the first valid estimate to the series start, (2) backward-fill
#' from the last valid estimate to the series end, and (3) linear
#' interpolation for any internal gaps.
#'
#' @export
#' @param data \[`ts`, `numeric()`\]\cr
#'   Univariate time series data. Accepts a numeric vector or a `ts` object.
#' @param window \[`integer(1)`: `50L`\]\cr
#'   Rolling window size in observations.
#' @param align \[`character(1)`: `"right"`\]\cr
#'   Window alignment. The available options are:
#'   * `"right"`: each output value uses the `window` observations
#'     ending at that time point.
#'   * `"center"`: each output value is centered on the time point.
#'   * `"left"`: each output value uses the `window` observations
#'     starting at that time point.
#' @param method \[`character(1)`: `"periodogram"`\]\cr
#'   Spectral estimation method. The available options are:
#'   * `"periodogram"`: raw periodogram via [stats::spectrum()].
#'   * `"ar"`: autoregressive spectral estimate via [stats::spectrum()].
#' @param detrend \[`character(1)`: `"linear"`\]\cr
#'   Detrending applied to each window segment before spectral estimation.
#'   The available options are:
#'   * `"none"`: no detrending.
#'   * `"linear"`: remove a least-squares linear trend.
#'   * `"diff"`: first-difference the segment.
#' @param states \[`logical(1)`: `TRUE`\]\cr
#'   If `TRUE`, classify each time point into a noise-color state based
#'   on the spectral exponent. If `FALSE`, return the tibble without the
#'   `state` column.
#' @param min_points \[`integer(1)`: `20L`\]\cr
#'   Minimum number of non-NA observations required in a window for
#'   spectral estimation.
#' @return An object of class `"spectral"` (inheriting from
#'   [tibble::tibble()]) with the following columns:
#'
#'   * `time`: integer or Date time index.
#'   * `value`: original series values.
#'   * `spectral_exponent`: rolling spectral exponent (beta).
#'   * `spectral_ratio`: rolling low-to-high frequency power ratio.
#'   * `r_squared`: R-squared of the log-log fit at each window.
#'   * `state`: (when `states = TRUE`) factor with levels
#'       `"white_noise"`, `"pink_noise"`, `"red_noise"`, `"brownian"`.
#'
#'   Attributes: `window`, `align`, `method`, `detrend`.
#'
#' @references
#' Held, H. & Kleinen, T. (2004). Detection of climate system bifurcations
#' by degenerate fingerprinting. \emph{Geophysical Research Letters}, 31,
#' L23207. \doi{10.1029/2004GL020972}
#'
#' Kleinen, T., Held, H., & Petschel-Held, G. (2003). The potential role
#' of spectral properties in detecting thresholds in the Earth system:
#' application to the thermohaline circulation. \emph{Ocean Dynamics}, 53,
#' 53--63. \doi{10.1007/s10236-002-0023-6}
#'
#' Dakos, V., Scheffer, M., van Nes, E.H., Brovkin, V., Petoukhov, V.,
#' & Held, H. (2008). Slowing down as an early warning signal for abrupt
#' climate change. \emph{Proceedings of the National Academy of Sciences},
#' 105(38), 14308--14312. \doi{10.1073/pnas.0802430105}
#'
#' @seealso [hurst()] for long-range dependence analysis;
#'   [compute_trend()] for trend classification;
#'   [detect_multivariate_warnings()] for multivariate EWS.
#' @examples
#' \donttest{
#' set.seed(42)
#' x <- cumsum(rnorm(500))
#' sp <- spectral_ews(x, window = 50L)
#' plot(sp)
#' summary(sp)
#' }
spectral_ews <- function(data, window = 50L, align = "right",
                         method = "periodogram", detrend = "linear",
                         states = TRUE, min_points = 20L) {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)
  method <- check_match(method, c("periodogram", "ar"))
  detrend <- check_match(detrend, c("none", "linear", "diff"))
  align <- check_match(align, c("right", "center", "left"))
  check_range(window, type = "integer", min = 2L, max = n)
  check_range(min_points, type = "integer", min = 4L, max = window)
  check_flag(states)
  beta_vec <- rep(NA_real_, n)
  ratio_vec <- rep(NA_real_, n)
  r2_vec <- rep(NA_real_, n)
  half_left <- (window - 1L) %/% 2L
  half_right <- window - 1L - half_left
  iter_start <- 1L
  iter_end <- n
  if (align == "center") {
    iter_start <- 1L + half_left
    iter_end <- n - half_right
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
    seg <- values[w_idx]
    n_valid <- sum(!is.na(seg))
    if (n_valid < min_points) next
    # Remove NAs for spectral computation
    seg_clean <- seg[!is.na(seg)]
    # Detrend
    seg_detrended <- spectral_detrend_(seg_clean, detrend)
    if (length(seg_detrended) < 10L) next
    # Compute spectral metrics
    result <- spectral_metrics_(seg_detrended, method)
    if (!is.null(result)) {
      beta_vec[k] <- result$beta
      ratio_vec[k] <- result$ratio
      r2_vec[k] <- result$r_squared
    }
  }
  first_valid <- which(!is.na(beta_vec))[1L]
  last_valid <- rev(which(!is.na(beta_vec)))[1L]
  if (!is.na(first_valid) && !is.na(last_valid)) {
    # Forward-fill beginning
    if (first_valid > 1L) {
      beta_vec[seq_len(first_valid - 1L)] <- beta_vec[first_valid]
      ratio_vec[seq_len(first_valid - 1L)] <- ratio_vec[first_valid]
      r2_vec[seq_len(first_valid - 1L)] <- r2_vec[first_valid]
    }
    # Backward-fill end
    if (last_valid < n) {
      beta_vec[(last_valid + 1L):n] <- beta_vec[last_valid]
      ratio_vec[(last_valid + 1L):n] <- ratio_vec[last_valid]
      r2_vec[(last_valid + 1L):n] <- r2_vec[last_valid]
    }
    # Linear interpolation for internal gaps
    na_pos <- which(is.na(beta_vec))
    internal_na <- na_pos[na_pos > first_valid & na_pos < last_valid]
    if (length(internal_na) > 0L) {
      valid_pos <- which(!is.na(beta_vec))
      beta_vec <- stats::approx(
        valid_pos, beta_vec[valid_pos], xout = seq_len(n), rule = 2
      )$y
      ratio_vec <- stats::approx(
        valid_pos, ratio_vec[valid_pos], xout = seq_len(n), rule = 2
      )$y
      r2_vec <- stats::approx(
        valid_pos, r2_vec[valid_pos], xout = seq_len(n), rule = 2
      )$y
    }
  }
  out <- tibble::tibble(
    time = time,
    value = data$values,
    spectral_exponent = beta_vec,
    spectral_ratio = ratio_vec,
    r_squared = r2_vec
  )
  if (states) {
    state_labels <- spectral_classify_(beta_vec)
    out$state <- factor(
      state_labels,
      levels = c("white_noise", "pink_noise", "red_noise", "brownian")
    )
  }
  structure(
    out,
    window = window,
    align = align,
    method = method,
    detrend = detrend,
    class = c("spectral", "tbl_df", "tbl", "data.frame")
  )
}

#' Detrend a numeric segment
#'
#' Applies detrending to a numeric vector prior to spectral estimation.
#' Supports no detrending, linear detrending (OLS residuals), and
#' first-differencing.
#'
#' @param x Numeric vector.
#' @param method One of `"none"`, `"linear"`, or `"diff"`.
#' @return Detrended numeric vector.
#' @noRd
spectral_detrend_ <- function(x, method) {
  switch(
    method,
    none = x,
    linear = {
      n <- length(x)
      t_idx <- seq_len(n)
      t_mean <- mean(t_idx)
      x_mean <- mean(x, na.rm = TRUE)
      xx_var <- sum((t_idx - t_mean)^2)
      if (xx_var > 0) {
        xy_cov <- sum((t_idx - t_mean) * (x - x_mean))
        slope <- xy_cov / xx_var
        intercept <- x_mean - slope * t_mean
        x - (intercept + slope * t_idx)
      } else {
        x - x_mean
      }
    },
    diff = diff(x)
  )
}

#' Compute spectral metrics for a single window
#'
#' Fits a log-log regression to the power spectrum and computes the spectral
#' exponent (beta), spectral ratio (low/high frequency power), and R-squared.
#'
#' @param x Numeric vector (window data, already detrended).
#' @param spec_method One of `"periodogram"` or `"ar"`.
#' @return A named list with elements `beta`, `ratio`, `r_squared`, or `NULL`
#'   on failure.
#' @noRd
spectral_metrics_ <- function(x, spec_method) {
  if (length(x) < 10L || stats::var(x, na.rm = TRUE) < 1e-12) {
    return(NULL)
  }
  spec_result <- if (spec_method == "periodogram") {
    try_(
      stats::spectrum(
        x, method = "pgram", plot = FALSE, detrend = FALSE, taper = 0.1
      )
    )
  } else {
    try_(stats::spectrum(x, method = "ar", plot = FALSE))
  }
  if (inherits(spec_result, "try-error") || is.null(spec_result)) {
    return(NULL)
  }
  freq <- spec_result$freq
  power <- spec_result$spec
  # Remove zero or negative frequencies/power for log transform
  valid <- freq > 0 & power > 0
  if (sum(valid) < 4L) {
    return(NULL)
  }
  freq <- freq[valid]
  power <- power[valid]
  log_freq <- log(freq)
  log_power <- log(power)
  # OLS: log(power) ~ log(freq)
  x_m <- mean(log_freq)
  y_m <- mean(log_power)
  xx_var <- sum((log_freq - x_m)^2)
  if (xx_var < 1e-12) return(NULL)
  xy_cov <- sum((log_freq - x_m) * (log_power - y_m))
  slope <- xy_cov / xx_var
  intercept <- y_m - slope * x_m
  y_pred <- intercept + slope * log_freq
  ss_tot <- sum((log_power - y_m)^2)
  ss_res <- sum((log_power - y_pred)^2)
  r_sq <- ifelse_(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
  # Spectral exponent: negate slope so positive beta = reddening
  beta <- -slope
  # Spectral ratio: sum(power at low freq) / sum(power at high freq)
  median_freq <- stats::median(freq)
  low_power <- sum(power[freq <= median_freq])
  high_power <- sum(power[freq > median_freq])
  ratio <- ifelse_(high_power > 0, low_power / high_power, NA_real_)
  list(beta = beta, ratio = ratio, r_squared = r_sq)
}

