#' Detect Multivariate Early Warning Signals
#'
#' Performs multivariate early warning signal (EWS) analysis on a system
#' represented by multiple time series. The function computes a suite of
#' metrics --- based on Min/Max Autocorrelation Factors (MAF), Principal
#' Component Analysis (PCA), and covariance structure --- over rolling or
#' expanding windows to detect rising instability before critical transitions.
#'
#' @details
#' Multivariate EWS extend univariate indicators to systems where several
#' variables interact. Individual-variable summaries (mean/max of standard
#' deviations and AR(1) coefficients) are complemented by dimension-reduction
#' indicators that capture coordinated changes across variables.
#'
#' The eleven available metrics are:
#'
#' * `meanSD`: Mean standard deviation across all variables in the
#'   window. Rising values indicate broadening fluctuations system-wide.
#' * `maxSD`: Maximum standard deviation among variables. Flags when any
#'   single variable becomes highly volatile.
#' * `meanAR`: Mean lag-1 autoregressive coefficient across variables.
#'   Increasing values signal system-wide critical slowing down.
#' * `maxAR`: Maximum AR(1) coefficient among variables. Detects
#'   critical slowing down concentrated in a single variable.
#' * `eigenMAF`: Ratio of the minimum eigenvalue of the MAF
#'   decomposition to the total. Declining ratios (increasingly negative
#'   after sign inversion) indicate rising spatial correlation
#'   (Dakos et al., 2012).
#' * `mafAR`: AR(1) coefficient of the first MAF component. The MAF
#'   transformation orders components by temporal smoothness, so the first
#'   component captures the slowest mode. Rising mafAR indicates critical
#'   slowing down in the most persistent system dimension.
#' * `mafSD`: Standard deviation of the first MAF component. Sign-inverted
#'   so that rising values indicate increasing fluctuation in the slowest
#'   mode.
#' * `pcaAR`: AR(1) coefficient of the first principal component. PCA
#'   captures the dominant variance direction; rising pcaAR signals that
#'   this direction is becoming more persistent.
#' * `pcaSD`: Standard deviation of the first principal component.
#'   Rising values indicate that the dominant variance mode is
#'   amplifying.
#' * `eigenCOV`: Dominant eigenvalue of the covariance matrix of
#'   standardized data. Increasing values indicate growing variance
#'   concentration along one axis.
#' * `maxCOV`: Maximum off-diagonal element of the covariance matrix.
#'   Rising values signal strengthening cross-variable coupling.
#'
#' **Rolling window** (`method = "rolling"`): a fixed-width window slides
#' across the series. Each metric is standardized (z-scored) across all
#' window positions for comparability. Kendall's tau trend statistic
#' is computed for each metric against time; metrics with
#' \eqn{\tau > 0.7} are flagged as showing significant upward trends.
#'
#' **Expanding window** (`method = "expanding"`): the window starts at
#' `min_window` and grows one observation at a time. A two-stage
#' standardization is applied (matching the EWSmethods implementation):
#' (1) each metric's raw value is z-scored against its own expanding
#' history up to that point; (2) the resulting z-scores are z-scored
#' again using an expanding cumulative mean and SD. This double
#' standardization puts all metrics on a common dimensionless scale
#' regardless of their original units. A warning is flagged when the
#' doubly-standardized strength exceeds `threshold` after the
#' `burn_in` period. Classification into system states (Stable,
#' Vulnerable, Warning, Critical, Failing) is based on
#' how many metrics simultaneously cross the threshold.
#'
#' **Dimension reduction series**: regardless of method, the function
#' also extracts the first MAF (MAF1) and first principal component
#' (PC1) time series, standardizes them, and flags points where the
#' absolute scaled value exceeds the threshold. These provide a
#' low-dimensional summary of the system's trajectory.
#'
#' @export
#' @param data \[`data.frame`\]\cr
#'   A data frame with a time column and at least two numeric time series
#'   columns. If `time_col` is not specified, the first column is assumed
#'   to be time.
#' @param time_col \[`character(1)` or `NULL`: `NULL`\]\cr
#'   Name of the time column. If `NULL`, the first column is used.
#' @param metrics \[`character()`: `"all"`\]\cr
#'   Metrics to compute. Pass `"all"` for all eleven, or a subset of:
#'   `"meanSD"`, `"maxSD"`, `"meanAR"`, `"maxAR"`, `"eigenMAF"`,
#'   `"mafAR"`, `"mafSD"`, `"pcaAR"`, `"pcaSD"`, `"eigenCOV"`,
#'   `"maxCOV"`.
#' @param method \[`character(1)`: `"rolling"`\]\cr
#'   Analysis approach: `"rolling"` or `"expanding"`.
#' @param window \[`numeric(1)`: `50`\]\cr
#'   For `"rolling"`: window size as a percentage of total series length.
#'   For `"expanding"`: not used directly (the window grows from
#'   `min_window` to `n`).
#' @param burn_in \[`integer(1)`: `10L`\]\cr
#'   Number of initial observations after which warnings begin to be
#'   flagged. Only used when `method = "expanding"`.
#' @param threshold \[`numeric(1)`: `2`\]\cr
#'   Number of standard deviations for warning detection. Only used when
#'   `method = "expanding"`.
#' @param tail_direction \[`character(1)`: `"one.tailed"`\]\cr
#'   Whether to test for warnings using `"one.tailed"` (positive
#'   deviations only) or `"two.tailed"` (both directions). Only used when
#'   `method = "expanding"`.
#'
#' @return An object of class `"multi_ews"`, a tibble (inheriting from
#'   `tbl_df`, `tbl`, `data.frame`) with metadata stored as attributes.
#'
#'   **Rolling** columns: `time`, `metric`, `score`, `std`.
#'   Attributes: `orig_data`, `cor` (Kendall tau per metric),
#'   `method` (`"rolling"`), `dimension_reduction`.
#'
#'   **Expanding** columns: `time`, `metric`, `score`, `z_score`,
#'   `detected` (integer 0/1).
#'   Attributes: `orig_data`, `threshold`, `tail_direction`,
#'   `classification` (tibble with `time`, `count`, `state`),
#'   `method` (`"expanding"`), `dimension_reduction`.
#'
#' @references
#' Dakos, V., Carpenter, S. R., Brock, W. A., Ellison, A. M., Guttal, V.,
#' Ives, A. R., ... & Scheffer, M. (2012). Methods for detecting early
#' warnings of critical transitions in time series illustrated using
#' simulated ecological data. *PLoS ONE*, 7(7), e41010.
#' \doi{10.1371/journal.pone.0041010}
#'
#' Scheffer, M., Bascompte, J., Brock, W. A., Brovkin, V., Carpenter,
#' S. R., Dakos, V., ... & Sugihara, G. (2009). Early-warning signals
#' for critical transitions. *Nature*, 461(7260), 53--59.
#' \doi{10.1038/nature08227}
#'
#' @seealso [generate_tipping_data()] for creating test data with known
#'   tipping points; [plot.multi_ews()] for visualization.
#' @examples
#' \donttest{
#' tip <- generate_tipping_data(n_time = 100, n_vars = 3, tipping_point = 60)
#' ews <- detect_multivariate_warnings(tip, method = "expanding", window = 50)
#' plot(ews)
#'
#' ews_roll <- detect_multivariate_warnings(tip, method = "rolling", window = 50)
#' plot(ews_roll)
#' }
detect_multivariate_warnings <- function(data, time_col = NULL, metrics = "all",
                                         method = "rolling", window = 50,
                                         burn_in = 10L, threshold = 2,
                                         tail_direction = "one.tailed") {
  check_missing(data)
  stopifnot_(
    is.data.frame(data) && ncol(data) >= 3L,
    "{.arg data} must be a data frame with a time column and at least
      two time series columns."
  )
  if (is.null(time_col)) {
    message_("No {.arg time_col} specified. Using the first column as time.")
    idx <- 1L
  } else {
    idx <- which(names(data) == time_col)
    stopifnot_(
      length(idx) == 1L,
      "Column {.val {time_col}} not found in {.arg data}."
    )
  }
  original_data <- data
  other_idx <- setdiff(seq_len(ncol(data)), idx)
  data_internal <- data[, c(idx, other_idx)]
  stopifnot_(
    all(vapply(data_internal[, -1L], is.numeric, logical(1L))),
    "All time series columns must be {.cls numeric}."
  )
  stopifnot_(
    all(!is.na(data_internal[, -1L])),
    "{.arg data} contains missing values.
      Interpolation before analysis is recommended."
  )
  method <- check_match(method, c("rolling", "expanding"))
  tail_direction <- check_match(tail_direction, c("one.tailed", "two.tailed"))
  available_metrics <- c(
    "meanSD", "maxSD", "meanAR", "maxAR", "eigenMAF", "mafAR",
    "mafSD", "pcaAR", "pcaSD", "eigenCOV", "maxCOV"
  )
  invalid <- setdiff(metrics, c(available_metrics, "all"))
  stopifnot_(
    length(invalid) == 0,
    c(
      "Invalid metrics: {.val {invalid}}.",
      `i` = "Available: {.val {available_metrics}}."
    )
  )
  if ("all" %in% metrics) {
    metrics <- available_metrics
  }
  ews_results <- mews_calculate(
    data = data_internal,
    metrics = metrics,
    window = window,
    method  = method,
    burn_in = burn_in
  )
  if (method == "rolling") {
    # Kendall tau per metric
    correlations <- vapply(
      ews_results$raw[, metrics, drop = FALSE],
      function(x) {
        if (sum(!is.na(x)) < 4L) {
          return(NA_real_)
        }
        suppressWarnings(
          stats::cor.test(x, seq_along(x), method = "kendall")$estimate
        )
      },
      numeric(1L)
    )
    long_data <- tidyr::pivot_longer(
      ews_results$raw,
      cols = -!!rlang::sym("time"),
      names_to = "metric",
      values_to = "score"
    )
    # Standardize scores per metric
    long_data$std <- NA_real_
    for (m in unique(long_data$metric)) {
      idx <- long_data$metric == m
      vals <- long_data$score[idx]
      if (sum(!is.na(vals)) > 1L &&
          stats::sd(vals, na.rm = TRUE) > 1e-10) {
        long_data$std[idx] <- as.numeric(scale(vals))
      } else {
        long_data$std[idx] <- 0
      }
    }

  } else {
    # Expanding: metric.score values are already z-scored in Stage 1
    # (inside mews_expanding_metrics). Here we apply Stage 2: a second
    # expanding z-score on those z-scores, matching EWSmethods::wMAF().
    long_data <- tidyr::pivot_longer(
      ews_results$raw,
      cols = -!!rlang::sym("time"),
      names_to = "metric",
      values_to = "score"
    )
    long_data <- long_data[order(long_data$metric, long_data$time), ]
    # Stage 2: expanding cumulative mean/sd of the Stage-1 z-scores
    long_data <- long_data |>
      dplyr::group_by(!!rlang::sym("metric")) |>
      dplyr::arrange(!!rlang::sym("time"), .by_group = TRUE) |>
      dplyr::mutate(
        rolling.mean = mews_rolling_mean(!!rlang::sym("score")),
        rolling.sd = mews_rolling_sd(!!rlang::sym("score"))
      ) |>
      dplyr::ungroup()
    long_data$z_score <- (long_data$score - long_data$rolling.mean) /
      long_data$rolling.sd
    burn_in_time <- data_internal[[1L]][burn_in]
    if (tail_direction == "two.tailed") {
      crossed <- long_data$z_score > threshold | long_data$z_score < -threshold
    } else {
      crossed <- long_data$z_score > threshold
    }
    crossed[is.na(crossed)] <- FALSE
    long_data$detected <- ifelse(
      long_data$time < burn_in_time, 0L,
      ifelse(crossed, 1L, 0L)
    )
    # Remove intermediate columns
    long_data$rolling.mean <- NULL
    long_data$rolling.sd <- NULL
  }
  out_tbl <- tibble::as_tibble(long_data)
  classification <- NULL
  if (method == "expanding") {
    n_metrics_count <- length(unique(out_tbl$metric))
    classification <- mews_classify(out_tbl, n_metrics_count)
  }
  dr <- mews_dimension_reduction(
    data = data_internal,
    window = onlyif(method == "rolling", window),
    method = method,
    threshold = ifelse_(method == "expanding", threshold, 2)
  )
  if (method == "rolling") {
    structure(
      out_tbl,
      orig_data = original_data,
      cor = correlations,
      method  = "rolling",
      dimension_reduction = dr,
      class = c("multi_ews", "tbl_df", "tbl", "data.frame")
    )
  } else {
    structure(
      out_tbl,
      orig_data = original_data,
      threshold = threshold,
      tail_direction = tail_direction,
      classification = classification,
      method = "expanding",
      dimension_reduction = dr,
      class = c("multi_ews", "tbl_df", "tbl", "data.frame")
    )
  }
}

mews_calculate <- function(data, metrics, window, method, burn_in = 10L) {
  time_vec <- data[, 1L]
  ts_data  <- as.matrix(data[, -1L])
  if (method == "rolling") {
    win_size_pts <- round(nrow(ts_data) * window / 100)
    min_window <- max(3L, ncol(ts_data) + 1L)
    if (win_size_pts < min_window) {
      warning_(
        "Window size too small ({win_size_pts} points).
        Using minimum of {min_window} points."
      )
      win_size_pts <- min_window
    }
    stopifnot_(
      nrow(data) >= win_size_pts,
      "Window size ({win_size_pts}) must not exceed data length ({nrow(ts_data)})."
    )
    n_obs <- nrow(ts_data)
    rolling_list <- lapply(
      seq(win_size_pts, n_obs),
      function(i) {
        w <- ts_data[(i - win_size_pts + 1L):i, , drop = FALSE]
        mews_window_metrics(w, metrics)
      }
    )
    rolling_metrics <- do.call(base::rbind, rolling_list)
    results_df <- as.data.frame(rolling_metrics)
    results_df$time <- time_vec[win_size_pts:length(time_vec)]
    list(raw = results_df)
  } else {
    # Expanding window: loop starts from burn_in (matching EWSmethods).
    # min_window is the structural minimum (need more obs than variables).
    min_window <- max(3L, ncol(ts_data) + 1L)
    start_from <- max(min_window, burn_in)
    stopifnot_(
      start_from < nrow(ts_data),
      "Not enough data points ({nrow(ts_data)}) for expanding analysis."
    )
    # Accumulators for incremental z-scoring
    env <- new.env(parent = emptyenv())
    env$roll_eigen <- list()
    env$roll_ar <- list()
    env$roll_sd <- list()
    env$roll_pca_ar <- list()
    env$roll_pca_sd <- list()
    env$roll_cov <- list()
    env$roll_mean_ar <- list()
    env$roll_max_ar <- list()
    env$roll_mean_sd <- list()
    env$roll_max_sd <- list()
    env$roll_cov_eigen <- list()
    results_list <- lapply(
      start_from:nrow(ts_data),
      function(i) {
      window_data <- ts_data[1L:i, , drop = FALSE]
      mews_expanding_metrics(window_data, metrics, i, env)
      }
    )
    results_df <- as.data.frame(do.call(base::rbind, results_list))
    results_df$time <- time_vec[start_from:nrow(ts_data)]
    list(raw = results_df)
  }
}

mews_window_metrics <- function(window_data, metrics) {
  if (nrow(window_data) < 3L) {
    return(rep(NA_real_, length(metrics)))
  }
  out <- rep(NA_real_, length(metrics))
  names(out) <- metrics
  for (metric in metrics) {
    if (metric %in% c("meanSD", "maxSD")) {
      sds <- apply(window_data, 2L, stats::sd, na.rm = TRUE)
      if (metric == "meanSD") {
        out[metric] <- mean(sds, na.rm = TRUE)
      }
      if (metric == "maxSD") {
        out[metric] <- max(sds, na.rm = TRUE)
      }
    }
    if (metric %in% c("meanAR", "maxAR")) {
      ars <- apply(window_data, 2L, mews_ar_robust)
      if (metric == "meanAR") {
        out[metric] <- mean(ars, na.rm = TRUE)
      }
      if (metric == "maxAR") {
        out[metric] <- max(ars, na.rm = TRUE)
      }
    }
    if (metric %in% c("pcaAR", "pcaSD", "eigenCOV", "maxCOV")) {
      window_scaled <- scale(window_data)
      V_scaled <- stats::cov(window_scaled)
      if (all(!is.na(V_scaled))) {
        if (metric %in% c("pcaAR", "pcaSD")) {
          pca_res <- stats::prcomp(
            window_scaled, scale. = FALSE, center = FALSE
          )
          pca_ts <- pca_res$x[, 1L]
          if (metric == "pcaAR") out[metric] <- mews_ar_robust(pca_ts)
          if (metric == "pcaSD") out[metric] <- stats::sd(pca_ts, na.rm = TRUE)
        }
        if (metric == "eigenCOV") {
          out[metric] <- eigen(V_scaled)$values[1L]
        }
        if (metric == "maxCOV") {
          out[metric] <- max(V_scaled[lower.tri(V_scaled, diag = FALSE)])
        }
      }
    }
    if (metric %in% c("eigenMAF", "mafAR", "mafSD")) {
      maf_res <- mews_maf(window_data)
      if (!is.null(maf_res)) {
        if (metric == "eigenMAF") {
          # Negate so that increasing values indicate rising spatial
          # correlation (approaching transition), matching EWSmethods
          out[metric] <- -min(maf_res$eigen_values) /
            sum(maf_res$eigen_values)
        }
        if (metric %in% c("mafAR", "mafSD")) {
          maf_ts_1 <- maf_res$mafs[, 1L]
          if (metric == "mafAR") {
            out[metric] <- mews_ar_robust(maf_ts_1)
          }
          if (metric == "mafSD") {
            out[metric] <- stats::sd(maf_ts_1, na.rm = TRUE)
          }
        }
      }
    }
  }
  out
}

mews_expanding_metrics <- function(window_data, metrics, i, env) {
  # Stage 1 z-score: each metric's raw value is z-scored against its own

  # expanding history, matching EWSmethods::wMAF() exactly.
  # All metrics use the same formula; only eigenMAF gets sign-inverted.
  result_vec <- rep(NA_real_, length(metrics))
  names(result_vec) <- metrics
  zscore_ <- function(val, acc_list) {
    vals <- unlist(acc_list)
    if (length(vals) > 1L) {
      (val - mean(vals, na.rm = TRUE)) / stats::sd(vals, na.rm = TRUE)
    } else {
      NA_real_
    }
  }
  for (metric in metrics) {
    if (metric %in% c("meanSD", "maxSD")) {
      sds <- apply(window_data, 2L, stats::sd, na.rm = TRUE)
      if (metric == "meanSD") {
        env$roll_mean_sd[[i]] <- mean(sds, na.rm = TRUE)
        result_vec[metric] <- zscore_(
          env$roll_mean_sd[[i]],
          env$roll_mean_sd
        )
      }
      if (metric == "maxSD") {
        env$roll_max_sd[[i]] <- max(sds, na.rm = TRUE)
        result_vec[metric] <- zscore_(
          env$roll_max_sd[[i]],
          env$roll_max_sd
        )
      }
    }
    if (metric %in% c("meanAR", "maxAR")) {
      ars <- apply(window_data, 2L, mews_ar_robust)
      if (metric == "meanAR") {
        env$roll_mean_ar[[i]] <- mean(ars, na.rm = TRUE)
        result_vec[metric] <- zscore_(
          env$roll_mean_ar[[i]],
          env$roll_mean_ar
        )
      }
      if (metric == "maxAR") {
        env$roll_max_ar[[i]] <- max(ars, na.rm = TRUE)
        result_vec[metric] <- zscore_(
          env$roll_max_ar[[i]],
          env$roll_max_ar
        )
      }
    }
    if (metric %in% c("eigenMAF", "mafAR", "mafSD")) {
      maf_res <- mews_maf(window_data)
      if (!is.null(maf_res)) {
        if (metric == "eigenMAF") {
          env$roll_eigen[[i]] <- min(
            maf_res$eigen_values / sum(maf_res$eigen_values)
          )
          result_vec[metric] <- -1 * zscore_(
            env$roll_eigen[[i]],
            env$roll_eigen
          )
        }
        if (metric == "mafAR") {
          env$roll_ar[[i]] <- mews_ar_robust(maf_res$mafs[, 1L])
          result_vec[metric] <- zscore_(env$roll_ar[[i]], env$roll_ar)
        }
        if (metric == "mafSD") {
          env$roll_sd[[i]] <- stats::sd(maf_res$mafs[, 1L], na.rm = TRUE)
          result_vec[metric] <- zscore_(env$roll_sd[[i]], env$roll_sd)
        }
      }
    }
    if (metric %in% c("pcaAR", "pcaSD")) {
      pca_res <- stats::prcomp(
        scale(window_data),
        scale. = FALSE,
        center = FALSE
      )
      if (metric == "pcaAR") {
        env$roll_pca_ar[[i]] <- mews_ar_robust(pca_res$x[, 1L])
        result_vec[metric] <- zscore_(env$roll_pca_ar[[i]], env$roll_pca_ar)
      }
      if (metric == "pcaSD") {
        env$roll_pca_sd[[i]] <- stats::sd(pca_res$x[, 1L], na.rm = TRUE)
        result_vec[metric] <- zscore_(env$roll_pca_sd[[i]], env$roll_pca_sd)
      }
    }
    if (metric %in% c("eigenCOV", "maxCOV")) {
      V_scaled <- stats::cov(scale(window_data))
      if (metric == "eigenCOV") {
        env$roll_cov_eigen[[i]] <- max(eigen(V_scaled)$values)
        result_vec[metric] <- zscore_(
          env$roll_cov_eigen[[i]],
          env$roll_cov_eigen
        )
      }
      if (metric == "maxCOV") {
        env$roll_cov[[i]] <- max(V_scaled[lower.tri(V_scaled, diag = FALSE)])
        result_vec[metric] <- zscore_(env$roll_cov[[i]], env$roll_cov)
      }
    }
  }
  result_vec
}

mews_maf <- function(window_data) {
  if (nrow(window_data) < 3L || ncol(window_data) < 2L) {
    return(NULL)
  }
  result <- try_({
    p <- ncol(window_data)
    n <- nrow(window_data)
    x <- scale(window_data, center = TRUE, scale = TRUE)
    if (any(is.na(x)) || any(!is.finite(x))) {
      return(NULL)
    }
    V <- stats::cov(x)
    svd_V <- svd(V)
    if (any(svd_V$d < .Machine$double.eps)) {
      return(NULL)
    }
    a <- svd_V$u %*% diag(svd_V$d^(-0.5)) %*% t(svd_V$u)
    y <- x %*% a
    dy <- apply(y, 2L, diff)
    cov_dy <- stats::cov(dy) * (n - 2L)
    svd_cov_dy <- svd(cov_dy / (n - 1L))
    u_dy <- svd_cov_dy$u[, p:1L, drop = FALSE]
    aa <- a %*% u_dy
    aa <- apply(aa, 2L, function(col) col / sqrt(sum(col^2)))
    maf <- x %*% aa
    # Sign correction: ensure positive correlation with time
    neg_time_cor <- diag(
      stats::cov(maf, matrix(rep(seq_len(n), p), n, p))
    ) < 0
    mafs <- t(
      apply(
        maf,
        1L,
        function(row) (-neg_time_cor + !neg_time_cor) * row
      )
    )
    d <- svd_cov_dy$d[p:1L]
    list(mafs = mafs, eigen_values = d)
  })
  if (inherits(result, "try-error")) {
    return(NULL)
  }
  result
}

mews_ar_robust <- function(x) {
  if (length(x) < 3L || all(is.na(x)) ||
      stats::sd(x, na.rm = TRUE) < 1e-10) {
    return(NA_real_)
  }
  result <- try_({
    ar_result <- stats::ar.ols(
      x,
      aic = FALSE,
      order.max = 1L,
      dmean = FALSE,
      intercept = FALSE
    )
    ifelse_(
      length(ar_result$ar) > 0L,
      ar_result$ar[1L],
      NA_real_
    )
  })
  if (inherits(result, "try-error")) {
    return(NA_real_)
  }
  result
}

mews_rolling_mean <- function(x) {
  k <- length(x)
  result <- numeric(k)
  for (i in seq_len(k)) {
    result[i] <- mean(x[1L:i], na.rm = TRUE)
  }
  result
}

mews_rolling_sd <- function(x) {
  k <- length(x)
  result <- numeric(k)
  for (i in seq_len(k)) {
    result[i] <- stats::sd(x[1L:i], na.rm = TRUE)
  }
  result
}

mews_classify <- function(ews_data, n_metrics) {
  ews_valid <- ews_data[!is.na(ews_data$z_score), ]
  if (nrow(ews_valid) == 0L) {
    return(NULL)
  }
  summary_by_time <- stats::aggregate(
    detected ~ time,
    data = ews_valid,
    FUN  = sum
  )
  names(summary_by_time) <- c("time", "count")
  if (n_metrics == 1L) {
    summary_by_time$state <- cut(
      summary_by_time$count,
      breaks = c(-Inf, 0, 1, Inf),
      labels = c("Stable", "Warning", "Critical"),
      right  = TRUE
    )
  } else {
    summary_by_time$state <- cut(
      summary_by_time$count,
      breaks = c(-Inf, 0, 1, 2, 0.5 * n_metrics, Inf),
      labels = c("Stable", "Vulnerable", "Warning",
                 "Critical", "Failing"),
      right  = TRUE
    )
  }
  tibble::as_tibble(summary_by_time)
}

mews_dimension_reduction <- function(data, window, method, threshold) {
  time_vec <- data[, 1L]
  ts_data  <- as.matrix(data[, -1L])
  maf_pc_for_window <- function(window_data) {
    if (nrow(window_data) < 3L || ncol(window_data) < 2L) {
      return(NULL)
    }
    result <- try_({
      p <- ncol(window_data)
      n <- nrow(window_data)
      x <- scale(window_data, center = TRUE, scale = TRUE)
      if (any(is.na(x)) || any(!is.finite(x))) {
        return(NULL)
      }
      V <- stats::cov(x)
      svd_V <- svd(V)
      if (any(svd_V$d < .Machine$double.eps)) {
        return(NULL)
      }
      a <- svd_V$u %*% diag(svd_V$d^(-0.5)) %*% t(svd_V$u)
      y <- x %*% a
      dy <- apply(y, 2L, diff)
      cov_dy <- stats::cov(dy) * (n - 2L)
      svd_cov_dy <- svd(cov_dy / (n - 1L))
      u_dy <- svd_cov_dy$u[, p:1L, drop = FALSE]
      aa <- a %*% u_dy
      aa <- apply(aa, 2L, function(col) col / sqrt(sum(col^2)))
      maf <- x %*% aa
      neg_time_cor <- diag(
        stats::cov(maf, matrix(rep(seq_len(n), p), n, p))
      ) < 0
      mafs <- t(
        apply(
          maf,
          1L,
          function(row) (-neg_time_cor + !neg_time_cor) * row
        )
      )
      pca_res <- stats::prcomp(x, scale. = FALSE, center = FALSE)
      pc1_ts <- pca_res$x[, 1L]
      list(maf1 = mafs[, 1L], pc1 = pc1_ts)
    })
    if (inherits(result, "try-error")) {
      return(NULL)
    }
    result
  }
  if (method == "rolling") {
    win_size_pts <- round(nrow(ts_data) * window / 100)
    min_window <- max(3L, ncol(ts_data) + 1L)
    win_size_pts <- max(min_window, win_size_pts)
    n_obs_dr <- nrow(ts_data)
    maf_pc_list <- lapply(
      seq(win_size_pts, n_obs_dr),
      function(i) {
        w <- ts_data[(i - win_size_pts + 1L):i, , drop = FALSE]
        res <- maf_pc_for_window(w)
        if (!is.null(res)) {
          c(maf1 = res$maf1[nrow(w)], pc1 = res$pc1[nrow(w)])
        } else {
          c(maf1 = NA_real_, pc1 = NA_real_)
        }
      }
    )
    maf_pc_results <- do.call(base::rbind, maf_pc_list)
    results_df <- as.data.frame(maf_pc_results)
    results_df$time <- time_vec[win_size_pts:length(time_vec)]
  } else {
    min_window <- max(3L, ncol(ts_data) + 1L)
    results_list <- lapply(
      min_window:nrow(ts_data),
      function(i) {
        window_data <- ts_data[1L:i, , drop = FALSE]
        res <- maf_pc_for_window(window_data)
        if (!is.null(res)) {
          c(time = time_vec[i], maf1 = res$maf1[i], pc1 = res$pc1[i])
        } else {
          c(time = time_vec[i], maf1 = NA_real_, pc1 = NA_real_)
        }
      }
    )
    results_df <- as.data.frame(do.call(base::rbind, results_list))
  }
  # Standardize and detect warnings
  for (col in c("maf1", "pc1")) {
    if (col %in% names(results_df)) {
      scaled_vals <- as.numeric(scale(results_df[[col]]))
      results_df[[paste0(col, "_scaled")]]  <- scaled_vals
      results_df[[paste0(col, "_warning")]] <- abs(scaled_vals) > threshold
    }
  }
  results_df
}
