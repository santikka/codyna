

#' Detect Changepoints in Time Series Data
#'
#' Identifies abrupt structural changes in a time series using one of three
#' algorithms: CUSUM (cumulative sum), binary segmentation, or PELT (Pruned
#' Exact Linear Time). The function detects shifts in mean, variance, or
#' both simultaneously, and returns a segmented time series with segment
#' statistics and changepoint indicators.
#'
#' @details
#' The three detection algorithms differ in computational strategy and
#' optimality guarantees:
#'
#' **CUSUM** (`method = "cusum"`). The cumulative sum of deviations from
#' the segment mean is computed as \eqn{S_k = \sum_{i=1}^{k}(x_i - \bar{x})}.
#' The point of maximum \eqn{|S_k|} is the candidate changepoint. If
#' splitting at this point reduces the penalised cost, the segment is split
#' and the algorithm recurses on both halves (binary segmentation driven by
#' CUSUM). This is fast and intuitive but yields approximate solutions.
#'
#' **Binary segmentation** (`method = "binseg"`). At each
#' recursion level, every valid split position within the segment is
#' evaluated via the likelihood ratio (cost of the combined segment minus
#' cost of the two sub-segments). The position maximising this gain is
#' selected. If the gain exceeds the penalty, the segment is split and
#' the algorithm recurses. This is more thorough than CUSUM-based splitting
#' but still approximate (greedy).
#'
#' **PELT** (`method = "pelt"`). A dynamic programming algorithm that
#' finds the exact optimal segmentation under the penalised cost criterion.
#' Pruning ensures that candidate changepoints that can provably never
#' be part of the optimal solution are discarded, giving an expected
#' \eqn{O(n)} running time for data with a linear number of changepoints
#' (Killick et al., 2012). PELT is recommended for most applications
#' because it combines exactness with efficiency.
#'
#' **Cost function.** All methods use a Gaussian log-likelihood cost.
#' For `type = "mean"`, the cost is \eqn{n \log(\hat{\sigma}^2)} where
#' \eqn{\hat{\sigma}^2} is the segment variance. For `type = "variance"`,
#' an additional constant term is added. For `type = "both"`, the full
#' normal log-likelihood is used, detecting simultaneous shifts in mean
#' and variance.
#'
#' **Penalty selection.** The penalty controls the trade-off between fit
#' and complexity. `"bic"` (\eqn{\log n}) provides a consistent estimator
#' of the true number of changepoints; `"aic"` (penalty = 2) is more
#' liberal and may over-segment; `"manual"` allows a user-specified value
#' via the `penalty_value` parameter.
#'
#' @export
#' @param data \[`numeric` or `ts`]\cr
#'   The time series to analyze. Accepts a numeric vector or a `ts` object.
#' @param method \[`character(1)`: `"cusum"`]\cr
#'   Changepoint detection algorithm. One of `"cusum"`,
#'   `"binseg"`, or `"pelt"`.
#' @param penalty \[`character(1)`: `"bic"`]\cr
#'   Penalty type for model selection. One of `"bic"`, `"aic"`, or
#'   `"manual"`.
#' @param min_segment \[`integer(1)`: `10L`]\cr
#'   Minimum number of observations between consecutive changepoints.
#'   Prevents spurious detections in short segments.
#' @param max_changepoints \[`integer(1)` or `NULL`: `NULL`]\cr
#'   Maximum number of changepoints to detect. If `NULL`, no limit is
#'   imposed (the penalty controls complexity).
#' @param type \[`character(1)`: `"mean"`]\cr
#'   Type of change to detect. One of `"mean"` (changes in mean only),
#'   `"variance"` (changes in variance only), or `"both"` (simultaneous
#'   changes in mean and variance).
#' @param penalty_value \[`numeric(1)` or `NULL`: `NULL`]\cr
#'   Manual penalty value. Only used when `penalty = "manual"`.
#' @return An object of class `"changepoint"` (inheriting from
#'   [tibble::tibble()]) with the following columns:
#'
#'     * `time` -- the time index.
#'     * `value` -- the original time series values.
#'     * `segment` -- integer sequential segment ID (1, 2, 3, ...).
#'     * `regime` -- integer smart regime ID. Segments with similar
#'       means share the same regime number (e.g., 1, 2, 1, 2 if the
#'       series alternates between two levels).
#'     * `segment_mean` -- mean of the segment containing this
#'       observation.
#'     * `segment_var` -- variance of the segment containing this
#'       observation.
#'     * `level` -- factor (`"low"`, `"medium"`, `"high"`). Classifies
#'       each segment's mean relative to the overall series: high if
#'       above mean + 0.5 SD, low if below mean - 0.5 SD, medium
#'       otherwise.
#'     * `direction` -- factor (`"lower"`, `"no_change"`, `"higher"`).
#'       Direction of the mean shift from the previous segment. `NA` for
#'       the first segment.
#'     * `magnitude` -- numeric difference from the previous segment
#'       mean. `NA` for the first segment.
#'     * `changepoint_type` -- character (`"change"` or `"return"`).
#'       `"change"` marks a transition to a never-before-seen regime;
#'       `"return"` marks a transition back to a previously visited
#'       regime. Only populated at changepoint rows (`changepoint = TRUE`);
#'       `NA` elsewhere.
#'     * `state` -- factor combining level and changepoint type
#'       (e.g., `"high_change"`, `"low_return"`, `"medium_initial"`).
#'       The first segment is always `"{level}_initial"`.
#'     * `changepoint` -- logical; `TRUE` at detected changepoint
#'       locations (the first observation of each new segment after the
#'       first).
#'
#'   The following attributes are stored on the returned object:
#'   `method`, `penalty`, `min_segment`, `max_changepoints`, `type`,
#'   `n_changepoints`, `changepoint_locations`.
#'
#' @references
#' Page, E. S. (1954). Continuous inspection schemes.
#' \emph{Biometrika}, 41(1/2), 100--115. \doi{10.2307/2333009}
#'
#' Scott, A. J. & Knott, M. (1974). A cluster analysis method for
#' grouping means in the analysis of variance. \emph{Biometrics},
#' 30(3), 507--512. \doi{10.2307/2529204}
#'
#' Killick, R., Fearnhead, P., & Eckley, I. A. (2012). Optimal detection
#' of changepoints with a linear computational cost. \emph{Journal of the
#' American Statistical Association}, 107(500), 1590--1598.
#' \doi{10.1080/01621459.2012.737745}
#'
#' @seealso [compute_trend()] for rolling trend classification;
#'   [resilience()] for rolling resilience metrics.
#' @examples
#' \donttest{
#' # Mean shift at t = 200
#' set.seed(42)
#' x <- c(rnorm(200, 0, 1), rnorm(200, 5, 1))
#' cpts <- detect_cpts(x, method = "pelt")
#' cpts
#' summary(cpts)
#' plot(cpts)
#'
#' # Multiple changepoints
#' set.seed(123)
#' y <- c(rnorm(100, 0), rnorm(100, 3), rnorm(100, -1), rnorm(100, 4))
#' cp2 <- detect_cpts(y, method = "binseg")
#' plot(cp2, type = "both")
#' }
detect_cpts <- function(data, method = "cusum", penalty = "bic",
                        min_segment = 10L, max_changepoints = NULL,
                        type = "mean", penalty_value = NULL) {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)
  method <- check_match(method, c("cusum", "binseg", "pelt"))
  penalty <- check_match(penalty, c("bic", "aic", "manual"))
  type <- check_match(type, c("mean", "variance", "both"))
  check_range(min_segment, type = "integer", min = 2L, max = floor(n / 2))
  if (!is.null(max_changepoints)) {
    check_range(
      max_changepoints,
      type = "integer",
      min = 1L,
      max = floor(n / min_segment) - 1L
    )
  }
  stopifnot_(
    penalty != "manual" || !is.null(penalty_value),
    "When {.arg penalty} is {.val manual},
    {.arg penalty_value} must be provided."
  )
  pen <- changepoint_penalty_(penalty, penalty_value, n, type)
  cpts <- switch(
    method,
    cusum = changepoint_cusum_(
      values, type, pen, min_segment, max_changepoints, 0L
    ),
    binseg = changepoint_binseg_(
      values, type, pen, min_segment, max_changepoints, 0L
    ),
    pelt = changepoint_pelt_(
      values, type, pen, min_segment, max_changepoints
    )
  )
  cpts <- sort(unique(as.integer(cpts)))
  if (length(cpts) > 0L) {
    valid <- cpts >= min_segment & cpts <= (n - min_segment)
    cps <- cpts[valid]
  }
  # Remove changepoints too close to each other
  # nocov start -- all algorithms already enforce min_segment spacing
  # if (length(cps) > 1L) {
  #   keep <- rep(TRUE, length(cps))
  #   for (i in seq(2L, length(cps))) {
  #     if ((cps[i] - cps[i - 1L]) < min_segment) {
  #       keep[i] <- FALSE
  #     }
  #   }
  #   cps <- cps[keep]
  # }
  # nocov end
  n_cpt <- length(cpts)
  seg_starts <- c(1L, cpts + 1L)
  seg_ends <- c(cpts, n)
  n_segments <- length(seg_starts)
  segment <- integer(n)
  seg_mean <- numeric(n)
  seg_var <- numeric(n)
  for (s in seq_len(n_segments)) {
    idx <- seg_starts[s]:seg_ends[s]
    segment[idx] <- s
    seg_vals <- values[idx]
    seg_mean[idx] <- mean(seg_vals, na.rm = TRUE)
    seg_var[idx] <- ifelse_(
      length(seg_vals) > 1L,
      stats::var(seg_vals, na.rm = TRUE),
      0
    )
  }
  cpt_logical <- rep(FALSE, n)
  if (n_cpt > 0L) {
    cpt_indicator_pos <- cpts + 1L
    cpt_indicator_pos <- cpt_indicator_pos[cpt_indicator_pos <= n]
    cpt_logical[cpt_indicator_pos] <- TRUE
  }
  unique_seg_means <- vapply(
    seq_len(n_segments),
    function(s) {
      mean(values[seg_starts[s]:seg_ends[s]], na.rm = TRUE)
    },
    numeric(1L)
  )
  cls <- changepoint_classify_(
    unique_seg_means,
    mean(values, na.rm = TRUE),
    stats::sd(values, na.rm = TRUE)
  )
  regime    <- cls$regime[segment]
  level     <- factor(cls$level[segment], levels = c("low", "medium", "high"))
  direction <- factor(
    cls$direction[segment],
    levels = c("lower", "no_change", "higher"),
    exclude = NULL
  )
  magnitude <- cls$magnitude[segment]
  cpt_type   <- cls$changepoint_type[segment]
  # changepoint_type only at actual changepoint rows, NA elsewhere
  cpt_type[!cpt_logical] <- NA_character_
  state <- factor(cls$state[segment], levels = unique(cls$state))
  out <- tibble::tibble(
    time = time,
    value = values,
    segment = segment,
    regime = regime,
    segment_mean = seg_mean,
    segment_var = seg_var,
    level = level,
    direction = direction,
    magnitude = magnitude,
    changepoint_type = cpt_type,
    state = state,
    changepoint = cpt_logical
  )
  structure(
    out,
    method = method,
    penalty = penalty,
    min_segment = min_segment,
    max_changepoints = max_changepoints,
    type = type,
    n_changepoints = n_cpt,
    changepoint_locations = cpts,
    class = c("changepoint", "tbl_df", "tbl", "data.frame")
  )
}

#' Compute segment cost under normal log-likelihood
#'
#' Returns the negative log-likelihood cost of modelling `x` as a single
#' segment. The cost depends on `type`: for "mean" only the variance term
#' matters; for "variance" the mean is assumed known; for "both" the full
#' normal log-likelihood is used.
#'
#' @param x Numeric vector (segment data).
#' @param type One of "mean", "variance", "both".
#' @return Single numeric cost value.
#' @noRd
changepoint_cost_ <- function(x, type) {
  n <- length(x)
  if (n < 2L) {
    return(0)
  }
  mu <- mean(x)
  switch(
    type,
    mean = {
      # RSS: residual sum of squares from segment mean.
      # Matches changepoint::cpt.mean() cost function.
      sum((x - mu)^2)
    },
    variance = {
      # n * log(var_MLE): Gaussian log-likelihood for variance detection.
      # Matches changepoint::cpt.var() cost function.
      v <- sum((x - mu)^2) / n
      if (v < .Machine$double.eps) v <- .Machine$double.eps
      n * log(v)
    },
    both = {
      # Full Normal log-likelihood: n * log(var_MLE) + n.
      # Matches changepoint::cpt.meanvar() cost function.
      v <- sum((x - mu)^2) / n
      if (v < .Machine$double.eps) v <- .Machine$double.eps
      n * log(v) + n
    },
    {
      sum((x - mu)^2)
    }
  )
}

#' Compute penalty value from penalty specification
#'
#' @param penalty Character: "bic", "aic", or "manual".
#' @param penalty_value Numeric value for manual penalty (ignored otherwise).
#' @param n Integer series length.
#' @return Single numeric penalty value.
#' @noRd
changepoint_penalty_ <- function(penalty, penalty_value, n, type) {
  # Number of parameters per changepoint: mean or variance = 2, both = 3
  # Matches the changepoint package (Killick et al.) penalty computation
  p <- ifelse_(type == "both", 3, 2)
  switch(
    penalty,
    bic = p * log(n),
    aic = 2 * p,
    manual = {
      if (is.null(penalty_value) || !is.numeric(penalty_value)) {
        stop_(
          "When {.arg penalty} is {.val manual},
          {.arg penalty_value} must be a numeric value."
        )
      }
      penalty_value
    },
    p * log(n)
  )
}

#' Find the single best changepoint in a segment using likelihood ratio
#'
#' Scans all candidate split positions within `x` and returns the position
#' that maximises the likelihood ratio (cost reduction). Returns NULL if no
#' valid split exists.
#'
#' @param x Numeric vector.
#' @param type Cost type.
#' @param min_segment Minimum segment length.
#' @return List with `pos` (best split position relative to x) and
#'   `gain` (cost improvement), or NULL.
#' @noRd
changepoint_best_split_ <- function(x, type, min_segment) {
  n <- length(x)
  if (n < 2L * min_segment) {
    return(NULL)
  }
  cost_full <- changepoint_cost_(x, type)
  best_gain <- -Inf
  best_pos <- NULL
  candidates <- seq(min_segment, n - min_segment)
  for (tau in candidates) {
    left <- x[seq_len(tau)]
    right <- x[(tau + 1L):n]
    cost_split <- changepoint_cost_(left, type) +
      changepoint_cost_(right, type)
    gain <- cost_full - cost_split
    if (gain > best_gain) {
      best_gain <- gain
      best_pos <- tau
    }
  }
  list(pos = best_pos, gain = best_gain)
}

#' CUSUM statistic for a segment
#'
#' Computes cumulative sum of deviations from the segment mean.
#' Returns the index of the maximum absolute CUSUM value.
#'
#' @param x Numeric vector.
#' @return Integer index of maximum |CUSUM|.
#' @noRd
changepoint_cusum_stat_ <- function(x) {
  n <- length(x)
  mu <- mean(x)
  cusum <- cumsum(x - mu)
  which.max(abs(cusum))
}

#' CUSUM-based changepoint detection via binary segmentation
#'
#' Uses CUSUM to identify the most likely changepoint, then tests whether
#' the split is worthwhile via a penalty-based likelihood ratio test.
#' Recurses on both halves if the split is accepted.
#'
#' @param x Numeric vector.
#' @param type Cost type.
#' @param pen Penalty value.
#' @param min_segment Minimum segment length.
#' @param max_cp Maximum changepoints remaining (NULL = unlimited).
#' @param offset Integer offset for converting local indices to global.
#' @return Integer vector of global changepoint positions.
#' @noRd
changepoint_cusum_ <- function(x, type, pen, min_segment, max_cp, offset) {
  n <- length(x)
  if (n < 2L * min_segment) {
    return(integer(0))
  }
  if (!is.null(max_cp) && max_cp <= 0L) {
    return(integer(0))
  }
  # Find CUSUM-based candidate
  tau <- changepoint_cusum_stat_(x)
  if (tau < min_segment || (n - tau) < min_segment) {
    return(integer(0))
  }
  # Check if split is worthwhile
  cost_full <- changepoint_cost_(x, type)
  left <- x[seq_len(tau)]
  right <- x[(tau + 1L):n]
  cost_split <- changepoint_cost_(left, type) +
    changepoint_cost_(right, type)
  gain <- cost_full - cost_split
  if (gain <= pen) return(integer(0))
  global_pos <- tau + offset
  new_max <- onlyif(!is.null(max_cp), max_cp - 1L)
  left_cp <- changepoint_cusum_(
    left, type, pen, min_segment, new_max, offset
  )
  if (!is.null(new_max)) {
    new_max <- new_max - length(left_cp)
  }
  right_cp <- changepoint_cusum_(
    right, type, pen, min_segment, new_max, offset + tau
  )
  sort(c(left_cp, global_pos, right_cp))
}

#' Binary segmentation changepoint detection
#'
#' Iteratively finds the best single split and recurses if the gain
#' exceeds the penalty.
#'
#' @param x Numeric vector.
#' @param type Cost type.
#' @param pen Penalty value.
#' @param min_segment Minimum segment length.
#' @param max_cp Maximum changepoints remaining.
#' @param offset Integer offset for global indexing.
#' @return Integer vector of global changepoint positions.
#' @noRd
changepoint_binseg_ <- function(x, type, pen, min_segment, max_cp, offset) {
  n <- length(x)
  if (n < 2L * min_segment) {
    return(integer(0))
  }
  if (!is.null(max_cp) && max_cp <= 0L) {
    return(integer(0))
  }
  result <- changepoint_best_split_(x, type, min_segment)
  if (result$gain <= pen) return(integer(0))
  tau <- result$pos
  global_pos <- tau + offset
  new_max <- onlyif(!is.null(max_cp), max_cp - 1L)
  left <- x[seq_len(tau)]
  left_cp <- changepoint_binseg_(
    left, type, pen, min_segment, new_max, offset
  )
  if (!is.null(new_max)) {
    new_max <- new_max - length(left_cp)
  }
  right <- x[(tau + 1L):n]
  right_cp <- changepoint_binseg_(
    right, type, pen, min_segment, new_max, offset + tau
  )
  sort(c(left_cp, global_pos, right_cp))
}

#' PELT changepoint detection
#'
#' Pruned Exact Linear Time algorithm. Uses dynamic programming with
#' pruning to find the exact optimal segmentation under a penalised
#' cost criterion.
#'
#' @param x Numeric vector.
#' @param type Cost type.
#' @param pen Penalty value.
#' @param min_segment Minimum segment length.
#' @param max_cp Maximum number of changepoints (NULL = unlimited).
#' @return Integer vector of changepoint positions.
#' @noRd
changepoint_pelt_ <- function(x, type, pen, min_segment, max_cpt) {

  n <- length(x)
  if (n < 2L * min_segment) {
    return(integer(0))
  }
  # F[t+1] = optimal cost for data 1..t (0-indexed: F[1] = cost for 0 data)
  f_cost <- rep(Inf, n + 1L)
  f_cost[1L] <- -pen  # F(0) = -beta so F(0) + C(0+1..t) + beta = C(1..t)
  cpt_trace <- vector("list", n + 1L)
  cpt_trace[[1L]] <- integer(0)
  # R: set of candidate last changepoints (0-indexed positions)
  candidates <- 0L
  for (t_star in seq_len(n)) {
    best_cost <- Inf
    best_tau <- 0L
    # First pass: find optimal cost
    cost_vals <- numeric(length(candidates))
    valid_mask <- logical(length(candidates))
    for (ci in seq_along(candidates)) {
      tau <- candidates[ci]
      seg_start <- tau + 1L
      seg_end <- t_star
      seg_len <- seg_end - seg_start + 1L
      if (seg_len < min_segment) {
        valid_mask[ci] <- FALSE
        cost_vals[ci] <- Inf
        next
      }
      valid_mask[ci] <- TRUE
      seg <- x[seg_start:seg_end]
      cost_vals[ci] <- f_cost[tau + 1L] + changepoint_cost_(seg, type) + pen
      if (cost_vals[ci] < best_cost) {
        best_cost <- cost_vals[ci]
        best_tau <- tau
      }
    }
    f_cost[t_star + 1L] <- best_cost
    # Reconstruct changepoint list
    prev_cpts <- cpt_trace[[best_tau + 1L]]
    if (best_tau > 0L) {
      cpt_trace[[t_star + 1L]] <- c(prev_cpts, best_tau)
    } else {
      cpt_trace[[t_star + 1L]] <- prev_cpts
    }
    # PELT pruning: keep candidates where F(tau) + C(tau+1..t) + pen <= F(t)
    # Candidates with segments too short are always kept (may become valid later)
    prune_keep <- integer(0)
    for (ci in seq_along(candidates)) {
      if (!valid_mask[ci] || cost_vals[ci] <= best_cost + pen) {
        prune_keep <- c(prune_keep, candidates[ci])
      }
    }
    # Add current t_star as a candidate
    candidates <- unique(c(prune_keep, t_star))
  }
  cpts <- cpt_trace[[n + 1L]]
  if (length(cpts) == 0L) {
    return(integer(0))
  }
  cpts <- sort(unique(cpts))
  # Enforce max_changepoints
  if (!is.null(max_cpt) && length(cpts) > max_cpt) {
    gains <- vapply(
      cpts,
      function(cpt) {
        # Estimate gain by removing this changepoint
        seg_before_start <- ifelse_(
          cpt == cpts[1L],
          1L,
          cpts[which(cpts == cpt) - 1L] + 1L
        )
        seg_after_end <- ifelse_(
          cpt == cpts[length(cpts)],
          n,
          cpts[which(cpts == cpt) + 1L]
        )
        left <- x[seg_before_start:cpt]
        right <- x[(cpt + 1L):seg_after_end]
        combined <- x[seg_before_start:seg_after_end]
        changepoint_cost_(combined, type) -
          (changepoint_cost_(left, type) + changepoint_cost_(right, type))
      },
      numeric(1)
    )
    keep_idx <- order(gains, decreasing = TRUE)[seq_len(max_cpt)]
    cpts <- sort(cpts[keep_idx])
  }
  cpts
}

#' Color palette for changepoint regime states
#'
#' Named colors for the 9 possible level x changepoint_type combinations.
#' @noRd
changepoint_colors_ <- function() {
  c(
    high_initial = "#E53935",
    high_change = "#C62828",
    high_return = "#EF5350",
    medium_initial = "#FFA726",
    medium_change = "#F57C00",
    medium_return = "#FFCC80",
    low_initial = "#43A047",
    low_change = "#2E7D32",
    low_return = "#66BB6A"
  )
}

#' Classify changepoint segments into meaningful regime labels
#'
#' Groups segments with similar means into shared regimes, classifies
#' each segment by level (high/medium/low) and direction (higher/lower),
#' and labels changepoints as "change" (new regime) or "return" (revisit).
#'
#' @param seg_means Numeric vector of per-segment means.
#' @param overall_mean Numeric, global series mean.
#' @param overall_sd Numeric, global series SD.
#' @return List with: regime (integer), level (character),
#'   direction (character), magnitude (numeric),
#'   changepoint_type (character), state (character).
#'   Each vector has length equal to number of segments.
#' @noRd
changepoint_classify_ <- function(seg_means, overall_mean, overall_sd) {
  n <- length(seg_means)
  threshold <- 0.5 * overall_sd
  regime <- integer(n)
  regime_means <- numeric(0)
  regime[1L] <- 1L
  regime_means[1L] <- seg_means[1L]
  if (n > 1L) {
    for (i in 2:n) {
      dists <- abs(seg_means[i] - regime_means)
      closest <- which.min(dists)
      if (dists[closest] <= threshold) {
        regime[i] <- closest
        regime_means[closest] <- mean(
          seg_means[regime[seq_len(i)] == closest]
        )
      } else {
        new_id <- length(regime_means) + 1L
        regime[i] <- new_id
        regime_means[new_id] <- seg_means[i]
      }
    }
  }
  level <- vapply(
    seg_means,
    function(m) {
      if (m > overall_mean + 0.5 * overall_sd) {
        return("high")
      }
      if (m < overall_mean - 0.5 * overall_sd) {
        return("low")
      }
      "medium"
    },
    character(1L)
  )
  direction <- rep(NA_character_, n)
  magnitude <- rep(NA_real_, n)
  if (n > 1L) {
    for (i in 2:n) {
      mag <- seg_means[i] - seg_means[i - 1L]
      magnitude[i] <- mag
      direction[i] <- if (mag > 0.1 * overall_sd) {
        "higher"
      } else if (mag < -0.1 * overall_sd) {
        "lower"
      } else {
        "no_change"
      }
    }
  }
  changepoint_type <- rep(NA_character_, n)
  seen_regimes <- regime[1L]
  if (n > 1L) {
    for (i in 2:n) {
      if (regime[i] %in% seen_regimes) {
        changepoint_type[i] <- "return"
      } else {
        changepoint_type[i] <- "change"
        seen_regimes <- c(seen_regimes, regime[i])
      }
    }
  }
  cpt_for_state <- changepoint_type
  cpt_for_state[is.na(cpt_for_state)] <- "initial"
  state <- paste0(level, "_", cpt_for_state)
  list(
    regime = regime,
    level = level,
    direction = direction,
    magnitude = magnitude,
    changepoint_type = changepoint_type,
    state = state
  )
}
