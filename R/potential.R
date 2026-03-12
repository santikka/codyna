#' Stability Landscape and Potential Analysis
#'
#' Estimates the quasi-potential landscape \eqn{U(x) = -\log P(x)} from
#' kernel density estimation, identifying stable attractors (potential
#' wells) and unstable equilibria (barriers). Supports global analysis of
#' the full series or rolling-window analysis to track how the number and
#' location of attractors change over time.
#'
#' @details
#' The potential function is derived from the Boltzmann relation between
#' probability density and energy: \eqn{P(x) \propto \exp(-U(x)/D)},
#' where \eqn{D} is a noise intensity parameter. Setting \eqn{D = 1}
#' gives \eqn{U(x) = -\log P(x)}. The probability density \eqn{P(x)} is
#' estimated via `stats::density()` with a Gaussian kernel.
#'
#' **Wells** (local minima of \eqn{U}) correspond to stable attractors --
#' states the system tends to occupy. **Barriers** (local maxima of
#' \eqn{U}) correspond to unstable equilibria separating attractors.
#'
#' Well depth measures how much "energy" is needed to escape an attractor:
#' deeper wells are more stable. Well width measures the basin of
#' attraction in state space. Together, depth and width characterise the
#' resilience of each attractor -- deep, wide wells indicate robust
#' stability, while shallow, narrow wells are easily escaped by
#' perturbations.
#'
#' **Detrending.** Non-stationary trends can bias the density estimate and
#' create spurious features. The `detrend` parameter offers three options:
#'
#' * `"none"` (default): No detrending. Appropriate for stationary
#'     or nearly stationary data.
#' * `"linear"`: Removes a linear trend via OLS before estimating
#'     the landscape. Appropriate when the data has a slow drift.
#' * `"diff"`: Uses first differences \eqn{\Delta x_t = x_t -
#'     x_{t-1}}. Appropriate for integrated processes (random walks).
#'     Note: reduces series length by one.
#'
#'
#' **Rolling-window mode.** When `window` is specified, the landscape is
#' re-estimated at each window position and the number of wells and the
#' location of the dominant (deepest) well are tracked over time. Changes
#' in the number of wells -- especially transitions from one to two wells
#' or vice versa -- can signal an approaching bifurcation (tipping point).
#'
#' @export
#' @param data \[`numeric` or `ts`\]\cr
#'   The time series to analyze. Accepts a numeric vector or a `ts` object.
#' @param window \[`integer(1)`: `NULL`\]\cr
#'   Rolling window size. If `NULL` (default), a single global landscape
#'   is estimated from the entire series. If specified, the landscape is
#'   re-estimated for each window position.
#' @param n_bins \[`integer(1)`: `50L`\]\cr
#'   Number of bins for the kernel density estimate grid. The actual grid
#'   resolution passed to `stats::density()` is `n_bins * 4` to ensure
#'   smooth interpolation.
#' @param bandwidth \[`numeric(1)`: `NULL`\]\cr
#'   Bandwidth for the kernel density estimate. If `NULL`, Silverman's
#'   rule of thumb is used (the `stats::density()` default).
#' @param detrend \[`character(1)`: `"none"`\]\cr
#'   Detrending method applied before estimating the landscape. One of
#'   `"none"`, `"linear"`, or `"diff"`.
#'
#' @return An object of class `"potential"` (a list) containing:
#'
#' * `landscape`: \[`tibble`\] Grid of the estimated landscape with
#'   columns `x` (state-space position), `density` (probability
#'   density), and `potential` (\eqn{-\log(\mathrm{density})}).
#'   Present only for global analysis (`window = NULL`).
#' * `wells`: \[`tibble`\] One row per detected well with columns
#'   `location` (x position), `depth` (barrier height minus well
#'   bottom), and `width` (distance between adjacent barriers).
#'   Present only for global analysis.
#' * `barriers`: \[`tibble`\] One row per detected barrier with
#'   columns `location` (x position) and `height` (potential value).
#'   Present only for global analysis.
#' * `n_wells`: \[`integer(1)`\] Number of detected wells.
#' * `values`: \[`numeric()`\] The (possibly detrended) data used
#'   for the analysis.
#' * `time`: \[`numeric()`\] The time index.
#' * `rolling`: \[`tibble`\] Present only for rolling analysis
#'   (`window` specified). Contains columns `time`, `value`,
#'   `n_wells`, and `dominant_well_location`.
#'
#'   Attributes: `n_bins`, `bandwidth`, `detrend`, `window`.
#'
#' @references
#' Livina, V.N., Kwasniok, F., & Lenton, T.M. (2010). Potential analysis
#' reveals changing number of climate states during the last 60 kyr.
#' \emph{Climate of the Past}, 6(1), 77--82.
#' \doi{10.5194/cp-6-77-2010}
#'
#' Scheffer, M., Carpenter, S.R., Lenton, T.M., Bascompte, J.,
#' Brock, W., Dakos, V., van de Koppel, J., van de Leemput, I.A.,
#' Levin, S.A., van Nes, E.H., Pascual, M., & Vandermeer, J. (2012).
#' Anticipating critical transitions.
#' \emph{Science}, 338(6105), 344--348.
#' \doi{10.1126/science.1225244}
#'
#' Strogatz, S.H. (2015). \emph{Nonlinear Dynamics and Chaos} (2nd ed.).
#' Westview Press.
#'
#' @seealso [resilience()] for rolling resilience metrics;
#'   [hurst()] for long-range dependence analysis;
#'   [compute_trend()] for trend classification.
#' @examples
#' \donttest{
#' # Single well: Ornstein-Uhlenbeck process (attractor at 0)
#' set.seed(42)
#' n <- 2000
#' x <- numeric(n)
#' for (i in 2:n) {
#'   x[i] <- x[i-1] - 0.5 * x[i-1] * 0.01 + 0.3 * rnorm(1)
#' }
#' pa <- potential_analysis(x)
#' pa
#' plot(pa)
#'
#' # Two wells: bimodal mixture
#' set.seed(42)
#' bimodal <- c(rnorm(1000, -2, 0.5), rnorm(1000, 2, 0.5))
#' pa2 <- potential_analysis(bimodal)
#' summary(pa2)
#' plot(pa2, type = "both")
#' }
potential_analysis <- function(data, window = NULL, n_bins = 50L,
                               bandwidth = NULL, detrend = "none") {
  check_missing(data)
  data <- prepare_timeseries_data(data)
  values <- data$values
  time <- data$time
  n <- length(values)
  detrend <- check_match(detrend, c("none", "linear", "diff"))
  check_range(n_bins, type = "integer", min = 10L, max = 10000L)
  if (!is.null(bandwidth)) {
    check_range(bandwidth, type = "numeric", scalar = TRUE, min = 0)
  }
  detrended <- potential_detrend_(values, detrend)
  if (detrend == "diff") {
    time <- time[-1L]
  }
  n_d <- length(detrended)
  if (!is.null(window)) {
    check_range(window, type = "integer", min = 30L, max = n_d)
  }
  if (is.null(window)) {
    result <- potential_landscape_(detrended, n_bins, bandwidth)
    out <- list(
      landscape = result$landscape,
      wells = result$wells,
      barriers = result$barriers,
      n_wells = result$n_wells,
      values = detrended,
      time = time
    )
  } else {
    n_positions <- n_d - window + 1L
    rolling_n_wells <- integer(n_d)
    rolling_dominant <- rep(NA_real_, n_d)
    rolling_n_wells[seq_len(window - 1L)] <- NA_integer_
    for (i in seq_len(n_positions)) {
      end_idx <- i + window - 1L
      w_data <- detrended[i:end_idx]
      w_clean <- w_data[!is.na(w_data)]
      if (length(w_clean) < 20L) {
        rolling_n_wells[end_idx] <- NA_integer_
        next
      }
      w_result <- try_(potential_landscape_(w_clean, n_bins, bandwidth))
      if (inherits(w_result, "try-error") || is.null(w_result)) {
        rolling_n_wells[end_idx] <- NA_integer_
        next
      }
      rolling_n_wells[end_idx] <- w_result$n_wells
      # Dominant well = deepest well
      if (nrow(w_result$wells) > 0L) {
        best <- which.max(w_result$wells$depth)
        rolling_dominant[end_idx] <- w_result$wells$location[best]
      }
    }
    rolling_tbl <- tibble::tibble(
      time = time,
      value = detrended,
      n_wells = rolling_n_wells,
      dominant_well_location = rolling_dominant
    )
    # Also compute a global landscape for reference
    global_result <- potential_landscape_(detrended, n_bins, bandwidth)
    out <- list(
      landscape = global_result$landscape,
      wells = global_result$wells,
      barriers = global_result$barriers,
      n_wells = global_result$n_wells,
      values = detrended,
      time = time,
      rolling = rolling_tbl
    )
  }
  structure(
    out,
    n_bins = n_bins,
    bandwidth = bandwidth,
    detrend = detrend,
    window = window,
    class = c("potential", "list")
  )
}

#' Find local minima in a numeric vector
#'
#' Identifies indices where the first difference changes from negative to
#' positive (valley). Boundary points are excluded.
#'
#' @param x Numeric vector.
#' @return Integer vector of indices.
#' @noRd
potential_find_minima_ <- function(x) {
  n <- length(x)
  if (n < 3L) {
    return(integer(0))
  }
  d <- diff(x)
  # Sign changes: negative -> positive
  sign_d <- sign(d)
  # Remove exact zeros by forward-filling
  for (i in seq_along(sign_d)) {
    if (sign_d[i] == 0L && i > 1L) {
      sign_d[i] <- sign_d[i - 1L]
    }
  }
  idx <- which(sign_d[-length(sign_d)] < 0 & sign_d[-1L] > 0) + 1L
  idx
}

#' Find local maxima in a numeric vector
#'
#' Identifies indices where the first difference changes from positive to
#' negative (peak). Boundary points are excluded.
#'
#' @param x Numeric vector.
#' @return Integer vector of indices.
#' @noRd
potential_find_maxima_ <- function(x) {
  n <- length(x)
  if (n < 3L) {
    return(integer(0))
  }
  d <- diff(x)
  sign_d <- sign(d)
  for (i in seq_along(sign_d)) {
    if (sign_d[i] == 0L && i > 1L) {
      sign_d[i] <- sign_d[i - 1L]
    }
  }
  idx <- which(sign_d[-length(sign_d)] > 0 & sign_d[-1L] < 0) + 1L
  idx
}

#' Compute potential landscape from a numeric vector
#'
#' Estimates kernel density, derives potential U(x) = -log(density), and
#' locates wells and barriers.
#'
#' @param x Numeric vector (cleaned, no NAs).
#' @param n_bins Integer grid resolution.
#' @param bandwidth Numeric bandwidth or NULL for Silverman default.
#' @return List with landscape, wells, barriers, n_wells.
#' @noRd
potential_landscape_ <- function(x, n_bins, bandwidth) {
  n_grid <- n_bins * 4L
  dens_args <- list(x = x, n = n_grid)
  if (!is.null(bandwidth)) {
    dens_args$bw <- bandwidth
  }
  dens <- do.call(stats::density, dens_args)
  # Potential: U(x) = -log(P(x)), with epsilon to avoid log(0)
  eps <- .Machine$double.eps
  density_vals <- pmax(dens$y, eps)
  potential_vals <- -log(density_vals)
  landscape <- tibble::tibble(
    x = dens$x,
    density = dens$y,
    potential = potential_vals
  )
  # Find wells (local minima of potential) and barriers (local maxima)
  well_idx <- potential_find_minima_(potential_vals)
  barrier_idx <- potential_find_maxima_(potential_vals)
  # Build wells tibble
  if (length(well_idx) == 0L) {
    # If no interior minimum found, use the global minimum
    global_min <- which.min(potential_vals)
    well_idx <- global_min
  }
  # Compute well properties
  well_locations <- dens$x[well_idx]
  well_depths <- numeric(length(well_idx))
  well_widths <- numeric(length(well_idx))
  # For each well, find depth and width relative to nearest barriers
  # Build ordered set of barrier + boundary positions for depth/width
  all_barriers <- sort(unique(c(1L, barrier_idx, length(potential_vals))))
  for (wi in seq_along(well_idx)) {
    w <- well_idx[wi]
    # Find left barrier (closest barrier index < w)
    left_candidates <- all_barriers[all_barriers < w]
    left_b <- ifelse_(
      length(left_candidates) > 0L,
      left_candidates[length(left_candidates)],
      1L
    )
    # Find right barrier (closest barrier index > w)
    right_candidates <- all_barriers[all_barriers > w]
    right_b <- ifelse_(
      length(right_candidates) > 0L,
      right_candidates[1L],
      length(potential_vals)
    )
    # Depth: min of (left_barrier_height - well_bottom, right_barrier_height - well_bottom)
    well_bottom <- potential_vals[w]
    left_height <- potential_vals[left_b]
    right_height <- potential_vals[right_b]
    well_depths[wi] <- min(left_height, right_height) - well_bottom
    # Width: distance in x between left and right barriers
    well_widths[wi] <- dens$x[right_b] - dens$x[left_b]
  }
  wells <- tibble::tibble(
    location = well_locations,
    depth = well_depths,
    width = well_widths
  )
  # Build barriers tibble
  if (length(barrier_idx) > 0L) {
    barriers <- tibble::tibble(
      location = dens$x[barrier_idx],
      height = potential_vals[barrier_idx]
    )
  } else {
    barriers <- tibble::tibble(
      location = numeric(0),
      height = numeric(0)
    )
  }
  list(
    landscape = landscape,
    wells = wells,
    barriers = barriers,
    n_wells = length(well_idx)
  )
}

#' Detrend a numeric vector
#'
#' @param x Numeric vector.
#' @param method One of "none", "linear", "diff".
#' @return Numeric vector (detrended).
#' @noRd
potential_detrend_ <- function(x, method) {
  switch(
    method,
    none = x,
    linear = {
      t_idx <- seq_along(x)
      fit <- stats::lm.fit(x = cbind(1, t_idx), y = x)
      fit$residuals
    },
    diff = {
      diff(x)
    }
  )
}
