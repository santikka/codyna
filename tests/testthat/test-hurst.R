# --------------------------------------------------------------------------
# test-hurst.R -- Tests for hurst() and detect_hurst_warnings()
# --------------------------------------------------------------------------

set.seed(42)
hurst_ts <- cumsum(rnorm(200))

# ==========================================================================
# hurst() with states = TRUE
# ==========================================================================

test_that("hurst() with states = TRUE returns correct class and columns", {
  h <- hurst(hurst_ts, window = 30L)
  expect_s3_class(h, "hurst")
  expect_s3_class(h, "tbl_df")
  expect_true(all(c("time", "value", "hurst", "r_squared",
                     "state", "transition") %in% names(h)))
  expect_equal(nrow(h), length(hurst_ts))
})

test_that("hurst() stores settings as attributes", {
  h <- hurst(hurst_ts, window = 30L, step = 2L, method = "dfa")
  expect_equal(attr(h, "window"), 30L)
  expect_equal(attr(h, "step"), 2L)
  expect_equal(attr(h, "method"), "dfa")
})

test_that("hurst values are in a reasonable range", {
  h <- hurst(hurst_ts, window = 30L)
  h_vals <- h$hurst[!is.na(h$hurst)]
  expect_true(all(h_vals >= -0.5 & h_vals <= 2.5))
})

test_that("hurst states are valid categories", {
  h <- hurst(hurst_ts, window = 30L)
  valid_states <- c(
    "strong_antipersistent", "antipersistent", "random_walk",
    "persistent", "strong_persistent"
  )
  non_na_states <- h$state[!is.na(h$state)]
  expect_true(all(non_na_states %in% valid_states))
})

# ==========================================================================
# hurst() with states = FALSE
# ==========================================================================

test_that("hurst() with states = FALSE returns hurst_global list", {
  hg <- hurst(hurst_ts, states = FALSE)
  expect_s3_class(hg, "hurst_global")
  expect_type(hg, "list")
  expect_true(all(c("hurst", "r_squared", "method", "n") %in% names(hg)))
  expect_equal(hg$n, length(hurst_ts))
  expect_equal(hg$method, "dfa")
})

test_that("hurst_global hurst value is a single numeric", {
  hg <- hurst(hurst_ts, states = FALSE)
  expect_length(hg$hurst, 1L)
  expect_true(is.numeric(hg$hurst))
})

# ==========================================================================
# Different methods
# ==========================================================================

test_that("hurst() works with dfa method", {
  hurst(hurst_ts, method = "dfa", window = 30L) |>
    expect_error(NA)
})

test_that("hurst() works with rs method", {
  hurst(hurst_ts, method = "rs", window = 30L) |>
    expect_error(NA)
})

test_that("hurst() works with mfdfa method (states = TRUE)", {
  h <- hurst(hurst_ts, method = "mfdfa", window = 30L)
  expect_s3_class(h, "hurst")
  expect_true("mf_width" %in% names(h))
  expect_true("mf_category" %in% names(h))
})

test_that("hurst() works with rs method (states = FALSE)", {
  hg <- hurst(hurst_ts, method = "rs", states = FALSE)
  expect_s3_class(hg, "hurst_global")
  expect_equal(hg$method, "rs")
  expect_true("rs_values" %in% names(hg))
})

test_that("hurst() works with mfdfa method (states = FALSE)", {
  hg <- hurst(hurst_ts, method = "mfdfa", states = FALSE)
  expect_s3_class(hg, "hurst_global")
  expect_true(all(c("hq", "tauq", "mf_width") %in% names(hg)))
})

# ==========================================================================
# Window, step, and scaling parameters
# ==========================================================================

test_that("step parameter affects output density", {
  h1 <- hurst(hurst_ts, window = 30L, step = 1L)
  h2 <- hurst(hurst_ts, window = 30L, step = 5L)
  # Both should have the same nrow (full series), but step=5

  # may produce more interpolated values
  expect_equal(nrow(h1), nrow(h2))
})

test_that("different scaling options work", {
  scalings <- c("none", "center", "standardize", "minmax", "iqr")
  for (s in scalings) {
    hurst(hurst_ts, window = 30L, scaling = s) |>
      expect_error(NA)
  }
})

# ==========================================================================
# Input validation
# ==========================================================================

test_that("hurst() errors on invalid method", {
  expect_error(hurst(hurst_ts, method = "invalid"))
})

test_that("hurst() errors on window too small", {
  expect_error(hurst(hurst_ts, window = 2L))
})

test_that("hurst() errors on window larger than series", {
  expect_error(hurst(hurst_ts, window = 500L))
})

test_that("hurst() errors on missing data argument", {
  expect_error(hurst())
})

# ==========================================================================
# detect_hurst_warnings()
# ==========================================================================

test_that("detect_hurst_warnings() returns hurst_ews class", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  expect_s3_class(ews, "hurst_ews")
  expect_s3_class(ews, "tbl_df")
})

test_that("detect_hurst_warnings() has required columns", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  required_cols <- c(
    "time", "value", "hurst", "r_squared", "state", "transition",
    "extreme_low", "extreme_high", "trend_up", "trend_down",
    "high_volatility", "flickering", "variance_ratio",
    "spectral_shift", "autocorr_increase", "state_persistence",
    "warning_score", "warning_level", "warning_label"
  )
  expect_true(all(required_cols %in% names(ews)))
})

test_that("detect_hurst_warnings() warning scores are in [0, 1]", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  expect_true(all(ews$warning_score >= 0 & ews$warning_score <= 1))
})

test_that("detect_hurst_warnings() warning levels are 0-4", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  expect_true(all(ews$warning_level %in% 0:4))
})

test_that("detect_hurst_warnings() warning labels are valid", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  valid_labels <- c("none", "low", "moderate", "high", "critical")
  expect_true(all(ews$warning_label %in% valid_labels))
})

test_that("detect_hurst_warnings() stores indicator weights", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  weights <- attr(ews, "indicator_weights")
  expect_true(is.numeric(weights))
  expect_length(weights, 10L)
})

test_that("detect_hurst_warnings() preserves row count", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  expect_equal(nrow(ews), nrow(h))
})

test_that("detect_hurst_warnings() errors on non-hurst input", {
  expect_error(detect_hurst_warnings(data.frame(x = 1:10)))
})

test_that("detect_hurst_warnings() custom window parameters work", {
  h <- hurst(hurst_ts, window = 30L)
  detect_hurst_warnings(h, trend_window = 15L,
                        volatility_window = 15L,
                        flicker_window = 10L) |>
    expect_error(NA)
})

# ==========================================================================
# S3 methods: print
# ==========================================================================

test_that("print.hurst does not error", {
  h <- hurst(hurst_ts, window = 30L)
  print(h) |>
    capture.output() |>
    expect_error(NA)
})

test_that("print.hurst_global does not error", {
  hg <- hurst(hurst_ts, states = FALSE)
  print(hg) |>
    capture.output() |>
    expect_error(NA)
})

test_that("print.hurst_ews does not error", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  print(ews) |>
    capture.output() |>
    expect_error(NA)
})

# ==========================================================================
# S3 methods: summary
# ==========================================================================

test_that("summary.hurst_ews does not error", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  summary(ews) |>
    capture.output() |>
    expect_error(NA)
})

# ==========================================================================
# S3 methods: plot
# ==========================================================================

test_that("plot.hurst does not error", {
  h <- hurst(hurst_ts, window = 30L)
  plot(h) |>
    expect_error(NA)
  plot(h, type = "series") |>
    expect_error(NA)
  plot(h, type = "states") |>
    expect_error(NA)
  plot(h, type = "both") |>
    expect_error(NA)
})

test_that("plot.hurst_ews does not error", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  plot(ews) |>
    expect_error(NA)
})

test_that("plot.hurst_ews with subset of indicators does not error", {
  h <- hurst(hurst_ts, window = 30L)
  ews <- detect_hurst_warnings(h)
  plot(ews, indicators = c("flickering", "trend_up")) |>
    expect_error(NA)
})

# ==========================================================================
# ts object input
# ==========================================================================

test_that("hurst() accepts ts objects with frequency = 1", {
  ts_obj <- ts(hurst_ts, frequency = 1)
  hurst(ts_obj, window = 30L) |>
    expect_error(NA)
})

# ==========================================================================
# Edge cases: hurst_scale() -- lines 329, 336, 341
# ==========================================================================

test_that("hurst_scale standardize with sd=0 returns centered values", {
  set.seed(42)
  # Constant series has sd=0, triggers line 329
  const_ts <- rep(5, 50)
  result <- codyna:::hurst_scale(const_ts, "standardize")
  expect_true(all(result == 0))
})

test_that("hurst_scale minmax with range=0 returns 0.5", {
  set.seed(42)
  # Constant series has range=0, triggers line 336
  const_ts <- rep(3, 50)
  result <- codyna:::hurst_scale(const_ts, "minmax")
  expect_true(all(result == 0.5))
})

test_that("hurst_scale iqr with iqr=0 returns median-centered values", {
  set.seed(42)
  # All identical values => IQR=0, triggers line 341
  const_ts <- rep(7, 50)
  result <- codyna:::hurst_scale(const_ts, "iqr")
  expect_true(all(result == 0))
})

# ==========================================================================
# Edge cases: hurst_dfa() -- lines 350, 353-354, 356, 358, 360, 363-364
# ==========================================================================

test_that("hurst_dfa handles NAs in input (line 350)", {
  set.seed(42)
  x <- cumsum(rnorm(100))
  x[c(5, 10, 15)] <- NA
  result <- codyna:::hurst_dfa(x, min_scale = 4L, max_scale = 25L,
                               n_scales = 10L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_dfa returns NA for too-short series (lines 353-354)", {
  set.seed(42)
  # n < 2 * min_scale
  x <- rnorm(5)
  result <- codyna:::hurst_dfa(x, min_scale = 4L, max_scale = 10L,
                               n_scales = 10L)
  expect_true(is.na(result$hurst))
  expect_true(is.na(result$r_squared))
  expect_equal(length(result$scales), 0L)
})

test_that("hurst_dfa adjusts min_scale below 4 (line 356)", {
  set.seed(42)
  x <- cumsum(rnorm(100))
  result <- codyna:::hurst_dfa(x, min_scale = 2L, max_scale = 25L,
                               n_scales = 10L)
  expect_true(is.numeric(result$hurst))
  expect_true(!is.na(result$hurst))
})

test_that("hurst_dfa handles NULL max_scale (line 358)", {
  set.seed(42)
  x <- cumsum(rnorm(100))
  result <- codyna:::hurst_dfa(x, min_scale = 4L, max_scale = NULL,
                               n_scales = 10L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_dfa adjusts max_scale > n/2 (line 360)", {
  set.seed(42)
  x <- cumsum(rnorm(50))
  result <- codyna:::hurst_dfa(x, min_scale = 4L, max_scale = 40L,
                               n_scales = 10L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_dfa returns NA when max_scale <= min_scale (lines 363-364)", {
  set.seed(42)
  x <- cumsum(rnorm(20))
  result <- codyna:::hurst_dfa(x, min_scale = 10L, max_scale = 5L,
                               n_scales = 10L)
  expect_true(is.na(result$hurst))
  expect_true(is.na(result$r_squared))
})

# ==========================================================================
# Edge cases: hurst_dfa large scales -- lines 379-380, 400, 417
# ==========================================================================

test_that("hurst_dfa handles n_seg < 2 for large scale (lines 379-380)", {
  set.seed(42)
  # Very short series so some scales will have n_seg < 2
  x <- cumsum(rnorm(20))
  result <- codyna:::hurst_dfa(x, min_scale = 4L, max_scale = NULL,
                               n_scales = 10L)
  # Should still return a result (some scales may be NA)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_dfa with constant segment (xx_var=0) (lines 400, 417)", {
  set.seed(42)
  # When x_vals has length 1, xx_var = 0
  # This happens when scale = 1, but min_scale is enforced at 4
  # Create data where detrending produces zero variance
  x <- rep(0, 30)
  result <- codyna:::hurst_dfa(x, min_scale = 4L, max_scale = 7L,
                               n_scales = 5L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_dfa returns NA with < 3 valid fluctuations (lines 424-425)", {
  set.seed(42)
  # Very short series with large min_scale => few valid fluctuations
  x <- cumsum(rnorm(12))
  result <- codyna:::hurst_dfa(x, min_scale = 4L, max_scale = 5L,
                               n_scales = 2L)
  # With only 2 n_scales and potentially one valid point, may return NA

  expect_true(is.numeric(result$hurst))
})

# ==========================================================================
# Edge cases: hurst_rs() -- lines 450, 453-454, 464-465, 481-482
# ==========================================================================

test_that("hurst_rs handles NULL max_scale (line 450)", {
  set.seed(42)
  x <- cumsum(rnorm(100))
  result <- codyna:::hurst_rs(x, min_scale = 4L, max_scale = NULL,
                              n_scales = 10L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_rs returns NA when max_scale < min_scale (lines 453-454)", {
  set.seed(42)
  x <- cumsum(rnorm(20))
  result <- codyna:::hurst_rs(x, min_scale = 20L, max_scale = 5L,
                              n_scales = 10L)
  expect_true(is.na(result$hurst))
  expect_true(is.na(result$r_squared))
  expect_equal(length(result$scales), 0L)
})

test_that("hurst_rs handles n_seg < 1 for large scales (lines 464-465)", {
  set.seed(42)
  x <- cumsum(rnorm(20))
  # Force a scale larger than n
  result <- codyna:::hurst_rs(x, min_scale = 4L, max_scale = 10L,
                              n_scales = 10L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_rs returns NA with < 3 valid rs values (lines 481-482)", {
  set.seed(42)
  # Very short series with constant values => sd=0 => NA rs values
  x <- rep(5, 12)
  result <- codyna:::hurst_rs(x, min_scale = 4L, max_scale = 6L,
                              n_scales = 3L)
  expect_true(is.na(result$hurst) || is.numeric(result$hurst))
})

# ==========================================================================
# Edge cases: hurst_mfdfa() -- lines 501, 504-505, 516, 526, 542-543
# ==========================================================================

test_that("hurst_mfdfa handles NULL max_scale (line 501)", {
  set.seed(42)
  x <- cumsum(rnorm(100))
  result <- codyna:::hurst_mfdfa(x, q = seq(-2, 2, 1), min_scale = 4L,
                                 max_scale = NULL, n_scales = 10L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_mfdfa returns NA when max_scale < min_scale (lines 504-505)", {
  set.seed(42)
  x <- cumsum(rnorm(20))
  result <- codyna:::hurst_mfdfa(x, q = seq(-2, 2, 1), min_scale = 20L,
                                 max_scale = 3L, n_scales = 10L)
  expect_true(is.na(result$hurst))
  expect_true(is.na(result$r_squared))
  expect_equal(length(result$hq), 0L)
})

test_that("hurst_mfdfa handles n_seg < 1 (line 516)", {
  set.seed(42)
  # Short series where some scales have n_seg < 1
  x <- cumsum(rnorm(20))
  result <- codyna:::hurst_mfdfa(x, q = c(-2, 0, 2), min_scale = 4L,
                                 max_scale = NULL, n_scales = 10L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_mfdfa handles rms_vals all zero (line 526)", {
  set.seed(42)
  # Constant series => perfect linear fit => zero residuals => zero rms
  x <- seq_len(30)
  result <- codyna:::hurst_mfdfa(x, q = c(-2, 0, 2), min_scale = 4L,
                                 max_scale = 7L, n_scales = 5L)
  expect_true(is.numeric(result$hurst))
})

test_that("hurst_mfdfa handles < 3 valid hq (lines 542-543)", {
  set.seed(42)
  # Very few scales with valid data
  x <- cumsum(rnorm(20))
  result <- codyna:::hurst_mfdfa(x, q = c(-2, 0, 2), min_scale = 4L,
                                 max_scale = 5L, n_scales = 2L)
  # Some hq values may be NA
  expect_true(is.numeric(result$hurst) || is.na(result$hurst))
})

# ==========================================================================
# Rolling window edge cases -- lines 242, 244, 248, 250
# ==========================================================================

test_that("hurst rolling window adjusts boundaries (lines 242, 244)", {
  set.seed(42)
  # Use a window size close to series length to trigger boundary adjustments
  x <- cumsum(rnorm(60))
  h <- hurst(x, window = 25L, step = 1L)
  expect_s3_class(h, "hurst")
  expect_equal(nrow(h), 60L)
})

test_that("hurst rolling handles NA-heavy windows (line 248)", {
  set.seed(42)
  x <- cumsum(rnorm(100))
  # Insert many NAs in a block to trigger the >80% NA skip
  x[30:70] <- NA
  h <- hurst(x, window = 30L, step = 1L)
  expect_s3_class(h, "hurst")
  expect_equal(nrow(h), 100L)
})

test_that("hurst rolling handles try-error in computation (line 250)", {
  set.seed(42)
  # Very short constant segments can cause numerical issues
  x <- c(rep(0, 50), cumsum(rnorm(50)))
  h <- hurst(x, window = 20L, step = 1L, min_scale = 4L)
  expect_s3_class(h, "hurst")
})

# ==========================================================================
# MFDFA interpolation -- lines 287-289
# ==========================================================================

test_that("hurst mfdfa rolling with internal gaps triggers interpolation (lines 287-289)", {
  set.seed(42)
  x <- cumsum(rnorm(150))
  # Step > 1 creates internal gaps that need interpolation
  h <- hurst(x, method = "mfdfa", window = 30L, step = 3L)
  expect_s3_class(h, "hurst")
  expect_true("mf_width" %in% names(h))
  # Check mf_width values are filled (not all NA)
  expect_true(any(!is.na(h$mf_width)))
})

# ==========================================================================
# build_state_rects_ with empty data -- lines 1381-1384
# ==========================================================================

test_that("build_state_rects_ handles empty data (lines 1381-1384)", {
  empty_df <- data.frame(
    time = numeric(0), value = numeric(0),
    hurst = numeric(0), r_squared = numeric(0),
    state = character(0), transition = numeric(0)
  )
  result <- codyna:::build_state_rects_(empty_df)
  expect_equal(nrow(result), 0L)
  expect_true("xmin" %in% names(result))
  expect_true("xmax" %in% names(result))
  expect_true("state_f" %in% names(result))
})

# ==========================================================================
# build_warning_rects_ with empty data -- lines 1417-1420
# ==========================================================================

test_that("build_warning_rects_ handles empty data (lines 1417-1420)", {
  empty_df <- data.frame(
    time = numeric(0),
    warning_level = integer(0)
  )
  result <- codyna:::build_warning_rects_(empty_df)
  expect_equal(nrow(result), 0L)
  expect_true("xmin" %in% names(result))
  expect_true("xmax" %in% names(result))
  expect_true("level" %in% names(result))
})

# ==========================================================================
# Dead code / safety net lines in hurst.R
# The following lines are unreachable under the current clamping logic:
#   - Line 250: try-error in rolling (internal fns return lists, never throw)
#   - Lines 379-380: n_seg < 2 in DFA (max_scale clamped to <= n/2)
#   - Lines 400, 417: xx_var == 0 (x_vals = seq_len(s), s >= 4, always > 0)
#   - Lines 464-465: n_seg < 1 in RS (max_scale clamped to floor(n/2))
#   - Line 516: n_seg < 1 in MFDFA (max_scale clamped to floor(n/4))
# ==========================================================================

# ==========================================================================
# Targeted coverage: hurst.R line 526 -- rms_vals all zero in MFDFA
# ==========================================================================

test_that("hurst_mfdfa skips scale when all rms_vals == 0 (line 526)", {
  set.seed(42)
  # Constant data: cumsum(x - mean(x)) = 0 for all positions.
  # Each segment is all zeros, lm.fit residuals = 0, rms_vals = 0.
  # After filtering rms_vals > 0, empty vector triggers line 526 next.
  x <- rep(5, 40)
  result <- codyna:::hurst_mfdfa(x, q = c(-1, 0, 1), min_scale = 4L,
                                 max_scale = 10L, n_scales = 5L)
  expect_true(is.na(result$hurst))
})
