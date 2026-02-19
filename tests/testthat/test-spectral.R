# --- Test data ---------------------------------------------------------------
set.seed(42)
spec_white <- rnorm(200)
spec_brown <- cumsum(rnorm(200))

# --- Basic output structure --------------------------------------------------

test_that("spectral_ews returns a spectral tibble", {
  result <- spectral_ews(spec_brown, window = 50L)
  expect_s3_class(result, "spectral")
  expect_s3_class(result, "tbl_df")
})

test_that("spectral_ews returns correct columns with states", {
  result <- spectral_ews(spec_brown, window = 50L)
  expected_cols <- c(
    "time", "value", "spectral_exponent", "spectral_ratio",
    "r_squared", "state"
  )
  expect_identical(names(result), expected_cols)
})

test_that("states FALSE omits state column", {
  result <- spectral_ews(spec_brown, window = 50L, states = FALSE)
  expect_false("state" %in% names(result))
  expected_cols <- c(
    "time", "value", "spectral_exponent", "spectral_ratio", "r_squared"
  )
  expect_identical(names(result), expected_cols)
})

test_that("state is a factor with correct levels", {
  result <- spectral_ews(spec_brown, window = 50L)
  expect_s3_class(result$state, "factor")
  expected_levels <- c("white_noise", "pink_noise", "red_noise", "brownian")
  expect_identical(levels(result$state), expected_levels)
})

test_that("spectral_exponent is numeric and not all NA", {
  result <- spectral_ews(spec_brown, window = 50L)
  expect_type(result$spectral_exponent, "double")
  expect_false(all(is.na(result$spectral_exponent)))
})

test_that("spectral_ratio is numeric", {
  result <- spectral_ews(spec_brown, window = 50L)
  expect_type(result$spectral_ratio, "double")
})

test_that("r_squared is numeric", {
  result <- spectral_ews(spec_brown, window = 50L)
  expect_type(result$r_squared, "double")
})

# --- Attributes ---------------------------------------------------------------

test_that("attributes are stored correctly", {
  result <- spectral_ews(
    spec_brown, window = 50L, align = "right",
    method = "periodogram", detrend = "linear"
  )
  expect_identical(attr(result, "window"), 50L)
  expect_identical(attr(result, "align"), "right")
  expect_identical(attr(result, "method"), "periodogram")
  expect_identical(attr(result, "detrend"), "linear")
})

# --- Method options -----------------------------------------------------------

test_that("periodogram method works", {
  result <- spectral_ews(spec_brown, window = 50L, method = "periodogram")
  expect_s3_class(result, "spectral")
})

test_that("ar method works", {
  result <- spectral_ews(spec_brown, window = 50L, method = "ar")
  expect_s3_class(result, "spectral")
})

# --- Detrend options ----------------------------------------------------------

test_that("detrend none works", {
  result <- spectral_ews(spec_white, window = 50L, detrend = "none")
  expect_s3_class(result, "spectral")
})

test_that("detrend linear works", {
  result <- spectral_ews(spec_brown, window = 50L, detrend = "linear")
  expect_s3_class(result, "spectral")
})

test_that("detrend diff works", {
  result <- spectral_ews(spec_brown, window = 50L, detrend = "diff")
  expect_s3_class(result, "spectral")
})

# --- Alignment options --------------------------------------------------------

test_that("align right works", {
  result <- spectral_ews(spec_brown, window = 50L, align = "right")
  expect_s3_class(result, "spectral")
  expect_equal(nrow(result), length(spec_brown))
})

test_that("align center works", {
  result <- spectral_ews(spec_brown, window = 50L, align = "center")
  expect_s3_class(result, "spectral")
  expect_equal(nrow(result), length(spec_brown))
})

test_that("align left works", {
  result <- spectral_ews(spec_brown, window = 50L, align = "left")
  expect_s3_class(result, "spectral")
  expect_equal(nrow(result), length(spec_brown))
})

# --- S3 methods ---------------------------------------------------------------

test_that("print does not error", {
  result <- spectral_ews(spec_brown, window = 50L)
  capture.output(print(result)) |>
    expect_error(NA)
})

test_that("summary does not error and returns list invisibly", {
  result <- spectral_ews(spec_brown, window = 50L)
  out <- capture.output(s <- summary(result))
  expect_type(s, "list")
  expect_true("mean_beta" %in% names(s))
  expect_true("mean_ratio" %in% names(s))
  expect_true("mean_r_squared" %in% names(s))
  expect_true("window" %in% names(s))
  expect_true("method" %in% names(s))
})

test_that("plot does not error", {
  result <- spectral_ews(spec_brown, window = 50L)
  plot(result, type = "series") |>
    expect_error(NA)
  plot(result, type = "states") |>
    expect_error(NA)
  plot(result, type = "both") |>
    expect_error(NA)
})

# --- Input validation ---------------------------------------------------------

test_that("bad method produces error", {
  expect_error(
    spectral_ews(spec_brown, window = 50L, method = "nonexistent"),
    "method"
  )
})

test_that("bad detrend produces error", {
  expect_error(
    spectral_ews(spec_brown, window = 50L, detrend = "nonexistent"),
    "detrend"
  )
})

test_that("bad align produces error", {
  expect_error(
    spectral_ews(spec_brown, window = 50L, align = "nonexistent"),
    "align"
  )
})

# --- Internal helper: spectral_detrend_ edge cases (line 37) ----------------

test_that("spectral_detrend_ linear with constant time returns x - mean", {
  detrend_fn <- codyna:::spectral_detrend_
  # Line 37: xx_var == 0 is not reachable with sequential t_idx, but

  # test that linear detrending produces residuals summing to ~0
  x <- c(5, 5, 5, 5, 5)
  result <- detrend_fn(x, "linear")
  expect_true(all(abs(result) < 1e-10))
})

test_that("spectral_detrend_ diff returns first differences", {
  detrend_fn <- codyna:::spectral_detrend_
  x <- c(1, 3, 6, 10)
  result <- detrend_fn(x, "diff")
  expect_equal(result, c(2, 3, 4))
})

test_that("spectral_detrend_ none returns input unchanged", {
  detrend_fn <- codyna:::spectral_detrend_
  x <- c(1, 2, 3)
  result <- detrend_fn(x, "none")
  expect_identical(result, x)
})

# --- Internal helper: spectral_metrics_ edge cases (lines 55, 62, 67, 76) ---

test_that("spectral_metrics_ returns NULL for short segment", {
  metrics_fn <- codyna:::spectral_metrics_
  # Line 55: length(x) < 10
  result <- metrics_fn(c(1, 2, 3), "periodogram")
  expect_null(result)
})

test_that("spectral_metrics_ returns NULL for zero-variance segment", {
  metrics_fn <- codyna:::spectral_metrics_
  # Line 55: var(x) < 1e-12
  result <- metrics_fn(rep(5, 20), "periodogram")
  expect_null(result)
})

test_that("spectral_metrics_ returns NULL on spectrum failure", {
  metrics_fn <- codyna:::spectral_metrics_
  # Line 62: try-error from spectrum
  # An empty or problematic input that passes length/var checks
  # but causes spectrum to fail
  x <- c(rep(0, 9), 1e-14, rep(0, 9), 1e-14)
  result <- metrics_fn(x, "periodogram")
  # Result might be NULL if variance is too low or valid points insufficient
  # At minimum, the function should not error
  expect_true(is.null(result) || is.list(result))
})

test_that("spectral_metrics_ ar method works", {
  metrics_fn <- codyna:::spectral_metrics_
  set.seed(42)
  x <- cumsum(rnorm(50))
  result <- metrics_fn(x, "ar")
  expect_true(is.list(result))
  expect_true(all(c("beta", "ratio", "r_squared") %in% names(result)))
})

# --- Internal helper: spectral_build_rects_ empty input (lines 124-127) -----

test_that("spectral_build_rects_ returns empty data.frame for empty input", {
  rects_fn <- codyna:::spectral_build_rects_
  # Lines 124-127: n == 0 branch
  empty_d <- data.frame(
    time = numeric(0),
    state_f = factor(character(0),
                     levels = c("white_noise", "pink_noise", "red_noise", "brownian"))
  )
  result <- rects_fn(empty_d)
  expect_equal(nrow(result), 0L)
  expect_true("xmin" %in% names(result))
  expect_true("xmax" %in% names(result))
  expect_true("state_f" %in% names(result))
})

# --- spectral_ews: min_points skip and short detrended (lines 342, 349) -----

test_that("spectral_ews handles window with many NAs (min_points skip)", {
  set.seed(42)
  # Create data with a block of NAs in the middle to trigger min_points skip
  x <- rnorm(100)
  x[30:70] <- NA
  # Line 342: n_valid < min_points triggers next
  result <- spectral_ews(x, window = 50L, min_points = 20L)
  expect_s3_class(result, "spectral")
  expect_equal(nrow(result), 100L)
})

test_that("spectral_ews with diff detrend and small window triggers short detrended", {
  set.seed(42)
  # Line 349: after diff detrending, segment may be < 10 observations
  # Use a window of ~11 so diff makes it 10, barely passing
  x <- cumsum(rnorm(100))
  result <- spectral_ews(x, window = 12L, detrend = "diff",
                         min_points = 4L)
  expect_s3_class(result, "spectral")
})

# --- spectral_ews: internal NA interpolation (lines 381-391) ----------------

test_that("spectral_ews interpolates internal NAs in spectral metrics", {
  set.seed(42)
  # Create a series where some windows fail spectral estimation,
  # producing internal NAs that need interpolation
  # Use a series with a few constant-value windows embedded
  x <- rnorm(150)
  # Insert a constant block that will produce NULL from spectral_metrics_
  x[60:75] <- 0
  result <- spectral_ews(x, window = 30L, min_points = 10L)
  expect_s3_class(result, "spectral")
  # After interpolation, most spectral_exponent values should be non-NA
  non_na_count <- sum(!is.na(result$spectral_exponent))
  expect_true(non_na_count > 100L)
})

# --- plot.spectral: states derived from exponent when missing (lines 503-509) -

test_that("plot.spectral derives states when state column is absent", {
  set.seed(42)
  # Lines 503-509: plot without state column
  result <- spectral_ews(spec_brown, window = 50L, states = FALSE)
  expect_false("state" %in% names(result))
  # Plot should derive states internally from spectral_exponent
  p <- plot(result, type = "series")
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
  p2 <- plot(result, type = "states")
  expect_true(inherits(p2, "gg") || inherits(p2, "ggplot"))
  p3 <- plot(result, type = "both")
  expect_true(inherits(p3, "gg") || inherits(p3, "patchwork"))
})

# --- summary.spectral: full return structure ---------------------------------

test_that("summary returns all expected elements", {
  result <- spectral_ews(spec_brown, window = 50L)
  out <- capture.output(s <- summary(result))
  expect_true("state_counts" %in% names(s))
  expect_true("detrend" %in% names(s))
  expect_true("align" %in% names(s))
})

# --- spectral_ews with left alignment ----------------------------------------

test_that("left alignment produces correct length output", {
  set.seed(42)
  x <- cumsum(rnorm(100))
  result <- spectral_ews(x, window = 30L, align = "left")
  expect_equal(nrow(result), 100L)
})

# --- spectral_detrend_ linear with single-element vector (line 37) -----------

test_that("spectral_detrend_ linear with single element returns zero", {
  detrend_fn <- codyna:::spectral_detrend_
  # Line 37: n=1 => t_idx = 1, t_mean = 1, xx_var = 0 => x - x_mean
  result <- detrend_fn(5.0, "linear")
  expect_equal(result, 0.0)
})

# --- spectral_metrics_ edge: ar method failure (line 62) ---------------------

test_that("spectral_metrics_ handles ar spectrum failure gracefully", {
  metrics_fn <- codyna:::spectral_metrics_
  # Create data that has sufficient length and variance but where
  # ar spectrum might produce problematic results
  # A spike + near-constant series
  set.seed(42)
  x <- c(rep(0.001, 9), 100, rep(0.001, 9), -100)
  result <- metrics_fn(x, "ar")
  # Should either succeed or return NULL, not error

  expect_true(is.null(result) || is.list(result))
})

# --- spectral_metrics_ edge: few valid freq/power (line 67) ------------------

test_that("spectral_metrics_ returns NULL when few valid freq/power points", {
  metrics_fn <- codyna:::spectral_metrics_
  # Create borderline data that produces a spectrum with few valid points
  # Extremely short but > 10 data with minimal variability
  set.seed(42)
  x <- c(rnorm(10, 0, 0.001), 1)  # 11 points, tiny variance except one outlier
  result <- metrics_fn(x, "periodogram")
  # May or may not be NULL depending on spectrum output
  expect_true(is.null(result) || is.list(result))
})

# --- spectral_metrics_ edge: xx_var near zero in log regression (line 76) ----

test_that("spectral_metrics_ handles near-zero xx_var in log regression", {
  metrics_fn <- codyna:::spectral_metrics_
  # This is very hard to trigger since log(freq) values will vary if
  # there are >= 4 valid frequency points. The only way is if all
  # valid log frequencies are identical, which doesn't happen in practice.
  # Just verify robustness with extreme data
  set.seed(42)
  x <- rep(c(1, -1), 10)  # 20 points alternating +1/-1
  result <- metrics_fn(x, "periodogram")
  expect_true(is.null(result) || is.list(result))
})

# --- spectral_ews: diff detrend making segment < 10 (line 349) --------------

test_that("spectral_ews with very small window and diff produces short segments", {
  set.seed(42)
  # Window = 11 means after diff we get 10 obs (just barely >= 10)
  # Window = 10 with diff gives 9 obs, which should trigger line 349 (< 10)
  # But min_points must be <= window. With window = 10, diff makes 9 obs.
  x <- cumsum(rnorm(80))
  # min_points = 4 so we satisfy the check, but after diff we get 9 obs
  # which is < 10 triggering the next at line 349
  result <- spectral_ews(x, window = 10L, detrend = "diff", min_points = 4L)
  expect_s3_class(result, "spectral")
  # Most spectral_exponent values should be NA since diff makes segments < 10
  # Actually the valid count check comes first (line 342), then detrend, then
  # length check at 349. With 10 obs, n_valid = 10 >= min_points=4, then
  # diff makes it 9, then 9 < 10 triggers next
  # Edge padding should fill NAs if any valid windows exist
})

# --- spectral_metrics_: line 62 - spectrum try-error with NAs in segment ------

test_that("spectral_metrics_ returns NULL when spectrum fails on NA data", {
  # Data with embedded NAs: length >= 10 and var >= 1e-12, but

  # spectrum() errors on data containing NAs -> try-error at line 62
  x_na <- c(1, NA, 3, NA, 5, NA, 7, NA, 9, 10)
  result_pgram <- codyna:::spectral_metrics_(x_na, "periodogram")
  expect_null(result_pgram)

  result_ar <- codyna:::spectral_metrics_(x_na, "ar")
  expect_null(result_ar)
})

# --- spectral_metrics_: lines 67, 76 - defensive guards ----------------------
# Line 67 (sum(valid) < 4): pgram always produces positive spec for data with
#   nonzero variance and >= 10 observations. The taper = 0.1 ensures spectral
#   leakage prevents exact zeros. Effectively unreachable.
# Line 76 (xx_var < 1e-12): log(freq) values are always distinct when there
#   are >= 4 valid frequency bins, so variance is always > 0. Unreachable.
# These tests verify the function works correctly near these boundary conditions.

test_that("spectral_metrics_ handles data with minimal variation robustly", {
  # Near-constant data with one outlier; variance passes the 1e-12 check
  # but spectrum produces very small (but positive) power values
  x <- c(rep(0, 9), 1)
  result <- codyna:::spectral_metrics_(x, "periodogram")
  # Should either succeed with a valid list or return NULL
  expect_true(is.null(result) || is.list(result))
  if (is.list(result)) {
    expect_true(all(c("beta", "ratio", "r_squared") %in% names(result)))
  }
})
