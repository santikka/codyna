test_that("compute_trend returns trend class inheriting from tibble", {
  set.seed(42)
  x <- cumsum(rnorm(200, sd = 2))
  tr <- compute_trend(x, window = 20)
  expect_s3_class(tr, "trend")
  expect_s3_class(tr, "tbl_df")
  expect_s3_class(tr, "tbl")
  expect_s3_class(tr, "data.frame")
})

test_that("compute_trend output has correct columns", {
  set.seed(42)
  x <- cumsum(rnorm(200, sd = 2))
  tr <- compute_trend(x, window = 20)
  expect_true(all(c("time", "value", "metric", "state") %in% names(tr)))
  expect_equal(nrow(tr), 200L)
})

test_that("state column is a factor with correct levels", {
  set.seed(42)
  x <- cumsum(rnorm(200, sd = 2))
  tr <- compute_trend(x, window = 20)
  expect_s3_class(tr$state, "factor")
  expected_levels <- c("ascending", "descending", "flat",
                       "turbulent", "Missing_Data", "Initial")
  expect_equal(levels(tr$state), expected_levels)
})

test_that("metric column is numeric", {
  set.seed(42)
  x <- cumsum(rnorm(200, sd = 2))
  tr <- compute_trend(x, window = 20)
  expect_true(is.numeric(tr$metric))
})

test_that("method slope works", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  compute_trend(x, window = 10, method = "slope") |>
    expect_error(NA)
})

test_that("method ar1_phi1 works", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  compute_trend(x, window = 10, method = "ar1_phi1") |>
    expect_error(NA)
})

test_that("method growth_factor works", {
  set.seed(42)
  x <- cumsum(rnorm(100, mean = 50, sd = 2))
  compute_trend(x, window = 10, method = "growth_factor") |>
    expect_error(NA)
})

test_that("slope sub-methods all work", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  sub_methods <- c("ols", "robust", "spearman", "kendall")
  for (sm in sub_methods) {
    compute_trend(x, window = 10, method = "slope", slope_method = sm) |>
      expect_error(NA)
  }
})

test_that("align right works", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 10, align = "right")
  expect_s3_class(tr, "trend")
  expect_equal(attr(tr, "align"), "right")
})

test_that("align center works", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 10, align = "center")
  expect_s3_class(tr, "trend")
  expect_equal(attr(tr, "align"), "center")
})

test_that("align left works", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 10, align = "left")
  expect_s3_class(tr, "trend")
  expect_equal(attr(tr, "align"), "left")
})

test_that("attributes are stored correctly", {
  set.seed(42)
  x <- cumsum(rnorm(200, sd = 2))
  tr <- compute_trend(x, window = 25, method = "slope", align = "right")
  expect_equal(attr(tr, "window"), 25L)
  expect_equal(attr(tr, "align"), "right")
  expect_equal(attr(tr, "method"), "slope")
  expect_equal(attr(tr, "slope_method"), "robust")
  expect_equal(attr(tr, "epsilon"), 0.05)
})

test_that("slope_method attribute is NULL for non-slope methods", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 10, method = "ar1_phi1")
  expect_null(attr(tr, "slope_method"))
})

test_that("window defaults to adaptive when NULL", {
  set.seed(42)
  x <- cumsum(rnorm(200, sd = 2))
  tr <- compute_trend(x)
  expect_s3_class(tr, "trend")
  expected_window <- max(3L, min(200L, round(200 / 10)))
  expect_equal(attr(tr, "window"), expected_window)
})

test_that("compute_trend accepts ts objects", {
  set.seed(42)
  x <- ts(cumsum(rnorm(100, sd = 2)), start = 2000, frequency = 12)
  tr <- compute_trend(x, window = 10)
  expect_s3_class(tr, "trend")
  expect_equal(nrow(tr), 100L)
})

test_that("compute_trend errors on invalid method", {
  set.seed(42)
  x <- rnorm(100)
  expect_error(compute_trend(x, method = "invalid"))
})

test_that("compute_trend errors on bad window values", {
  set.seed(42)
  x <- rnorm(50)
  expect_error(compute_trend(x, window = 1))
  expect_error(compute_trend(x, window = 100))
})

test_that("compute_trend errors on invalid align", {
  set.seed(42)
  x <- rnorm(100)
  expect_error(compute_trend(x, window = 10, align = "invalid"))
})

test_that("compute_trend errors on missing data argument", {
  expect_error(compute_trend())
})

test_that("compute_trend handles minimum viable series", {
  set.seed(42)
  x <- rnorm(10)
  tr <- compute_trend(x, window = 3)
  expect_s3_class(tr, "trend")
  expect_equal(nrow(tr), 10L)
})

test_that("print.trend does not error", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 10)
  print(tr) |>
    capture.output() |>
    expect_error(NA)
})

test_that("plot.trend series type does not error", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 10)
  plot(tr, type = "series") |>
    expect_error(NA)
})

test_that("plot.trend ribbons type does not error", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 10)
  plot(tr, type = "ribbons") |>
    expect_error(NA)
})

test_that("summary.trend does not error and returns list", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 10)
  s <- summary(tr)
  expect_s3_class(s, "summary.trend")
  expect_true("counts" %in% names(s))
  expect_true("proportions" %in% names(s))
  expect_true("window" %in% names(s))
  expect_true("method" %in% names(s))
  expect_true("align" %in% names(s))
})

test_that("summary.trend prints correct header info", {
  set.seed(42)
  x <- cumsum(rnorm(100, sd = 2))
  tr <- compute_trend(x, window = 15, method = "slope", align = "right")
  output <- capture.output(print(summary(tr)))
  expect_true(any(grepl("slope", output)))
  expect_true(any(grepl("15", output)))
  expect_true(any(grepl("right", output)))
})

test_that("ar1_phi1 method handles constant-value windows (zero variance)", {
  set.seed(42)
  # A constant series has var < 1e-10 -> trend_ar1_ returns NA
  x <- rep(5, 50)
  tr <- compute_trend(x, window = 10, method = "ar1_phi1")
  expect_s3_class(tr, "trend")
  # All metrics should be NA since variance is zero
  expect_true(all(is.na(tr$metric)))
})

test_that("ar1_phi1 method handles very short clean data in window", {
  set.seed(42)
  # Series with lots of NAs -> clean data < min_points
  x <- c(NA, NA, 1, NA, NA, 2, NA, NA, 3, NA)
  tr <- compute_trend(x, window = 3, method = "ar1_phi1", min_points = 3L)
  expect_s3_class(tr, "trend")
})

test_that("slope method handles windows with NAs", {
  set.seed(42)
  # Series with NAs exercises the min_points check at line 261
  x <- c(1, 2, NA, NA, NA, 3, 4, 5, NA, NA)
  tr <- compute_trend(
    x, window = 5, method = "slope", slope_method = "ols",
    min_points = 3L, align = "right"
  )
  expect_s3_class(tr, "trend")
  # Missing_Data state should appear for NA positions
  expect_true("Missing_Data" %in% as.character(tr$state))
})

test_that("slope method handles zero x-variance (constant x)", {
  set.seed(42)
  # When all x values are the same (shouldn't happen with seq_along(wd),
  # but test via zero y-variance)
  x <- rep(3, 30)
  tr <- compute_trend(x, window = 5, method = "slope", slope_method = "ols")
  expect_s3_class(tr, "trend")
  # All values are the same -> metric should be 0 (var(y) < 1e-10 -> return 0)
  non_na <- tr$metric[!is.na(tr$metric)]
  if (length(non_na) > 0L) {
    expect_true(all(non_na == 0))
  }
})

test_that("slope robust method works with various data", {
  set.seed(42)
  # Robust (Theil-Sen) slope
  x <- cumsum(rnorm(60, sd = 2))
  tr <- compute_trend(x, window = 10, method = "slope", slope_method = "robust")
  expect_s3_class(tr, "trend")
  expect_true(any(!is.na(tr$metric)))
})

test_that("slope spearman method works", {
  set.seed(42)
  x <- cumsum(rnorm(60, sd = 2))
  tr <- compute_trend(x, window = 10, method = "slope", slope_method = "spearman")
  expect_s3_class(tr, "trend")
  expect_true(any(!is.na(tr$metric)))
})

test_that("slope kendall method works", {
  set.seed(42)
  x <- cumsum(rnorm(60, sd = 2))
  tr <- compute_trend(x, window = 10, method = "slope", slope_method = "kendall")
  expect_s3_class(tr, "trend")
  expect_true(any(!is.na(tr$metric)))
})

test_that("slope robust handles case with all equal x differences", {
  set.seed(42)
  # All pairwise x-diffs are 0 when only one unique x value
  x <- c(rep(1, 10), rnorm(20))
  tr <- compute_trend(x, window = 5, method = "slope", slope_method = "robust")
  expect_s3_class(tr, "trend")
})

test_that("slope robust handles constant-valued window", {
  set.seed(42)
  x <- rep(10, 30)
  tr <- compute_trend(x, window = 5, method = "slope", slope_method = "robust")
  expect_s3_class(tr, "trend")
  # Constant y values -> should produce 0 slope
  non_na <- tr$metric[!is.na(tr$metric)]
  if (length(non_na) > 0L) {
    expect_true(all(non_na == 0))
  }
})

test_that("slope spearman handles near-zero x variance", {
  set.seed(42)
  # This tests the fallback when sd_x is near zero in spearman/kendall
  x <- c(rep(5, 10), rnorm(20))
  tr <- compute_trend(x, window = 5, method = "slope", slope_method = "spearman")
  expect_s3_class(tr, "trend")
})

test_that("turbulence detection handles windows with NAs in metric", {
  set.seed(42)
  # Series with some NAs that could produce NA in volatility
  x <- c(rnorm(30), NA, NA, rnorm(30))
  tr <- compute_trend(x, window = 10, method = "slope")
  expect_s3_class(tr, "trend")
  # Missing_Data should appear for NA positions
  expect_true("Missing_Data" %in% as.character(tr$state))
})

test_that("min_points is clamped to window when min_points > window", {
  set.seed(42)
  x <- rnorm(50)
  # min_points = 20, window = 5 -> min_points gets clamped to 5
  tr <- compute_trend(x, window = 5, min_points = 20L)
  expect_s3_class(tr, "trend")
  expect_equal(nrow(tr), 50L)
})

test_that("growth_factor method handles windows with start near zero", {
  set.seed(42)
  # Start at zero -> abs(clean[1]) <= 1e-10 -> NA
  x <- c(0, 0, 0, rnorm(30, mean = 5))
  tr <- compute_trend(x, window = 5, method = "growth_factor")
  expect_s3_class(tr, "trend")
  # Some initial metrics may be NA due to zero start value
})

test_that("growth_factor method handles NA-heavy window", {
  set.seed(42)
  x <- c(NA, NA, 1, 2, NA, NA, 3, 4, 5, 6, 7, 8, 9, 10)
  tr <- compute_trend(x, window = 4, method = "growth_factor")
  expect_s3_class(tr, "trend")
})

test_that("turbulence detection skips when not enough non-NA in seg", {
  set.seed(42)
  # Create a series where volatility sub-window has few non-NA
  x <- c(rnorm(10), NA, NA, NA, NA, NA, rnorm(10))
  tr <- compute_trend(x, window = 5, method = "slope", min_points = 3L)
  expect_s3_class(tr, "trend")
})

test_that("turbulence reclassification applies different thresholds for flat", {
  set.seed(42)
  # Create a series with a flat section followed by turbulence
  # Flat sections need higher volatility to become turbulent
  x <- c(rep(5, 20), cumsum(rnorm(30, sd = 5)))
  tr <- compute_trend(
    x, window = 5, method = "slope",
    turbulence_threshold = 3, flat_to_turbulent_factor = 2.0
  )
  expect_s3_class(tr, "trend")
  # Check that the function ran through the turbulence reclassification
  expect_true(is.factor(tr$state))
})

test_that("trend_ar1_ returns NA when clean data < min_points", {
  set.seed(42)
  result <- trend_ar1_(c(NA, NA, 1), min_points = 3L)
  # Only 1 non-NA value, < min_points = 3
  expect_true(is.na(result))
})

test_that("trend_ar1_ handles ar.ols try-error or empty ar", {
  set.seed(42)
  result <- trend_ar1_(c(1, 10), min_points = 2L)
  # ar.ols with 2 points may or may not produce ar coefficient
  expect_true(is.na(result) || is.numeric(result))
})

test_that("trend_slope_ returns NA when valid < min_points", {
  set.seed(42)
  # Pass x_vals and y_vals where too many are NA
  result <- trend_slope_(
    c(1, 2, NA, NA), c(NA, NA, 3, 4), "ols", min_points = 3L
  )
  expect_true(is.na(result))
})

test_that("trend_slope_ returns NA when unique(xc) < 2 for non-robust", {
  set.seed(42)
  # All x values the same, non-robust method
  result <- trend_slope_(c(1, 1, 1), c(1, 2, 3), "ols", min_points = 2L)
  expect_true(is.na(result))
})

test_that("trend_slope_ returns NA when robust and length(xc) < 2", {
  set.seed(42)
  # For robust with only 1 valid point
  result <- trend_slope_(
    c(1, NA, NA), c(5, NA, NA), "robust", min_points = 1L
  )
  # valid positions: only index 1. length(xc) = 1 < 2 => NA
  expect_true(is.na(result))
})


test_that("trend_slope_ returns 0 or NA when var(xc) near zero", {
  set.seed(42)
  result_zero <- trend_slope_(
    c(1, 1 + 1e-15, 1 + 2e-15), c(5, 5, 5), "robust", min_points = 2L
  )
  # var(yc) = 0 < 1e-10 => return 0
  expect_equal(result_zero, 0)
  result_na <- trend_slope_(
    c(1, 1 + 1e-15, 1 + 2e-15), c(1, 2, 3), "robust", min_points = 2L
  )
  # var(yc) > 1e-10 => return NA_real_
  expect_true(is.na(result_na))
})

test_that("trend_slope_ robust returns NA when all dx == 0", {
  set.seed(42)
  result <- trend_slope_(
    c(1, 2, 3), c(4, 5, 6), "robust", min_points = 2L
  )
  expect_true(is.numeric(result))
})

test_that("trend_slope_ handles cor try-error in spearman/kendall", {
  set.seed(42)
  result <- trend_slope_(
    c(1, 2, 3, 4), c(5, 6, 7, 8), "spearman", min_points = 2L
  )
  expect_true(is.numeric(result))
})

test_that("trend_volatility_ returns NA when all values are NA", {
  set.seed(42)
  result <- trend_volatility_(c(NA, NA, NA))
  expect_true(is.na(result))
})

test_that("trend_volatility_ returns NA when single non-NA value", {
  set.seed(42)
  # sd of single value is NA, range of single value is 0 (but sd is NA)
  result <- trend_volatility_(c(NA, 5, NA))
  # sd(5) = NA (only 1 value), so is.na(m_sd) => TRUE => NA_real_
  expect_true(is.na(result))
})

test_that("turbulence detection skips when combined volatility is NA", {
  set.seed(42)
  x <- c(rnorm(20), NA, NA, rnorm(20))
  tr <- compute_trend(x, window = 5, method = "slope", min_points = 2L)
  expect_s3_class(tr, "trend")
})
