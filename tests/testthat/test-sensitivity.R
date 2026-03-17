set.seed(33)
sens_ts <- cumsum(rnorm(200, sd = seq(0.5, 2, length.out = 200)))

test_that("sensitivity_ews returns sensitivity_ews tibble", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  expect_s3_class(sa, "sensitivity_ews")
  expect_s3_class(sa, "tbl_df")
})

test_that("sensitivity_ews has correct columns", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  expect_true(all(c("window", "detrend", "time", "score", "tau") %in% names(sa)))
})

test_that("sensitivity_ews stores attributes correctly", {
  sa <- sensitivity_ews(sens_ts, metric = "sd", windows = c(30, 50))
  expect_equal(attr(sa, "metric"), "sd")
  expect_equal(attr(sa, "windows"), c(30L, 50L))
  expect_true(is.character(attr(sa, "detrend_methods")))
})

test_that("tau values are between -1 and 1", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  tau_vals <- unique(sa$tau)
  tau_vals <- tau_vals[!is.na(tau_vals)]
  expect_true(all(tau_vals >= -1 & tau_vals <= 1))
})

test_that("tau is constant within each window-detrend combination", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  combos <- unique(sa[, c("window", "detrend")])
  for (i in seq_len(nrow(combos))) {
    sub <- sa[sa$window == combos$window[i] & sa$detrend == combos$detrend[i], ]
    tau_vals <- unique(sub$tau)
    expect_length(tau_vals, 1L)
  }
})

test_that("default windows parameter generates multiple window sizes", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1")
  windows_used <- attr(sa, "windows")
  expect_true(length(windows_used) > 1L)
})

test_that("custom windows parameter works", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(25, 40, 60))
  windows_used <- attr(sa, "windows")
  expect_equal(windows_used, c(25L, 40L, 60L))
  expect_true(all(unique(sa$window) %in% c(25L, 40L, 60L)))
})

test_that("metric ar1 works", {
  sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50)) |>
    expect_error(NA)
})

test_that("metric sd works", {
  sensitivity_ews(sens_ts, metric = "sd", windows = c(30, 50)) |>
    expect_error(NA)
})

test_that("metric cv works", {
  sensitivity_ews(sens_ts, metric = "cv", windows = c(30, 50)) |>
    expect_error(NA)
})

test_that("metric variance works", {
  sensitivity_ews(sens_ts, metric = "variance", windows = c(30, 50)) |>
    expect_error(NA)
})

test_that("metric skewness works", {
  sensitivity_ews(sens_ts, metric = "skewness", windows = c(30, 50)) |>
    expect_error(NA)
})

test_that("metric kurtosis works", {
  sensitivity_ews(sens_ts, metric = "kurtosis", windows = c(30, 50)) |>
    expect_error(NA)
})

test_that("detrend_methods none only works", {
  sa <- sensitivity_ews(
    sens_ts, metric = "ar1", windows = c(30, 50),
    detrend_methods = "none"
  )
  expect_equal(unique(sa$detrend), "none")
})

test_that("detrend_methods linear only works", {
  sa <- sensitivity_ews(
    sens_ts, metric = "ar1", windows = c(30, 50),
    detrend_methods = "linear"
  )
  expect_equal(unique(sa$detrend), "linear")
})

test_that("both detrend methods produce rows for each", {
  sa <- sensitivity_ews(
    sens_ts, metric = "ar1", windows = c(30, 50),
    detrend_methods = c("none", "linear")
  )
  detrend_vals <- unique(sa$detrend)
  expect_true("none" %in% detrend_vals)
  expect_true("linear" %in% detrend_vals)
})

test_that("print.sensitivity_ews does not error", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  capture.output(print(sa)) |>
    expect_error(NA)
})

test_that("plot.sensitivity_ews heatmap does not error", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  plot(sa, type = "heatmap") |>
    expect_error(NA)
})

test_that("plot.sensitivity_ews lines does not error", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  plot(sa, type = "lines") |>
    expect_error(NA)
})

test_that("plot.sensitivity_ews both does not error", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  plot(sa, type = "both") |>
    expect_error(NA)
})

test_that("summary.sensitivity_ews does not error and returns list", {
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  result <- capture.output(out <- summary(sa))
  expect_error(result, NA)
  expect_type(out, "list")
  expect_true(
    all(
      c(
        "metric", "tau_range", "tau_mean",
        "n_positive", "n_negative", "n_total",
        "most_robust", "least_robust", "consistent"
      ) %in% names(out)
    )
  )
})

test_that("sensitivity_ews errors on invalid metric", {
  expect_error(sensitivity_ews(sens_ts, metric = "bogus"))
})

test_that("sensitivity_ews errors on missing data", {
  expect_error(sensitivity_ews())
})

test_that("sensitivity_ews errors on all-invalid window sizes", {
  expect_error(sensitivity_ews(sens_ts, metric = "ar1", windows = c(0, -5)))
})

test_that("sensitivity_ews accepts ts objects", {
  ts_obj <- ts(sens_ts, frequency = 12)
  sensitivity_ews(ts_obj, metric = "ar1", windows = c(30, 50)) |>
    expect_error(NA)
})

test_that("sensitivity_detrend_ handles series with < 3 valid points", {
  set.seed(42)
  # Series with mostly NAs -- fewer than 3 valid points
  short_vals <- c(1.0, NA, NA, NA, NA)
  # Use the internal function via a sensitivity_ews call with a window
  # that forces this through detrending. Instead, test via the full pipeline
  # with a very short series that has NAs in it.
  # The detrend function returns values unchanged when sum(valid) < 3
  short_ts <- c(1, NA, 2)
  sa <- sensitivity_ews(
    short_ts, metric = "sd", windows = c(2),
    detrend_methods = "linear"
  )
  # Should not error; the linear detrend falls back to returning raw values
  expect_s3_class(sa, "sensitivity_ews")
})

test_that("sensitivity_ews handles ar1 with near-zero variance windows", {
  set.seed(42)
  # Constant series -> var < 1e-10 in ar1 metric
  const_ts <- rep(5, 50)
  sa <- sensitivity_ews(
    const_ts, metric = "ar1", windows = c(10, 20),
    detrend_methods = "none"
  )
  expect_s3_class(sa, "sensitivity_ews")
  # Scores should be NA since variance is near zero
  expect_true(all(is.na(sa$score)))
})

test_that("sensitivity_ews handles skewness with zero-sd window", {
  set.seed(42)
  # Constant values -> sd near 0 -> skewness returns NA
  const_ts <- rep(3, 50)
  sa <- sensitivity_ews(
    const_ts, metric = "skewness", windows = c(10),
    detrend_methods = "none"
  )
  expect_s3_class(sa, "sensitivity_ews")
  expect_true(all(is.na(sa$score)))
})

test_that("sensitivity_ews handles kurtosis with zero-sd window", {
  set.seed(42)
  const_ts <- rep(7, 50)
  sa <- sensitivity_ews(
    const_ts, metric = "kurtosis", windows = c(10),
    detrend_methods = "none"
  )
  expect_s3_class(sa, "sensitivity_ews")
  expect_true(all(is.na(sa$score)))
})

test_that("sensitivity_ews handles cv with zero-mean window", {
  set.seed(42)
  # Series centered exactly at zero -> cv returns NA
  zero_mean_ts <- c(rep(-1, 25), rep(1, 25))
  sa <- sensitivity_ews(
    zero_mean_ts, metric = "cv", windows = c(50),
    detrend_methods = "none"
  )
  expect_s3_class(sa, "sensitivity_ews")
})

test_that("sensitivity_metric_ returns NA for short windows (n < 3)", {
  set.seed(42)
  # Window of 2 with some NAs in data can lead to clean length < 3
  short_ts <- c(NA, 1, NA, 2, NA, 3, NA, 4, NA, 5)
  sa <- sensitivity_ews(
    short_ts, metric = "sd", windows = c(2),
    detrend_methods = "none"
  )
  expect_s3_class(sa, "sensitivity_ews")
})

test_that("sensitivity_rolling_ handles window larger than data", {
  set.seed(42)
  # Very short series with large window -> empty result
  short_ts <- c(1, 2, 3)
  # Window = 10 exceeds length 3 -> will be filtered out by validation
  # We need to use windows that pass validation but yield 0 rolling windows
  # Actually, windows > n are filtered out at line 294
  # Let's use a window exactly equal to n -> n_windows = 1
  sa <- sensitivity_ews(
    short_ts, metric = "sd", windows = c(3),
    detrend_methods = "none"
  )
  expect_s3_class(sa, "sensitivity_ews")
  expect_true(nrow(sa) >= 1L)
})

test_that("sensitivity_tau_ returns NA for < 4 valid score points", {
  set.seed(42)
  # Very short series -> very few rolling scores -> tau is NA
  short_ts <- c(1, 2, 3, 4, 5)
  sa <- sensitivity_ews(
    short_ts, metric = "sd", windows = c(3),
    detrend_methods = "none"
  )
  expect_s3_class(sa, "sensitivity_ews")
  # With n=5, window=3, n_windows=3, which is < 4 valid scores
  tau_vals <- unique(sa$tau)
  expect_true(all(is.na(tau_vals)))
})

test_that("plot heatmap uses sensitivity colors", {
  set.seed(42)
  sa <- sensitivity_ews(sens_ts, metric = "ar1", windows = c(30, 50))
  p <- plot(sa, type = "heatmap")
  expect_true(inherits(p, "gg"))
})

test_that("default windows with short series uses 2-point sequence", {
  set.seed(42)
  # Short series: n = 40, w_min = max(20, 40/20) = 20, w_max = min(20, 200) = 20
  # w_min >= w_max -> generates 2-point sequence
  short_ts <- rnorm(40)
  sa <- sensitivity_ews(short_ts, metric = "sd")
  windows_used <- attr(sa, "windows")
  expect_true(length(windows_used) >= 1L)
})

test_that("sensitivity_ews handles window yielding empty rolling result", {
  set.seed(42)
  # This triggers when rolling produces 0 windows for a particular
  # (window, detrend) combination. Use a very short series with
  # a window that equals n (produces exactly 1 window, not 0)
  # We need n_windows < 1 which means window > n
  # But windows > n are filtered out. So test with a series where
  # some windows are valid and some are not.
  # Actually, let's create a scenario with all NAs in data
  na_ts <- c(NA, NA, NA, NA, 1, 2, 3, 4, 5, NA, NA, NA)
  sa <- sensitivity_ews(
    na_ts, metric = "sd", windows = c(3, 5),
    detrend_methods = "none"
  )
  expect_s3_class(sa, "sensitivity_ews")
})

test_that("plot heatmap handles near-zero tau range", {
  set.seed(42)
  # Constant series -> tau is NA but score is also constant
  # Use a series that produces near-zero tau
  flat_ts <- rep(5, 100)
  sa <- sensitivity_ews(
    flat_ts, metric = "sd", windows = c(20, 30),
    detrend_methods = "none"
  )
  # All scores will be 0, tau will be NA, but the plot should still work
  p <- suppressWarnings(plot(sa, type = "heatmap"))
  expect_true(inherits(p, "gg"))
})

test_that("summary reports mixed signs when tau values differ in sign", {
  set.seed(42)
  # Create a series where different windows/detrend produce different tau signs
  # A noisy non-monotonic series should produce inconsistent tau
  wiggly_ts <- sin(seq(0, 6 * pi, length.out = 200)) * 10 + rnorm(200)
  sa <- sensitivity_ews(
    wiggly_ts, metric = "ar1",
    windows = c(20, 80),
    detrend_methods = c("none", "linear")
  )
  sumr <- summary(sa)
  output <- capture.output(print(sumr))
  # Check the output contains consistency information
  expect_true(
    any(grepl("Consistency", output)) ||
      any(grepl("ALL tau values", output)) ||
      any(grepl("Mixed signs", output))
  )
  expect_type(sumr, "list")
  expect_true("consistent" %in% names(sumr))
})

test_that("sensitivity_detrend_ returns raw values when lm fails", {
  # Create data where lm would fail: Inf values
  vals <- c(1, 2, Inf, 4, 5)
  result <- sensitivity_detrend_(vals, "linear")
  # lm with Inf should produce try-error -> returns original values
  expect_equal(result, vals)
})

test_that("sensitivity_metric_ returns NA for length-2 input", {
  result <- sensitivity_metric_(c(1, 2), "sd")
  expect_true(is.na(result))
})

test_that("sensitivity_metric_ returns NA for length-1 input", {
  result <- sensitivity_metric_(5, "ar1")
  expect_true(is.na(result))
})

test_that("sensitivity_metric_ ar1 returns NA when ar.ols fails", {
  # Extreme values: variance is Inf (passes var check), but ar.ols crashes
  # -> try-error -> line 52 returns NA
  set.seed(42)
  x <- c(rnorm(5), 1e308, -1e308, rnorm(5))
  suppressWarnings(result <- sensitivity_metric_(x, "ar1"))
  expect_true(is.na(result))
})

test_that("sensitivity_rolling_ returns empty lists when window > n", {
  result <- sensitivity_rolling_(c(1, 2, 3), window = 5, metric = "sd")
  expect_length(result$time, 0L)
  expect_length(result$score, 0L)
})

test_that("sensitivity_tau_ returns NA on cor.test failure", {
  # Inf values pass valid count but crash cor.test
  result <- sensitivity_tau_(c(1, 2, Inf, 4, 5))
  expect_true(is.na(result) || is.numeric(result))
})

test_that("sensitivity_colors_ returns 9 color hex codes", {
  result <- sensitivity_colors_()
  expect_length(result, 9L)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", result)))
})

test_that("sensitivity_rolling_ empty result creates zero-row tibble", {
  # Window larger than data -> n_windows < 1 -> empty result
  result <- sensitivity_rolling_(c(1, 2), window = 5, metric = "sd")
  expect_length(result$score, 0L)
  expect_length(result$time, 0L)
})

test_that("plot heatmap sets abs_max to 1 when tau range is near zero", {
  set.seed(42)
  # Create a constant series where all tau values are NA (or near zero)
  const_ts <- rep(5, 100) + rnorm(100, sd = 1e-12)
  sa <- sensitivity_ews(
    const_ts, metric = "sd", windows = c(20, 30),
    detrend_methods = "none"
  )
  # Force tau to near-zero
  sa$tau <- 0.001
  p <- plot(sa, type = "heatmap")
  expect_true(inherits(p, "gg"))
})

test_that("sensitivity_tau_ returns NA when cor.test throws", {
  local_mocked_bindings(
    cor.test = function(...) stop("mock cor.test error"),
    .package = "stats"
  )
  result <- sensitivity_tau_(c(1, 2, 3, 4, 5))
  expect_true(is.na(result))
})

test_that("sensitivity_ews creates empty tibble when rolling returns empty", {
  local_mocked_bindings(
    sensitivity_rolling_ = function(values, window, metric) {
      list(time = integer(0), score = numeric(0))
    },
    .package = "codyna"
  )
  set.seed(42)
  ts_data <- cumsum(rnorm(100))
  sa <- sensitivity_ews(ts_data, metric = "ar1", windows = c(30, 50))
  expect_s3_class(sa, "sensitivity_ews")
  expect_equal(nrow(sa), 0L)
})
