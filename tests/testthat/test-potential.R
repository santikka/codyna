set.seed(33)
bimodal <- c(rnorm(100, -2), rnorm(100, 2))

test_that("potential_analysis returns class 'potential'", {
  pa <- potential_analysis(bimodal)
  expect_s3_class(pa, "potential")
  expect_true(inherits(pa, "list"))
})

test_that("result contains required elements", {
  pa <- potential_analysis(bimodal)
  expected <- c("landscape", "wells", "barriers", "n_wells", "values", "time")
  vapply(expected, function(nm) {
    expect_true(nm %in% names(pa), info = paste("missing element:", nm))
    TRUE
  }, logical(1))
})

test_that("landscape is a tibble with correct columns", {
  pa <- potential_analysis(bimodal)
  expect_s3_class(pa$landscape, "tbl_df")
  expect_true(all(c("x", "density", "potential") %in% names(pa$landscape)))
})

test_that("wells is a tibble with correct columns", {
  pa <- potential_analysis(bimodal)
  expect_s3_class(pa$wells, "tbl_df")
  expect_true(all(c("location", "depth", "width") %in% names(pa$wells)))
})

test_that("barriers is a tibble with correct columns", {
  pa <- potential_analysis(bimodal)
  expect_s3_class(pa$barriers, "tbl_df")
  expect_true(all(c("location", "height") %in% names(pa$barriers)))
})

test_that("n_wells is a positive integer", {
  pa <- potential_analysis(bimodal)
  expect_true(is.numeric(pa$n_wells))
  expect_true(pa$n_wells >= 1L)
  expect_equal(pa$n_wells, as.integer(pa$n_wells))
})

test_that("detrend 'none' works", {
  potential_analysis(bimodal, detrend = "none") |>
    expect_error(NA)
})

test_that("detrend 'linear' works", {
  potential_analysis(bimodal, detrend = "linear") |>
    expect_error(NA)
})

test_that("detrend 'diff' works and shortens series by one", {
  pa <- potential_analysis(bimodal, detrend = "diff")
  expect_error(pa, NA)
  expect_equal(length(pa$values), length(bimodal) - 1L)
})

test_that("rolling window adds 'rolling' element", {
  pa <- potential_analysis(bimodal, window = 50L)
  expect_true("rolling" %in% names(pa))
  expect_s3_class(pa$rolling, "tbl_df")
  expect_true(
    all(c("time", "value", "n_wells", "dominant_well_location") %in%
          names(pa$rolling))
  )
})

test_that("global analysis does not have 'rolling' element", {
  pa <- potential_analysis(bimodal)
  expect_null(pa$rolling)
})

test_that("n_bins affects landscape grid resolution", {
  pa_small <- potential_analysis(bimodal, n_bins = 20L)
  pa_large <- potential_analysis(bimodal, n_bins = 100L)
  expect_true(nrow(pa_small$landscape) < nrow(pa_large$landscape))
})

test_that("print.potential does not error", {
  pa <- potential_analysis(bimodal)
  print(pa) |>
    capture.output() |>
    expect_error(NA)
})

test_that("plot.potential does not error", {
  pa <- potential_analysis(bimodal)
  plot(pa) |>
    expect_error(NA)
  plot(pa, type = "density") |>
    expect_error(NA)
  plot(pa, type = "both") |>
    expect_error(NA)
})

test_that("plot.potential with rolling does not error", {
  pa <- potential_analysis(bimodal, window = 50L)
  plot(pa) |>
    expect_error(NA)
})

test_that("bad detrend produces error", {
  expect_error(potential_analysis(bimodal, detrend = "quadratic"))
})

test_that("bad n_bins produces error", {
  expect_error(potential_analysis(bimodal, n_bins = 3L))
  expect_error(potential_analysis(bimodal, n_bins = 20000L))
})

test_that("potential values are numeric", {
  pa <- potential_analysis(bimodal)
  expect_true(is.numeric(pa$landscape$potential))
  expect_false(any(is.na(pa$landscape$potential)))
})

test_that("bimodal data detects at least 2 wells", {
  set.seed(42)
  clear_bimodal <- c(rnorm(200, -3, 0.5), rnorm(200, 3, 0.5))
  pa <- potential_analysis(clear_bimodal)
  expect_true(pa$n_wells >= 2L)
})

test_that("wells have non-negative depth", {
  pa <- potential_analysis(bimodal)
  expect_true(all(pa$wells$depth >= 0))
})

test_that("barriers have positive height when present", {
  set.seed(42)
  clear_bimodal <- c(rnorm(200, -3, 0.5), rnorm(200, 3, 0.5))
  pa <- potential_analysis(clear_bimodal)
  if (nrow(pa$barriers) > 0L) {
    expect_true(all(pa$barriers$height > 0))
  }
})

test_that("attributes store analysis settings", {
  pa <- potential_analysis(bimodal, n_bins = 30L, detrend = "linear")
  expect_equal(attr(pa, "n_bins"), 30L)
  expect_equal(attr(pa, "detrend"), "linear")
  expect_null(attr(pa, "window"))
})

test_that("potential_find_minima_ returns empty for short vector", {
  # n < 3 returns integer(0)
  expect_length(potential_find_minima_(c(1, 2)), 0L)
  expect_length(potential_find_minima_(numeric(0)), 0L)
  expect_length(potential_find_minima_(5), 0L)
})

test_that("potential_find_minima_ handles zero sign forward-filling", {
  # sign_d[i] == 0 && i > 1 => forward-fill
  # Create a sequence with a flat portion: 3, 2, 1, 1, 1, 2, 3
  # diff = -1, -1, 0, 0, 1, 1 => sign = -1, -1, 0, 0, 1, 1
  # After forward-fill: -1, -1, -1, -1, 1, 1
  # sign change at index 4->5: minima at position 5
  x <- c(3, 2, 1, 1, 1, 2, 3)
  result <- potential_find_minima_(x)
  expect_true(length(result) >= 1L)
})

test_that("potential_find_maxima_ returns empty for short vector", {
  # n < 3 returns integer(0)
  expect_length(potential_find_maxima_(c(1, 2)), 0L)
  expect_length(potential_find_maxima_(numeric(0)), 0L)
})

test_that("potential_find_maxima_ handles zero sign forward-filling", {
  # sign_d[i] == 0 && i > 1 => forward-fill
  # Create: 1, 2, 3, 3, 3, 2, 1
  # diff = 1, 1, 0, 0, -1, -1 => sign = 1, 1, 0, 0, -1, -1
  # After forward-fill: 1, 1, 1, 1, -1, -1
  # sign change at index 4->5: maxima at position 5
  x <- c(1, 2, 3, 3, 3, 2, 1)
  result <- potential_find_maxima_(x)
  expect_true(length(result) >= 1L)
})

test_that("potential_landscape_ respects bandwidth parameter", {
  set.seed(33)
  x <- rnorm(200)
  # bandwidth is not NULL
  result <- potential_landscape_(x, n_bins = 50L, bandwidth = 0.5)
  expect_true(is.list(result))
  expect_true("landscape" %in% names(result))
  expect_true("wells" %in% names(result))
})

test_that("potential_landscape_ uses global min when no interior minimum", {
  # no interior minimum found, fallback to global minimum
  # Unimodal data with very tight distribution = single potential well
  set.seed(42)
  x <- rnorm(500, 0, 0.3)
  result <- potential_landscape_(x, n_bins = 20L, bandwidth = NULL)
  expect_true(result$n_wells >= 1L)
})

test_that("potential_analysis accepts bandwidth parameter", {
  set.seed(42)
  pa <- potential_analysis(bimodal, bandwidth = 0.5)
  expect_equal(attr(pa, "bandwidth"), 0.5)
  expect_s3_class(pa, "potential")
})

test_that("potential_analysis rejects negative bandwidth", {
  # check_range on bandwidth
  expect_error(potential_analysis(bimodal, bandwidth = -1))
})

test_that("rolling window produces non-NA n_wells for valid windows", {
  set.seed(42)
  x <- rnorm(150)
  # windows with too few clean observations get NA
  # Use a small window so first positions are NA (pre-fill)
  pa <- potential_analysis(x, window = 40L)
  expect_true("rolling" %in% names(pa))
  # First (window-1) positions should be NA_integer_
  expect_true(all(is.na(pa$rolling$n_wells[1:39])))
  # Later positions should have valid n_wells
  expect_true(any(!is.na(pa$rolling$n_wells[40:150])))
})

test_that("rolling window with detrend diff works", {
  set.seed(42)
  x <- cumsum(rnorm(200))
  pa <- potential_analysis(x, window = 50L, detrend = "diff")
  expect_true("rolling" %in% names(pa))
  expect_s3_class(pa$rolling, "tbl_df")
})

test_that("rolling window dominant well location is tracked", {
  set.seed(42)
  # Bimodal data should show wells moving in rolling analysis
  x <- c(rnorm(100, -3, 0.5), rnorm(100, 3, 0.5))
  pa <- potential_analysis(x, window = 50L)
  rolling <- pa$rolling
  # dominant_well_location should have non-NA values for valid windows
  valid_dom <- rolling$dominant_well_location[!is.na(rolling$n_wells)]
  expect_true(any(!is.na(valid_dom)))
})

test_that("plot.potential landscape type with rolling appends panel", {
  set.seed(42)
  pa <- potential_analysis(bimodal, window = 50L)
  # landscape type with rolling present does NOT early-return
  p <- plot(pa, type = "landscape")
  # Should produce a patchwork (landscape + rolling panels)
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
  # Force rendering to execute lazy closures (e.g., breaks function)
  grDevices::pdf(tempfile())
  suppressWarnings(print(p))
  grDevices::dev.off()
})

test_that("plot.potential density type returns single ggplot", {
  set.seed(42)
  pa <- potential_analysis(bimodal)
  p <- plot(pa, type = "density")
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

test_that("plot.potential both type with rolling includes all panels", {
  set.seed(42)
  pa <- potential_analysis(bimodal, window = 50L)
  p <- plot(pa, type = "both")
  expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
})

test_that("print.potential shows rolling info when present", {
  set.seed(42)
  pa <- potential_analysis(bimodal, window = 50L)
  # rolling section in print output
  out <- capture.output(print(pa))
  out_str <- paste(out, collapse = "\n")
  expect_true(
    grepl("Rolling", out_str) || grepl("rolling", out_str) ||
      grepl("n_wells range", out_str)
  )
})

test_that("print.potential shows bandwidth info", {
  set.seed(42)
  pa <- potential_analysis(bimodal, bandwidth = 0.5)
  out <- capture.output(print(pa))
  out_str <- paste(out, collapse = "\n")
  expect_true(grepl("0.5", out_str))
})

test_that("print.potential shows auto bandwidth when NULL", {
  set.seed(42)
  pa <- potential_analysis(bimodal)
  out <- capture.output(print(pa))
  out_str <- paste(out, collapse = "\n")
  expect_true(grepl("Silverman", out_str))
})

test_that("print.potential shows barrier details when present", {
  set.seed(42)
  clear_bimodal <- c(rnorm(200, -3, 0.5), rnorm(200, 3, 0.5))
  pa <- potential_analysis(clear_bimodal)
  out <- capture.output(s <- print(pa))
  out_str <- paste(out, collapse = "\n")
  if (nrow(pa$barriers) > 0L) {
    expect_true(grepl("Barriers", out_str))
  }
})

test_that("potential_analysis accepts ts object", {
  set.seed(42)
  x <- ts(rnorm(100), start = 1, frequency = 1)
  pa <- potential_analysis(x)
  expect_s3_class(pa, "potential")
})

test_that("potential_landscape_ falls back to global min when no interior minimum", {
  # when potential_find_minima_ returns integer(0),
  # the code falls back to which.min(potential_vals)
  # KDE always produces a smooth bell shape with at least 1 interior minimum,
  # so we mock potential_find_minima_ to return empty to exercise this branch
  local_mocked_bindings(
    potential_find_minima_ = function(x) integer(0),
    .package = "codyna"
  )
  set.seed(33)
  x <- rnorm(100)
  result <- potential_landscape_(x, n_bins = 50L, bandwidth = NULL)
  # Should find exactly 1 well via global minimum fallback
  expect_equal(result$n_wells, 1L)
  expect_true(nrow(result$wells) == 1L)
})

#  w_clean < 20 in rolling loop
# These lines require window data with NAs so that fewer than 20 clean values
# remain. The exported function's global landscape call (line 402) also gets
# NAs, so we mock potential_landscape_ to handle NAs gracefully.

test_that("rolling window marks NA when w_clean < 20 via mocked landscape", {
  original_fn <- potential_landscape_
  local_mocked_bindings(
    potential_landscape_ = function(x, n_bins, bandwidth) {
      x_clean <- x[!is.na(x)]
      original_fn(x_clean, n_bins, bandwidth)
    },
    .package = "codyna"
  )
  set.seed(42)
  # First 25 values are NA, then 80 real values
  x <- c(rep(NA_real_, 25), rnorm(80))
  # Window = 30: windows covering the NA section have w_clean < 20
  pa <- potential_analysis(x, window = 30L)
  expect_s3_class(pa, "potential")
  rolling <- pa$rolling
  # Early windows covering the NA section should have NA n_wells
  expect_true(any(is.na(rolling$n_wells)))
})

test_that("rolling window marks NA when potential_landscape_ errors", {
  original_fn <- codyna:::potential_landscape_
  call_count <- 0L
  local_mocked_bindings(
    potential_landscape_ = function(x, n_bins, bandwidth) {
      call_count <<- call_count + 1L
      # Fail the first 3 rolling-window calls
      if (call_count <= 3L) stop("simulated landscape error")
      original_fn(x, n_bins, bandwidth)
    },
    .package = "codyna"
  )
  set.seed(33)
  x <- rnorm(100)
  pa <- potential_analysis(x, window = 30L)
  expect_s3_class(pa, "potential")
  rolling <- pa$rolling
  # First 3 rolling-window calls errored -> those positions get NA
  expect_true(any(is.na(rolling$n_wells)))
})
