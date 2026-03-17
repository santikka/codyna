set.seed(11)
tna_ts <- cumsum(rnorm(200))

test_that("to_tna is a generic function", {
  expect_true(is.function(to_tna))
  expect_true(isGeneric("to_tna") || length(methods("to_tna")) > 0L)
})

test_that("default method errors on unsupported class", {
  expect_error(to_tna(data.frame(x = 1:10)))
  expect_error(to_tna(1:10))
  expect_error(to_tna(list(a = 1)))
})

test_that("to_tna.trend returns wide-format tibble", {
  tr <- compute_trend(tna_ts, window = 20L)
  result <- to_tna(tr)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_true(all(grepl("^T[0-9]+$", names(result))))
})

test_that("to_tna.trend values are character state labels", {
  tr <- compute_trend(tna_ts, window = 20L)
  result <- to_tna(tr)
  valid_states <- c("ascending", "descending", "flat", "turbulent",
                    "Missing_Data", "Initial")
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(all(vals %in% valid_states))
})

test_that("to_tna.trend column count matches non-NA states", {
  tr <- compute_trend(tna_ts, window = 20L)
  result <- to_tna(tr)
  n_non_na <- sum(!is.na(tr$state))
  expect_equal(ncol(result), n_non_na)
})

test_that("to_tna.hurst returns wide-format tibble", {
  h <- hurst(tna_ts, window = 30L, method = "dfa")
  result <- to_tna(h)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_true(all(grepl("^T[0-9]+$", names(result))))
})

test_that("to_tna.hurst values are valid hurst state labels", {
  h <- hurst(tna_ts, window = 30L, method = "dfa")
  result <- to_tna(h)
  valid_states <- c(
    "strong_antipersistent", "antipersistent", "random_walk",
    "persistent", "strong_persistent"
  )
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(all(vals %in% valid_states))
})

test_that("to_tna.hurst_global errors with informative message", {
  hg <- hurst(tna_ts, states = FALSE)
  expect_error(to_tna(hg), "rolling")
})

test_that("to_tna.spectral returns wide-format tibble", {
  sp <- spectral_ews(tna_ts, window = 30L, states = TRUE)
  result <- to_tna(sp)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_true(all(grepl("^T[0-9]+$", names(result))))
})

test_that("to_tna.spectral values are valid noise-colour states", {
  sp <- spectral_ews(tna_ts, window = 30L, states = TRUE)
  result <- to_tna(sp)
  valid_states <- c("white_noise", "pink_noise", "red_noise", "brownian")
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(all(vals %in% valid_states))
})

test_that("to_tna.spectral errors when states = FALSE", {
  sp <- spectral_ews(tna_ts, window = 30L, states = FALSE)
  expect_error(to_tna(sp), "state")
})

test_that("to_tna.changepoint returns wide-format tibble", {
  cp <- detect_cpts(tna_ts, method = "pelt")
  result <- to_tna(cp)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_true(all(grepl("^T[0-9]+$", names(result))))
})

test_that("to_tna.changepoint values are character", {
  cp <- detect_cpts(tna_ts, method = "pelt")
  result <- to_tna(cp)
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(is.character(vals))
})

test_that("to_tna.changepoint column count matches series length", {
  cp <- detect_cpts(tna_ts, method = "pelt")
  result <- to_tna(cp)
  # state column should have no NAs, so ncol = nrow(cp)
  n_non_na <- sum(!is.na(cp$state))
  expect_equal(ncol(result), n_non_na)
})

test_that("to_tna.changepoint state parameter selects column", {
  cp <- detect_cpts(tna_ts, method = "pelt")
  result_state <- to_tna(cp, state = "state")
  result_level <- to_tna(cp, state = "level")
  expect_s3_class(result_state, "tbl_df")
  expect_s3_class(result_level, "tbl_df")
  # Values from level column should be low/medium/high
  level_vals <- unlist(result_level[1, ], use.names = FALSE)
  expect_true(all(level_vals %in% c("low", "medium", "high")))
})

test_that("to_tna.multi_ews expanding returns wide-format tibble", {
  mews_df <- data.frame(
    Time = seq_len(100),
    V1 = cumsum(rnorm(100)),
    V2 = cumsum(rnorm(100)),
    V3 = cumsum(rnorm(100))
  )
  ews <- detect_multivariate_warnings(
    mews_df, method = "expanding", window = 50, time_col = "Time"
  )
  result <- to_tna(ews)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_true(all(grepl("^T[0-9]+$", names(result))))
})

test_that("to_tna.multi_ews rolling errors", {
  mews_df <- data.frame(
    Time = seq_len(100),
    V1 = cumsum(rnorm(100)),
    V2 = cumsum(rnorm(100)),
    V3 = cumsum(rnorm(100))
  )
  ews <- detect_multivariate_warnings(
    mews_df, method = "rolling", window = 50, time_col = "Time"
  )
  expect_error(to_tna(ews), "expanding")
})

test_that("to_tna.surrogate_test errors with informative message", {
  # Create a mock surrogate_test object
  mock_st <- structure(
    list(observed_tau = 0.5, p_value = 0.01),
    class = "surrogate_test"
  )
  expect_error(to_tna(mock_st), "surrogate_test")
})

test_that("to_tna.sensitivity_ews errors with informative message", {
  sa <- sensitivity_ews(tna_ts, metric = "ar1", windows = c(30, 50))
  expect_error(to_tna(sa), "sensitivity_ews")
})

test_that("to_tna output is always single-row tibble with T-columns", {
  tr <- compute_trend(tna_ts, window = 20L)
  h <- hurst(tna_ts, window = 30L)
  sp <- spectral_ews(tna_ts, window = 30L, states = TRUE)
  cpt <- detect_cpts(tna_ts, method = "pelt")
  results <- list(
    to_tna(tr),
    to_tna(h),
    to_tna(sp),
    to_tna(cpt)
  )
  for (r in results) {
    expect_s3_class(r, "tbl_df")
    expect_equal(nrow(r), 1L)
    expect_true(all(grepl("^T[0-9]+$", names(r))))
    # All values must be character
    vals <- unlist(r[1, ], use.names = FALSE)
    expect_true(is.character(vals))
  }
})

test_that("seq_to_wide_ errors when all states are NA", {
  # Create a mock trend object with only NA states
  mock_trend <- tibble::tibble(
    time = 1:5,
    value = rnorm(5),
    metric = rep(NA_real_, 5),
    state = rep(NA_character_, 5)
  )
  class(mock_trend) <- c("trend", "tbl_df", "tbl", "data.frame")
  expect_error(to_tna(mock_trend), "only NA")
})

test_that("to_tna.resilience returns wide-format tibble from classified object", {
  set.seed(42)
  res <- resilience(tna_ts, window = 30L)
  cls <- classify_resilience(res)
  result <- to_tna(cls)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_true(all(grepl("^T[0-9]+$", names(result))))
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(is.character(vals))
})

test_that("to_tna.resilience extracts per-metric state with state argument", {
  set.seed(42)
  res <- resilience(tna_ts, window = 30L)
  cls <- classify_resilience(res)
  # Extract a specific metric category (vsi)
  result_vsi <- to_tna(cls, state = "vsi")
  expect_s3_class(result_vsi, "tbl_df")
  expect_equal(nrow(result_vsi), 1L)
  vals <- unlist(result_vsi[1, ], use.names = FALSE)
  valid_cats <- c(
    "Excellent", "Solid", "Fair", "Vulnerable", "Failing", "Troubled"
  )
  expect_true(all(vals %in% valid_cats))
})

test_that("to_tna.resilience errors on unclassified resilience object", {
  set.seed(42)
  res <- resilience(tna_ts, window = 30L)
  # res has no _category columns
  expect_error(to_tna(res), "classify_resilience")
})

test_that("to_tna.resilience errors on invalid state name", {
  set.seed(42)
  res <- resilience(tna_ts, window = 30L)
  cls <- classify_resilience(res)
  expect_error(to_tna(cls, state = "nonexistent"), "not found")
})

test_that("to_tna.resilience validates state argument type", {
  set.seed(42)
  res <- resilience(tna_ts, window = 30L)
  cls <- classify_resilience(res)
  expect_error(to_tna(cls, state = 123))
  expect_error(to_tna(cls, state = c("vsi", "cv")))
})

test_that("to_tna.hurst errors when no state columns exist", {
  # Create mock hurst object without state or mf_category columns
  mock_hurst <- tibble::tibble(
    time = 1:10,
    value = rnorm(10),
    hurst = runif(10, 0.3, 0.8),
    r_squared = runif(10, 0.8, 1.0)
  )
  class(mock_hurst) <- c("hurst", "tbl_df", "tbl", "data.frame")
  expect_error(to_tna(mock_hurst), "states = TRUE")
})

test_that("to_tna.hurst errors when requested state column not found", {
  # Create mock hurst object with only state (no mf_category)
  mock_hurst <- tibble::tibble(
    time = 1:10,
    value = rnorm(10),
    hurst = runif(10, 0.3, 0.8),
    r_squared = runif(10, 0.8, 1.0),
    state = sample(c("persistent", "random_walk"), 10, replace = TRUE)
  )
  class(mock_hurst) <- c("hurst", "tbl_df", "tbl", "data.frame")
  expect_error(to_tna(mock_hurst, state = "mf_category"), "not found")
})

test_that("to_tna.hurst_ews returns wide-format tibble", {
  set.seed(42)
  h <- hurst(tna_ts, window = 30L, method = "dfa")
  ews <- detect_hurst_warnings(h)
  result <- to_tna(ews)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_true(all(grepl("^T[0-9]+$", names(result))))
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(is.character(vals))
})

test_that("to_tna.hurst_ews extracts underlying hurst state", {
  set.seed(42)
  h <- hurst(tna_ts, window = 30L, method = "dfa")
  ews <- detect_hurst_warnings(h)
  result <- to_tna(ews, state = "state")
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  vals <- unlist(result[1, ], use.names = FALSE)
  valid_states <- c(
    "strong_antipersistent", "antipersistent", "random_walk",
    "persistent", "strong_persistent"
  )
  expect_true(all(vals %in% valid_states))
})

test_that("to_tna.hurst_ews validates state argument", {
  set.seed(42)
  h <- hurst(tna_ts, window = 30L, method = "dfa")
  ews <- detect_hurst_warnings(h)
  expect_error(to_tna(ews, state = 123))
  expect_error(to_tna(ews, state = c("state", "warning_label")))
})

test_that("to_tna.hurst_ews errors when requested state not found", {
  # Create mock hurst_ews object with only warning_label (no state)
  mock_ews <- tibble::tibble(
    time = 1:10,
    value = rnorm(10),
    warning_label = sample(c("none", "low", "moderate"), 10, replace = TRUE)
  )
  class(mock_ews) <- c("hurst_ews", "tbl_df", "tbl", "data.frame")
  expect_error(to_tna(mock_ews, state = "state"), "not found")
})

test_that("to_tna.hurst_ews errors when no state columns exist", {
  mock_ews <- tibble::tibble(
    time = 1:10,
    value = rnorm(10),
    hurst = runif(10)
  )
  class(mock_ews) <- c("hurst_ews", "tbl_df", "tbl", "data.frame")
  expect_error(to_tna(mock_ews), "no state columns")
})

test_that("to_tna.multi_ews errors when classification is missing", {
  # Create mock expanding multi_ews without classification attribute
  mock_ews <- tibble::tibble(
    time = 1:10,
    metric = rep("test", 10),
    score = runif(10)
  )
  attr(mock_ews, "method") <- "expanding"
  attr(mock_ews, "classification") <- NULL
  class(mock_ews) <- c("multi_ews", "tbl_df", "tbl", "data.frame")
  expect_error(to_tna(mock_ews), "No classification")
})

test_that("to_tna.multi_ews errors when classification has no state column", {
  mock_ews <- tibble::tibble(
    time = 1:10,
    metric = rep("test", 10),
    score = runif(10)
  )
  attr(mock_ews, "method") <- "expanding"
  attr(mock_ews, "classification") <- tibble::tibble(
    time = 1:5,
    count = rep(1, 5)
  )
  class(mock_ews) <- c("multi_ews", "tbl_df", "tbl", "data.frame")
  expect_error(to_tna(mock_ews), "No classification")
})

test_that("to_tna.trend errors when state column is missing", {
  mock_trend <- tibble::tibble(
    time = 1:10,
    value = rnorm(10),
    metric = runif(10)
  )
  class(mock_trend) <- c("trend", "tbl_df", "tbl", "data.frame")
  expect_error(to_tna(mock_trend), "state")
})

test_that("to_tna.changepoint extracts regime states", {
  cp <- detect_cpts(tna_ts, method = "pelt")
  result <- to_tna(cp, state = "regime")
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(all(grepl("^regime_", vals)))
})

test_that("to_tna.changepoint extracts segment states", {
  cp <- detect_cpts(tna_ts, method = "pelt")
  result <- to_tna(cp, state = "segment")
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(all(grepl("^segment_", vals)))
})

test_that("to_tna.changepoint errors on invalid column name", {
  cp <- detect_cpts(tna_ts, method = "pelt")
  expect_error(to_tna(cp, state = "nonexistent"), "not found")
})

test_that("to_tna.potential returns wide-format tibble from rolling analysis", {
  set.seed(42)
  pot <- potential_analysis(tna_ts, window = 50L)
  result <- to_tna(pot)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  expect_true(all(grepl("^T[0-9]+$", names(result))))
  vals <- unlist(result[1, ], use.names = FALSE)
  expect_true(is.character(vals))
  valid_states <- c("flat", "unimodal", "bimodal", "multimodal")
  expect_true(all(vals %in% valid_states))
})

test_that("to_tna.potential errors on global (non-rolling) analysis", {
  set.seed(42)
  pot_global <- potential_analysis(tna_ts)
  expect_error(to_tna(pot_global), "rolling")
})

test_that("to_tna.potential handles all n_wells categories including 0", {
  # Mock a potential object with rolling data that includes 0 wells
  mock_pot <- list(
    rolling = tibble::tibble(
      time = 1:6,
      value = rnorm(6),
      n_wells = c(NA_integer_, 0L, 1L, 2L, 3L, 1L),
      dominant_well_location = rep(NA_real_, 6)
    )
  )
  class(mock_pot) <- c("potential", "list")
  result <- to_tna(mock_pot)
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 1L)
  vals <- unlist(result[1, ], use.names = FALSE)
  # NA is stripped, so we should have 5 columns: flat, unimodal, bimodal, multimodal, unimodal
  expect_equal(ncol(result), 5L)
  expect_equal(vals, c("flat", "unimodal", "bimodal", "multimodal", "unimodal"))
})
