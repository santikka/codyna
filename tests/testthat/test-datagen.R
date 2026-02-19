# --------------------------------------------------------------------------
# Tests for generate_ts_data(), generate_tipping_data(), compare_ts()
# --------------------------------------------------------------------------

# --- generate_ts_data() ---------------------------------------------------

test_that("generate_ts_data returns list with data and plot elements", {
  set.seed(42)
  out <- generate_ts_data(
    n_individuals = 1, data_length = 200,
    data_type = "resilience", generate_plot = TRUE
  )
  expect_type(out, "list")
  expect_named(out, c("data", "plot"))
  expect_s3_class(out$data, "data.frame")
  expect_s3_class(out$plot, "gg")
})

test_that("generate_ts_data returns data.frame when generate_plot = FALSE", {
  set.seed(42)
  df <- generate_ts_data(
    n_individuals = 1, data_length = 200,
    data_type = "resilience", generate_plot = FALSE
  )
  expect_s3_class(df, "data.frame")
  expect_true("time" %in% names(df))
  expect_true("value" %in% names(df))
  expect_true("true_phase" %in% names(df))
  expect_equal(nrow(df), 200L)
})

test_that("generate_ts_data resilience mode has correct length and columns", {
  set.seed(42)
  df <- generate_ts_data(
    n_individuals = 2, data_length = 300,
    data_type = "resilience", generate_plot = FALSE
  )
  expect_equal(nrow(df), 600L)
  expect_true("id" %in% names(df))
  expect_s3_class(df$true_phase, "factor")
  expect_s3_class(df$id, "factor")
})

test_that("generate_ts_data clustered mode returns correct output", {
  set.seed(42)
  df <- generate_ts_data(
    data_length = 200, data_type = "clustered",
    n_levels = 3, generate_plot = FALSE
  )
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 200L)
  expect_true("quantile_class" %in% names(df))
  expect_s3_class(df$true_phase, "factor")
  expect_equal(length(levels(df$true_phase)), 3L)
})

test_that("generate_ts_data different seeds produce different data", {
  df1 <- generate_ts_data(
    data_length = 100, data_type = "clustered",
    n_levels = 3, seed = 1, generate_plot = FALSE
  )
  df2 <- generate_ts_data(
    data_length = 100, data_type = "clustered",
    n_levels = 3, seed = 99, generate_plot = FALSE
  )
  expect_false(identical(df1$value, df2$value))
})

test_that("generate_ts_data rejects invalid inputs", {
  expect_error(
    generate_ts_data(data_length = -10, generate_plot = FALSE)
  )
  expect_error(
    generate_ts_data(n_individuals = 0, generate_plot = FALSE)
  )
  expect_error(
    generate_ts_data(data_type = "invalid", generate_plot = FALSE)
  )
  expect_error(
    generate_ts_data(
      data_type = "clustered", n_levels = 10, generate_plot = FALSE
    )
  )
})

test_that("generate_ts_data errors when results_df provided without results_col", {
  expect_error(
    generate_ts_data(
      data_length = 100, generate_plot = TRUE,
      results_df = data.frame(x = 1:10)
    )
  )
})

# --- generate_tipping_data() ----------------------------------------------

test_that("generate_tipping_data returns data.frame with correct dimensions", {
  tp <- generate_tipping_data(n_time = 100, n_vars = 3, tipping_point = 50)
  expect_s3_class(tp, "data.frame")
  expect_equal(nrow(tp), 100L)
  expect_equal(ncol(tp), 4L)
  expect_true("Time" %in% names(tp))
  expect_true(all(paste0("VAR", 1:3) %in% names(tp)))
})

test_that("generate_tipping_data ou model works", {
  generate_tipping_data(
    n_time = 50, n_vars = 2, tipping_point = 25, model = "ou"
  ) |>
    expect_error(NA)
})

test_that("generate_tipping_data ar model works", {
  generate_tipping_data(
    n_time = 50, n_vars = 2, tipping_point = 25, model = "ar"
  ) |>
    expect_error(NA)
})

test_that("generate_tipping_data respects custom tipping and saturation points", {
  tp <- generate_tipping_data(
    n_time = 200, n_vars = 5,
    tipping_point = 50, saturation_point = 150
  )
  expect_equal(nrow(tp), 200L)
  expect_equal(ncol(tp), 6L)
  expect_true(all(vapply(tp[, paste0("VAR", 1:5)], is.numeric, logical(1L))))
})

test_that("generate_tipping_data values are numeric", {
  tp <- generate_tipping_data(n_time = 80, n_vars = 4, tipping_point = 40)
  var_cols <- paste0("VAR", 1:4)
  is_num <- vapply(tp[, var_cols], is.numeric, logical(1L))
  expect_true(all(is_num))
})

test_that("generate_tipping_data rejects invalid inputs", {
  expect_error(
    generate_tipping_data(n_time = 5, n_vars = 2, tipping_point = 3)
  )
  expect_error(
    generate_tipping_data(n_time = 100, n_vars = 2, model = "invalid")
  )
  expect_error(
    generate_tipping_data(n_time = 100, n_vars = 0, tipping_point = 50)
  )
})

# --- compare_ts() ---------------------------------------------------------

test_that("compare_ts single panel returns ggplot", {
  set.seed(42)
  df <- generate_ts_data(
    data_length = 200, data_type = "clustered",
    n_levels = 3, generate_plot = FALSE
  )
  p <- compare_ts(
    data = df, ts_col = "value", time_col = "time",
    state1 = "true_phase", title1 = "Test Panel"
  )
  expect_s3_class(p, "gg")
})

test_that("compare_ts two panel returns patchwork", {
  set.seed(42)
  df <- generate_ts_data(
    data_length = 200, data_type = "clustered",
    n_levels = 3, generate_plot = FALSE
  )
  p <- compare_ts(
    data = df, ts_col = "value", time_col = "time",
    state1 = "true_phase", state2 = "quantile_class",
    title1 = "True", title2 = "Quantile"
  )
  expect_s3_class(p, "patchwork")
})

test_that("compare_ts infers time index when time_col is NULL", {
  set.seed(42)
  df <- generate_ts_data(
    data_length = 100, data_type = "clustered",
    n_levels = 2, generate_plot = FALSE
  )
  compare_ts(
    data = df, ts_col = "value", time_col = NULL,
    state1 = "true_phase"
  ) |>
    expect_error(NA)
})

test_that("compare_ts errors on missing columns", {
  df <- data.frame(time = 1:10, value = rnorm(10), state = rep("A", 10))
  expect_error(
    compare_ts(data = df, ts_col = "value", time_col = "time",
               state1 = "nonexistent")
  )
})

test_that("compare_ts errors when data is not a data.frame", {
  expect_error(
    compare_ts(data = 1:10, ts_col = "value", state1 = "state")
  )
})

test_that("compare_ts errors when data2 row count differs", {
  df1 <- data.frame(time = 1:10, value = rnorm(10), s1 = rep("A", 10))
  df2 <- data.frame(s2 = rep("B", 5))
  expect_error(
    compare_ts(
      data = df1, data2 = df2, ts_col = "value", time_col = "time",
      state1 = "s1", state2 = "s2"
    )
  )
})

# --------------------------------------------------------------------------
# generate_ts_data: results_df path (lines 162-172)
# --------------------------------------------------------------------------

test_that("generate_ts_data with results_df and results_col produces plot", {
  set.seed(42)
  # First generate data
  df <- generate_ts_data(
    n_individuals = 1, data_length = 200,
    data_type = "resilience", generate_plot = FALSE
  )
  # Create an external results data.frame with matching rows
  results <- data.frame(
    time = df$time,
    value = df$value,
    predicted = sample(
      c("Low", "Medium", "High"), nrow(df), replace = TRUE
    )
  )
  out <- generate_ts_data(
    n_individuals = 1, data_length = 200,
    data_type = "resilience", generate_plot = TRUE,
    results_df = results, results_col = "predicted"
  )
  expect_type(out, "list")
  expect_named(out, c("data", "plot"))
  expect_s3_class(out$plot, "patchwork")
})

test_that("generate_ts_data errors when results_col not found in results_df", {
  set.seed(42)
  df_dummy <- data.frame(x = 1:200)
  expect_error(
    generate_ts_data(
      n_individuals = 1, data_length = 200,
      data_type = "resilience", generate_plot = TRUE,
      results_df = df_dummy, results_col = "nonexistent"
    )
  )
})

# --------------------------------------------------------------------------
# generate_ts_data: clustered mode with plot (lines 174-180)
# --------------------------------------------------------------------------

test_that("generate_ts_data clustered mode with generate_plot = TRUE", {
  set.seed(42)
  out <- generate_ts_data(
    n_individuals = 1, data_length = 200,
    data_type = "clustered", n_levels = 3,
    generate_plot = TRUE
  )
  expect_type(out, "list")
  expect_named(out, c("data", "plot"))
  expect_s3_class(out$plot, "patchwork")
  expect_true("quantile_class" %in% names(out$data))
})

# --------------------------------------------------------------------------
# generate_ts_data: multiple individuals with clustered plot and facets
# --------------------------------------------------------------------------

test_that("generate_ts_data clustered with multiple individuals includes facets", {
  set.seed(42)
  out <- generate_ts_data(
    n_individuals = 2, data_length = 200,
    data_type = "clustered", n_levels = 3,
    generate_plot = TRUE
  )
  expect_type(out, "list")
  expect_s3_class(out$plot, "patchwork")
  expect_true("id" %in% names(out$data))
})

# --------------------------------------------------------------------------
# generate_tipping_data: AR model with saturation (line 329)
# --------------------------------------------------------------------------

test_that("generate_tipping_data ar model reaches saturation branch", {
  # Set saturation_point < n_time so some iterations have i > saturation_point
  tp <- generate_tipping_data(
    n_time = 100, n_vars = 3, tipping_point = 30,
    stability_strength = 0.8, forcing_strength = 0.01,
    saturation_point = 60, model = "ar"
  )
  expect_s3_class(tp, "data.frame")
  expect_equal(nrow(tp), 100L)
  # Verify all values are finite numeric
  var_cols <- paste0("VAR", 1:3)
  expect_true(all(vapply(tp[, var_cols], function(x) all(is.finite(x)),
                         logical(1L))))
})

# --------------------------------------------------------------------------
# compare_ts: data2 is not a data.frame (line 460)
# --------------------------------------------------------------------------

test_that("compare_ts errors when data2 is not a data.frame", {
  set.seed(42)
  df <- data.frame(time = 1:10, value = rnorm(10), s1 = rep("A", 10))
  expect_error(
    compare_ts(
      data = df, data2 = list(s2 = rep("B", 10)),
      ts_col = "value", time_col = "time",
      state1 = "s1", state2 = "s2"
    )
  )
})

# --------------------------------------------------------------------------
# compare_ts: state2 not found in data2 (line 462)
# --------------------------------------------------------------------------

test_that("compare_ts errors when state2 column not in data2", {
  set.seed(42)
  df1 <- data.frame(time = 1:10, value = rnorm(10), s1 = rep("A", 10))
  df2 <- data.frame(other_col = rep("B", 10))
  expect_error(
    compare_ts(
      data = df1, data2 = df2, ts_col = "value", time_col = "time",
      state1 = "s1", state2 = "s2"
    )
  )
})

# --------------------------------------------------------------------------
# compare_ts: valid data2 merge path (line 467)
# --------------------------------------------------------------------------

test_that("compare_ts merges state2 from data2 successfully", {
  set.seed(42)
  df1 <- data.frame(time = 1:10, value = rnorm(10), s1 = rep("A", 10))
  df2 <- data.frame(s2 = rep("B", 10))
  p <- compare_ts(
    data = df1, data2 = df2, ts_col = "value", time_col = "time",
    state1 = "s1", state2 = "s2"
  )
  expect_s3_class(p, "patchwork")
})

# --------------------------------------------------------------------------
# compare_ts: state2 not in data when data2 is NULL (line 470)
# --------------------------------------------------------------------------

test_that("compare_ts errors when state2 not found in data (no data2)", {
  set.seed(42)
  df <- data.frame(time = 1:10, value = rnorm(10), s1 = rep("A", 10))
  expect_error(
    compare_ts(
      data = df, ts_col = "value", time_col = "time",
      state1 = "s1", state2 = "nonexistent"
    )
  )
})

# --------------------------------------------------------------------------
# compare_ts: facet_by parameter (lines 507-510)
# --------------------------------------------------------------------------

test_that("compare_ts with facet_by adds facets", {
  set.seed(42)
  df <- data.frame(
    time = rep(1:10, 2),
    value = rnorm(20),
    state = rep("A", 20),
    group = rep(c("G1", "G2"), each = 10)
  )
  p <- compare_ts(
    data = df, ts_col = "value", time_col = "time",
    state1 = "state", facet_by = "group"
  )
  expect_s3_class(p, "gg")
  # Verify facet was applied by checking ggplot layers
  expect_true(length(p$facet$params$facets) > 0L ||
                !is.null(p$facet$params$facets))
})

test_that("compare_ts with facet_by and two panels works", {
  set.seed(42)
  df <- data.frame(
    time = rep(1:10, 2),
    value = rnorm(20),
    s1 = rep("A", 20),
    s2 = rep("B", 20),
    group = rep(c("G1", "G2"), each = 10)
  )
  p <- compare_ts(
    data = df, ts_col = "value", time_col = "time",
    state1 = "s1", state2 = "s2",
    facet_by = "group"
  )
  expect_s3_class(p, "patchwork")
})

# --------------------------------------------------------------------------
# compare_ts: missing columns error including facet_by (line 477)
# --------------------------------------------------------------------------

test_that("compare_ts errors when facet_by column is missing", {
  set.seed(42)
  df <- data.frame(time = 1:10, value = rnorm(10), s1 = rep("A", 10))
  expect_error(
    compare_ts(
      data = df, ts_col = "value", time_col = "time",
      state1 = "s1", facet_by = "nonexistent_group"
    )
  )
})

# --------------------------------------------------------------------------
# datagen_resilience_: n_stable_levels <= 0 (lines 554-555)
# --------------------------------------------------------------------------

test_that("datagen_resilience_ errors for n_stable_levels = 0", {
  set.seed(42)
  expect_error(
    codyna:::datagen_resilience_(200, 0, 100, 20)
  )
})

# --------------------------------------------------------------------------
# datagen_resilience_: n_stable_levels = 1 (lines 558-560)
# --------------------------------------------------------------------------

test_that("datagen_resilience_ works with n_stable_levels = 1", {
  set.seed(42)
  df <- generate_ts_data(
    n_individuals = 1, data_length = 200,
    data_type = "resilience", n_stable_levels = 1,
    generate_plot = FALSE
  )
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 200L)
  # With 1 stable level, we should have "Stable" in the phases
  expect_true("Stable" %in% levels(df$true_phase))
})

# --------------------------------------------------------------------------
# datagen_resilience_: n_stable_levels = 2 (lines 563-566)
# --------------------------------------------------------------------------

test_that("datagen_resilience_ works with n_stable_levels = 2", {
  set.seed(42)
  df <- generate_ts_data(
    n_individuals = 1, data_length = 200,
    data_type = "resilience", n_stable_levels = 2,
    generate_plot = FALSE
  )
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 200L)
  # With 2 stable levels, should have "Stable Low" and "Stable High"
  expect_true("Stable Low" %in% levels(df$true_phase))
  expect_true("Stable High" %in% levels(df$true_phase))
})

# --------------------------------------------------------------------------
# datagen_resilience_: n_stable_levels >= 4 (lines 572-575)
# --------------------------------------------------------------------------

test_that("datagen_resilience_ works with n_stable_levels = 4", {
  set.seed(42)
  df <- generate_ts_data(
    n_individuals = 1, data_length = 300,
    data_type = "resilience", n_stable_levels = 4,
    generate_plot = FALSE
  )
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 300L)
  # With 4 stable levels, should have "Stable 1" through "Stable 4"
  stable_levels <- paste("Stable", 1:4)
  expect_true(all(stable_levels %in% levels(df$true_phase)))
})

# --------------------------------------------------------------------------
# datagen_cluster_proportions_: explicit proportions (lines 531-538)
# --------------------------------------------------------------------------

test_that("datagen_cluster_proportions_ uses provided proportions", {
  set.seed(42)
  props <- c(0.3, 0.3, 0.4)
  result <- codyna:::datagen_cluster_proportions_(3, props, 0.1, TRUE)
  expect_equal(result, props)
})

test_that("datagen_cluster_proportions_ errors on wrong length proportions", {
  set.seed(42)
  expect_error(
    codyna:::datagen_cluster_proportions_(3, c(0.5, 0.5), 0.1, TRUE)
  )
})

test_that("datagen_cluster_proportions_ errors when proportions don't sum to 1", {
  set.seed(42)
  expect_error(
    codyna:::datagen_cluster_proportions_(3, c(0.3, 0.3, 0.3), 0.1, TRUE)
  )
})

# --------------------------------------------------------------------------
# datagen_cluster_proportions_: non-randomize path (line 540)
# --------------------------------------------------------------------------

test_that("datagen_cluster_proportions_ returns equal proportions when not randomized", {
  set.seed(42)
  result <- codyna:::datagen_cluster_proportions_(4, NULL, 0.1, FALSE)
  expect_equal(result, rep(0.25, 4))
})

# --------------------------------------------------------------------------
# datagen_cluster_proportions_: warning path (lines 546-547)
# --------------------------------------------------------------------------

test_that("datagen_cluster_proportions_ warns when constraints unsatisfiable", {
  set.seed(42)
  # min_cluster_size = 0.99 is impossible for multiple clusters:
  # each cluster must be >= 99% of total, which is impossible for n > 1
  expect_warning(
    codyna:::datagen_cluster_proportions_(3, NULL, 0.99, TRUE)
  )
})

# --------------------------------------------------------------------------
# generate_ts_data: results_df with multi-individual and facet_by
# --------------------------------------------------------------------------

test_that("generate_ts_data with results_df and multiple individuals", {
  set.seed(42)
  df <- generate_ts_data(
    n_individuals = 2, data_length = 200,
    data_type = "resilience", generate_plot = FALSE
  )
  results <- data.frame(
    time = df$time,
    value = df$value,
    id = df$id,
    predicted = sample(
      c("Low", "High"), nrow(df), replace = TRUE
    )
  )
  out <- generate_ts_data(
    n_individuals = 2, data_length = 200,
    data_type = "resilience", generate_plot = TRUE,
    results_df = results, results_col = "predicted"
  )
  expect_type(out, "list")
  expect_s3_class(out$plot, "patchwork")
})
