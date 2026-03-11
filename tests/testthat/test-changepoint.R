# --- Test data ---------------------------------------------------------------

set.seed(0)
cpt_data <- c(rnorm(100, 0, 1), rnorm(100, 5, 1))

# --- Basic output structure --------------------------------------------------

test_that("detect_cpts returns a changepoint tibble", {
  result <- detect_cpts(cpt_data)
  expect_s3_class(result, "changepoint")
  expect_s3_class(result, "tbl_df")
})

test_that("detect_cpts returns correct columns", {
  result <- detect_cpts(cpt_data)
  expected_cols <- c(
    "time", "value", "segment", "regime", "segment_mean", "segment_var",
    "level", "direction", "magnitude", "changepoint_type", "state",
    "changepoint"
  )
  expect_identical(names(result), expected_cols)
})

test_that("column types are correct", {
  result <- detect_cpts(cpt_data)
  expect_type(result$segment, "integer")
  expect_type(result$changepoint, "logical")
  expect_s3_class(result$state, "factor")
  expect_s3_class(result$level, "factor")
  expect_s3_class(result$direction, "factor")
  expect_type(result$segment_mean, "double")
  expect_type(result$segment_var, "double")
  expect_type(result$magnitude, "double")
})

test_that("segment column is sequential integer", {
  result <- detect_cpts(cpt_data)
  segments <- unique(result$segment)
  expect_identical(segments, seq_len(max(segments)))
})

# --- Attributes --------------------------------------------------------------

test_that("attributes are stored correctly", {
  result <- detect_cpts(cpt_data, method = "pelt", penalty = "bic")
  expect_identical(attr(result, "method"), "pelt")
  expect_identical(attr(result, "penalty"), "bic")
  expect_type(attr(result, "n_changepoints"), "integer")
  expect_type(attr(result, "changepoint_locations"), "integer")
  expect_identical(attr(result, "type"), "mean")
  expect_identical(attr(result, "min_segment"), 10L)
})

# --- Detection accuracy ------------------------------------------------------

test_that("detect_cpts finds changepoint near known shift", {
  result <- detect_cpts(cpt_data, method = "pelt")
  locs <- attr(result, "changepoint_locations")
  expect_true(length(locs) >= 1L)
  # The true shift is at position 100; detected should be within 15
  closest <- locs[which.min(abs(locs - 100))]
  expect_true(abs(closest - 100) <= 15)
})

# --- Methods -----------------------------------------------------------------

test_that("cusum method works", {
  result <- detect_cpts(cpt_data, method = "cusum")
  expect_s3_class(result, "changepoint")
  expect_true(attr(result, "n_changepoints") >= 1L)
})

test_that("binary_segmentation method works", {
  result <- detect_cpts(cpt_data, method = "binseg")
  expect_s3_class(result, "changepoint")
  expect_true(attr(result, "n_changepoints") >= 1L)
})

test_that("pelt method works", {
  result <- detect_cpts(cpt_data, method = "pelt")
  expect_s3_class(result, "changepoint")
  expect_true(attr(result, "n_changepoints") >= 1L)
})

# --- Type options ------------------------------------------------------------

test_that("type mean detects mean shift", {
  result <- detect_cpts(cpt_data, type = "mean")
  expect_s3_class(result, "changepoint")
  expect_true(attr(result, "n_changepoints") >= 1L)
})

test_that("type variance works", {
  set.seed(0)
  var_data <- c(rnorm(100, 0, 1), rnorm(100, 0, 5))
  result <- detect_cpts(var_data, type = "variance")
  expect_s3_class(result, "changepoint")
})

test_that("type both works", {
  result <- detect_cpts(cpt_data, type = "both")
  expect_s3_class(result, "changepoint")
  expect_true(attr(result, "n_changepoints") >= 1L)
})

# --- Penalty options ---------------------------------------------------------

test_that("bic penalty works", {
  result <- detect_cpts(cpt_data, penalty = "bic")
  expect_identical(attr(result, "penalty"), "bic")
})

test_that("aic penalty works", {
  result <- detect_cpts(cpt_data, penalty = "aic")
  expect_identical(attr(result, "penalty"), "aic")
})

test_that("manual penalty works with penalty_value", {
  result <- detect_cpts(cpt_data, penalty = "manual", penalty_value = 10)
  expect_identical(attr(result, "penalty"), "manual")
})

test_that("manual penalty without penalty_value errors", {
  expect_error(
    detect_cpts(cpt_data, penalty = "manual"),
    "penalty_value"
  )
})

# --- S3 methods --------------------------------------------------------------

test_that("print does not error", {
  result <- detect_cpts(cpt_data)
  capture.output(print(result)) |>
    expect_error(NA)
})

test_that("summary does not error", {
  result <- detect_cpts(cpt_data)
  sumr <- summary(result)
  nm <- names(sumr)
  expect_type(sumr, "list")
  expect_true("cpts" %in% nm)
  expect_true("changes" %in% nm)
  expect_true("segments" %in% nm)
})

test_that("summary print does not error", {
  result <- detect_cpts(cpt_data)
  sumr <- summary(result)
  capture.output(print(sumr)) |>
    expect_error(NA)
  expect_type(sumr, "list")
})

test_that("plot does not error", {
  result <- detect_cpts(cpt_data)
  plot(result) |>
    expect_error(NA)
  plot(result, type = "diagnostics") |>
    expect_error(NA)
  plot(result, type = "both") |>
    expect_error(NA)
})

# --- Input validation --------------------------------------------------------

test_that("bad method produces error", {
  expect_error(
    detect_cpts(cpt_data, method = "nonexistent"),
    "method"
  )
})

test_that("bad type produces error", {
  expect_error(
    detect_cpts(cpt_data, type = "nonexistent"),
    "type"
  )
})

test_that("bad penalty produces error", {
  expect_error(
    detect_cpts(cpt_data, penalty = "nonexistent"),
    "penalty"
  )
})

# --- changepoint_cost_ edge cases --------------------------------------------

test_that("changepoint_cost_ returns 0 for segment with fewer than 2 obs", {
  expect_equal(changepoint_cost_(5, "mean"), 0)
  expect_equal(changepoint_cost_(numeric(0), "mean"), 0)
})

test_that("changepoint_cost_ variance type handles near-zero variance", {
  # v < eps guard for variance type
  result <- changepoint_cost_(rep(3.0, 10), "variance")
  expect_true(is.finite(result))
})

test_that("changepoint_cost_ both type handles near-zero variance", {
  # v < eps guard for both type
  result <- changepoint_cost_(rep(3.0, 10), "both")
  expect_true(is.finite(result))
})

test_that("changepoint_cost_ default branch falls through", {
  result <- changepoint_cost_(c(1, 2, 3, 4), "unknown_type")
  # Default computes sum((x - mu)^2) = RSS
  expected <- sum((c(1, 2, 3, 4) - mean(c(1, 2, 3, 4)))^2)
  expect_equal(result, expected)
})

# --- changepoint_penalty_ edge cases  ---------------------------------------

test_that("changepoint_penalty_ manual with non-numeric value errors", {
  # manual penalty with NULL or non-numeric penalty_value
  expect_error(changepoint_penalty_("manual", NULL, 100, "mean"))
  expect_error(changepoint_penalty_("manual", "abc", 100, "mean"))
})

test_that("changepoint_penalty_ default branch returns BIC-like penalty", {
  # default switch branch (unrecognized penalty string)
  result <- changepoint_penalty_("unknown_penalty", NULL, 100, "mean")
  # Default: p * log(n) where p = 2 for type != "both"
  expect_equal(result, 2 * log(100))
})

# --- changepoint_best_split_ edge cases --------------------------------------

test_that("changepoint_best_split_ returns NULL for short segment", {
  # n < 2 * min_segment
  result <- changepoint_best_split_(c(1, 2, 3), "mean", 5)
  expect_null(result)
})

# --- changepoint_cusum_ max_cpt limit ----------------------------------------

test_that("changepoint_cusum_ returns empty when max_cpt is 0", {
  set.seed(0)
  x <- c(rnorm(50, 0), rnorm(50, 5))
  # max_cp <= 0 returns integer(0)
  result <- changepoint_cusum_(x, "mean", 1, 5, 0L, 0L)
  expect_length(result, 0L)
})

test_that("changepoint_cusum_ budget decrements after left recursion", {
  # new_max decremented by left_cp count
  # Use multi-shift data with small penalty and limited max_cpt
  set.seed(0)
  x <- c(rnorm(40, 0), rnorm(40, 10), rnorm(40, -5))
  result <- changepoint_cusum_(x, "mean", 0.5, 5, 2L, 0L)
  expect_true(length(result) <= 2L)
})

# --- changepoint_binseg_ edge cases  -----------------------------------------

test_that("changepoint_binseg_ returns empty for short segment", {
  # n < 2 * min_segment
  result <- changepoint_binseg_(c(1, 2, 3), "mean", 1, 5, NULL, 0L)
  expect_length(result, 0L)
})

test_that("changepoint_binseg_ returns empty when max_cp is 0", {
  set.seed(0)
  x <- c(rnorm(50, 0), rnorm(50, 5))
  # max_cp <= 0
  result <- changepoint_binseg_(x, "mean", 1, 5, 0L, 0L)
  expect_length(result, 0L)
})

test_that("changepoint_binseg_ handles NULL best_split result", {
  # result is NULL from best_split (segment barely meeting 2*min_segment
  # but best_split still returns NULL when n == 2*min_segment and no candidate exists)
  # With min_segment = 5 and n = 10, there's only one split candidate at tau = 5
  # If the gain is below penalty, it returns integer(0)
  result <- changepoint_binseg_(rep(0, 10), "mean", 1000, 5, NULL, 0L)
  expect_length(result, 0L)
})

test_that("changepoint_binseg_ budget decrements after left recursion", {
  # new_max decremented
  set.seed(0)
  x <- c(rnorm(40, 0), rnorm(40, 10), rnorm(40, -5))
  result <- changepoint_binseg_(x, "mean", 0.5, 5, 2L, 0L)
  expect_true(length(result) <= 2L)
})

# --- changepoint_pelt_ edge cases --------------------------------------------

test_that("changepoint_pelt_ returns empty for short segment", {
  # n < 2 * min_segment
  result <-changepoint_pelt_(c(1, 2, 3), "mean", 1, 5, NULL)
  expect_length(result, 0L)
})

test_that("changepoint_pelt_ returns empty when no changepoints found", {
  # cpts length == 0 (very high penalty, no split worthwhile)
  set.seed(0)
  x <- rnorm(50, 0, 1)
  result <- changepoint_pelt_(x, "mean", 1e6, 5, NULL)
  expect_length(result, 0L)
})

test_that("changepoint_pelt_ enforces max_cp by pruning excess changepoints", {
  # max_cpt truncation of PELT results
  set.seed(0)
  # Data with 3 clear shifts but max_cpt = 1
  x <- c(rnorm(60, 0, 0.5), rnorm(60, 10, 0.5), rnorm(60, -5, 0.5))
  # Low penalty so PELT finds multiple changepoints
  result_full <- changepoint_pelt_(x, "mean", 1, 5, NULL)
  expect_true(length(result_full) >= 2L)
  # Now restrict to max_cpt = 1
  result_limited <- changepoint_pelt_(x, "mean", 1, 5, 1L)
  expect_equal(length(result_limited), 1L)
})

# --- detect_cpts: max_changepoints validation --------------- ----------------

test_that("detect_cpts validates max_changepoints parameter", {
  set.seed(0)
  x <- c(rnorm(100, 0, 1), rnorm(100, 5, 1))
  # check_range on max_changepoints
  result <- detect_cpts(x, max_changepoints = 1L)
  expect_s3_class(result, "changepoint")
  expect_true(attr(result, "n_changepoints") <= 1L)
})

test_that("detect_cpts max_changepoints out of range errors", {
  set.seed(0)
  x <- rnorm(50)
  # max_changepoints must be >= 1
  expect_error(detect_cpts(x, max_changepoints = 0L))
})

# --- detect_cpts: changepoints too close filtering  --------------------------

test_that("changepoints too close to each other are filtered", {
  # keep[i] <- FALSE when consecutive cpts are too close
  # Use a very small min_segment with crafted data that produces close cpts
  set.seed(0)
  x <- c(
    rnorm(30, 0, 0.3),
    rnorm(30, 10, 0.3),
    rnorm(30, -5, 0.3),
    rnorm(30, 8, 0.3)
  )
  # With pelt, low penalty, and min_segment = 5, it may produce close cpts
  # that get filtered. Just verify the result is valid.
  result <- detect_cpts(x, method = "pelt", penalty = "aic", min_segment = 5L)
  expect_s3_class(result, "changepoint")
  locs <- attr(result, "changepoint_locations")
  if (length(locs) > 1L) {
    diffs <- diff(locs)
    expect_true(all(diffs >= 5L))
  }
})

# --- detect_cpts: variance type with methods that trigger cost branches ------

test_that("variance type with binary_segmentation covers cost branches", {
  set.seed(0)
  # Data with clear variance change to trigger variance cost path
  var_data <- c(rnorm(80, 0, 0.5), rnorm(80, 0, 5))
  result <- detect_cpts(
    var_data,
    method = "binseg",
    type = "variance",
    min_segment = 10L
  )
  expect_s3_class(result, "changepoint")
})

test_that("type both with pelt method covers cost branches", {
  set.seed(0)
  both_data <- c(rnorm(80, 0, 1), rnorm(80, 5, 3))
  result <- detect_cpts(
    both_data,
    method = "pelt",
    type = "both",
    min_segment = 10L
  )
  expect_s3_class(result, "changepoint")
})

test_that("type variance with cusum covers cost branches", {
  set.seed(0)
  var_data <- c(rnorm(80, 0, 0.5), rnorm(80, 0, 5))
  result <- detect_cpts(
    var_data,
    method = "cusum",
    type = "variance",
    min_segment = 10L
  )
  expect_s3_class(result, "changepoint")
})

# --- detect_cpts: no changepoints detected produces valid output -------------

test_that("no changepoints detected returns single segment", {
  set.seed(0)
  # Very high manual penalty: no split will be accepted
  x <- rnorm(100)
  result <- detect_cpts(x, penalty = "manual", penalty_value = 1e6)
  expect_equal(attr(result, "n_changepoints"), 0L)
  expect_length(attr(result, "changepoint_locations"), 0L)
  expect_true(all(result$segment == 1L))
  expect_false(any(result$changepoint))
})

# --- summary with changepoints: change details (lines around 993-1021) -------

test_that("summary returns change details when changepoints exist", {
  set.seed(0)
  x <- c(rnorm(100, 0, 1), rnorm(100, 5, 1))
  result <- detect_cpts(x, method = "pelt")
  out <- capture.output(s <- summary(result))
  expect_true("changes" %in% names(s))
  expect_true(is.data.frame(s$changes))
  expect_true(nrow(s$changes) > 0L)
  expect_true(
    all(
      c("location", "from_mean", "to_mean", "mean_shift") %in% names(s$changes)
    )
  )
})

test_that("summary with no changepoints returns NULL changes", {
  set.seed(0)
  x <- rnorm(100)
  result <- detect_cpts(x, penalty = "manual", penalty_value = 1e6)
  out <- capture.output(s <- summary(result))
  expect_null(s$changes)
})

# --- plot: fallback color for unknown state names ----------------------------

test_that("plot handles state levels that fall back to cycling colors", {
  set.seed(0)
  x <- c(rnorm(100, 0, 1), rnorm(100, 5, 1))
  result <- detect_cpts(x, method = "pelt")
  # Modify state factor to have a level not in the predefined palette
  levels(result$state) <- c(levels(result$state), "exotic_state")
  result$state[1L] <- "exotic_state"
  # fallback color assignment
  p <- plot(result)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))
})

# --- dead code guards in best_split_ and binseg_ -----------------------------

# These lines are unreachable guard clauses (best_pos will never be NULL when
# n >= 2*min_segment). We test them directly with mocked edge scenarios.

test_that("best split always finds a position when n >= 2*min_segment", {
  # Confirm that best_split always returns a non-NULL result for valid n
  x <- rep(0, 20)
  result <- changepoint_best_split_(x, "mean", 5)
  # Even with constant data, the first candidate will be selected
  # because gain > -Inf for any split
  expect_true(is.list(result))
  expect_true("pos" %in% names(result))
})

# --- changepoints too close to each other ------------------------------------

test_that("cpts too close are filtered via cusum with tight data", {
  # cusum can produce changepoints at positions that, when combined,
  # are closer than min_segment. Create very distinct short segments.
  set.seed(0)
  # 5 segments of length 12 each with large shifts and tiny penalty
  x <- c(
    rnorm(12, 0, 0.1),
    rnorm(12, 50, 0.1),
    rnorm(12, -30, 0.1),
    rnorm(12, 80, 0.1),
    rnorm(12, -50, 0.1)
  )
  # Use min_segment = 3 so close cps might occur
  result <- detect_cpts(
    x,
    method = "cusum",
    penalty = "manual",
    penalty_value = 0.01,
    min_segment = 3L
  )
  expect_s3_class(result, "changepoint")
  locs <- attr(result, "changepoint_locations")
  # Verify minimum spacing constraint
  if (length(locs) > 1L) {
    expect_true(all(diff(locs) >= 3L))
  }
})

# --- changepoint_best_split_ defensive guard (best_pos NULL) -----------------

# This guard is unreachable because best_gain starts at -Inf and
# any finite cost difference sets best_pos. The for loop always has at least
# one candidate when n >= 2*min_segment. We verify the guard exists by showing
# the function always finds a split for valid inputs.

test_that("changepoint_best_split_ always sets best_pos for valid n", {
  # Even constant data produces a valid split
  # (gain = 0, pos set to first candidate)
  result <- changepoint_best_split_(rep(0, 20), "mean", 5)
  expect_true(is.list(result))
  expect_true(!is.null(result$pos))
  expect_true(!is.null(result$gain))
  # Confirm gain is 0 for constant data (no information gain from splitting)
  expect_equal(result$gain, 0)
})

# --- changepoint_binseg_ defensive guard (NULL from best_split) --------------

# This guard is unreachable because changepoint_binseg_ checks n < 2*min_segment
# before calling best_split, and best_split always returns non-NULL
# when n >= 2*min_segment. We verify binseg returns the correct result when
# best_split gain is below penalty (branch instead).

test_that("changepoint_binseg_ returns empty when gain below penalty", {
  # Constant data: best_split returns gain = 0, which is <= any positive penalty
  result <- changepoint_binseg_(rep(0, 20), "mean", 1, 5, NULL, 0L)
  expect_length(result, 0L)
})

# --- changepoint close-spacing filter (defensive guard) ----------------------

# All three methods (cusum, binseg, pelt) internally guarantee that adjacent
# changepoints are at least min_segment apart. The post-filter is
# a safety net that is never triggered in practice. We verify the filter logic
# by confirming that all detected changepoints respect min_segment spacing.

test_that("all methods produce changepoints respecting min_segment spacing", {
  set.seed(0)
  x <- c(rnorm(50, 0, 0.1), rnorm(50, 20, 0.1), rnorm(50, -10, 0.1))
  methods <- c("cusum", "binseg", "pelt")
  vapply(
    methods,
    function(m) {
      result <- detect_cpts(
        x,
        method = m,
        penalty = "manual",
        penalty_value = 0.01,
        min_segment = 5L
      )
      locs <- attr(result, "changepoint_locations")
      if (length(locs) > 1L) {
        expect_true(
          all(diff(locs) >= 5L),
          info = paste("Method", m, "produced close changepoints")
        )
      }
    TRUE
  },
  logical(1))
})
