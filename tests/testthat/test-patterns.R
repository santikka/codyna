test_that("n-gram patterns can be discovered", {
  discover_patterns(mock_sequence, type = "ngram") |>
    expect_error(NA)
})

test_that("gapped pattern can be discovered", {
  discover_patterns(mock_sequence, type = "gapped") |>
    expect_error(NA)
})

test_that("repeated pattern can be discovered", {
  discover_patterns(mock_sequence, type = "repeated") |>
    expect_error(NA)
})

test_that("custom patterns can be discovered", {
  discover_patterns(mock_sequence, pattern = "A->*") |>
    expect_error(NA)
  discover_patterns(mock_sequence, pattern = "*->B") |>
    expect_error(NA)
})

test_that("patterns can be filtered based on states", {
  patterns <- discover_patterns(mock_sequence, start = "B", end = "C")
  expect_true("B->A->C" %in% patterns$pattern)
  expect_true("B->C->C" %in% patterns$pattern)
  patterns <- discover_patterns(mock_sequence, contain = "C|B")
  expect_true("A->C" %in% patterns$pattern)
  expect_true("A->B" %in% patterns$pattern)
  expect_true("C->C" %in% patterns$pattern)
  expect_false("A->A" %in% patterns$pattern)
})

test_that("patterns can be filtered based on frequency and support", {
  patterns <- discover_patterns(mock_sequence, min_freq = 3)
  expect_true(all(patterns$frequency >= 3))
  patterns <- discover_patterns(mock_sequence, min_support = 0.2)
  expect_true(all(patterns$support >= 0.2))
})

test_that("output has zero rows if no patterns are found", {
  patterns <- discover_patterns(mock_sequence, pattern = "A->******->B")
  expect_true(nrow(patterns) == 0L)
  patterns <- discover_patterns(mock_sequence, min_freq = 6)
  expect_true(nrow(patterns) == 0L)
})

test_that("pattern discovery supports different input formats", {
  discover_patterns(engagement[1:100, ]) |>
    expect_error(NA)
  discover_patterns(mock_tna) |>
    expect_error(NA)
  discover_patterns(mock_group_tna) |>
    expect_error(NA)
  discover_patterns(mock_tna_data) |>
    expect_error(NA)
})

test_that("pattern discovery supports last observation as group", {
  mock_sequence_out <- mock_sequence
  mock_sequence_out$T7 <- rep(c("D", "F"), length.out = nrow(mock_sequence))
  patterns <- discover_patterns(mock_sequence_out, group = "last_obs") |>
    expect_error(NA)
  all(c("count_D", "count_F", "chisq") %in% names(patterns)) |>
    expect_true()
})

test_that("pattern discovery supports groups", {
  patterns <- discover_patterns(mock_group_tna) |>
    expect_error(NA)
  all(c("count_Group 1", "count_Group 2") %in% names(patterns)) |>
    expect_true()
})

test_that("discovered patterns can be printed", {
  patterns <- discover_patterns(mock_sequence)
  print(patterns) |>
    capture.output() |>
    expect_error(NA)
})

test_that("discovered patterns can be plotted", {
  patterns1 <- discover_patterns(mock_sequence)
  plot(patterns1) |>
    expect_error(NA)
  grp <- rep(1:2, length.out = nrow(mock_sequence))
  patterns2 <- discover_patterns(mock_sequence, group = grp)
  plot(patterns2) |>
    expect_error(NA)
  plot(patterns2, group = "1", global = FALSE) |>
    expect_error(NA)
  plot(patterns2, group = "1", global = TRUE) |>
    expect_error(NA)
})

test_that("outcomes can be analyzed", {
  fit <- analyze_outcome(
    engagement[idx, ],
    outcome = rep(1:2, length.out = length(idx)),
    len = 1
  ) |>
    expect_error(NA)
})

test_that("outcome analysis supports different input formats", {
  out <-  rep(1:2, length.out = nrow(mock_sequence))
  analyze_outcome(mock_tna, outcome = out, len = 2) |>
    expect_error(NA)
  analyze_outcome(mock_group_tna, outcome = out, len = 2, mixed = FALSE) |>
    expect_error(NA)
  analyze_outcome(mock_tna_data, outcome = out, len = 2) |>
    expect_error(NA)
})

test_that("outcome analysis supports mixed effects", {
  out <-  rep(1:2, length.out = nrow(mock_sequence))
  analyze_outcome(mock_tna_data, outcome = out, group = "group", len = 2) |>
    suppressMessages() |> # singular
    expect_error(NA)
  analyze_outcome(mock_group_tna, outcome = out, len = 2) |>
    suppressMessages() |> # singular
    expect_error(NA)
})

test_that("outcome analysis works with last_obs outcome", {
  # Lines 555-556: resp$last = TRUE path in analyze_outcome
  # Create data where last column has exactly 2 unique outcomes
  mock_seq_out <- mock_sequence
  mock_seq_out$T7 <- rep(c("D", "F"), length.out = nrow(mock_sequence))
  fit <- analyze_outcome(
    mock_seq_out,
    outcome = "last_obs",
    len = 1,
    min_freq = 1,
    min_support = 0.01
  ) |>
    expect_error(NA)
})

test_that("outcome analysis supports reference parameter", {
  # Lines 567-568: reference parameter changes factor levels
  out <- rep(c("A", "B"), length.out = length(idx))
  fit <- analyze_outcome(
    engagement[idx, ],
    outcome = out,
    reference = "B",
    len = 1
  ) |>
    expect_error(NA)
})

# ============================================================================
# Dead code documentation: patterns.R lines 563-564
#
# Lines 563-564 in analyze_outcome() are the `else` branch at line 562:
#   } else {
#     response <- resp$var
#     data[[response]] <- factor(data[[response]])
#   }
#
# ROOT CAUSE: Bug in extract_outcome() at data.R line 192. The code reads:
#   return(list(last = FALSE, outcome = x[[outcome]]), var = outcome)
# The closing parenthesis of list() is BEFORE `var = outcome`, making
# `var = outcome` a second argument to return(). R does not allow
# multi-argument returns and throws an error.
#
# This means the only code path that would set resp$var to a non-NULL value
# (when outcome is a column name) CRASHES instead. All other return paths
# from extract_outcome() set var = NULL:
#   - Missing outcome: var = NULL
#   - "last_obs": var = NULL
#   - External vector (length == nrow): var = NULL
#
# Therefore is.null(resp$var) on line 558 is always TRUE, and the else
# branch (lines 563-564) is never reached.
# ============================================================================

test_that("extract_outcome crashes when outcome is a column name", {
  # Verify the bug in extract_outcome that makes lines 563-564 unreachable
  mock_df <- data.frame(
    T1 = c("A", "B", "A"),
    T2 = c("B", "A", "B"),
    outcome_col = c("X", "Y", "X")
  )
  # extract_outcome with a column name should crash due to misplaced paren
  expect_error(
    codyna:::extract_outcome(mock_df, "outcome_col"),
    "multi-argument returns"
  )
})

test_that("analyze_outcome crashes when outcome is a column name", {
  # Confirm the bug propagates to the public API
  mock_df <- data.frame(
    T1 = c("A", "B", "A", "B", "A"),
    T2 = c("B", "A", "B", "A", "B"),
    T3 = c("A", "B", "A", "B", "A"),
    outcome_col = c("X", "Y", "X", "Y", "X")
  )
  expect_error(
    analyze_outcome(mock_df, outcome = "outcome_col", len = 1),
    "multi-argument returns"
  )
})
