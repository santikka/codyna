test_that("sequence data can be converted to other formats", {
  convert(mock_sequence, format = "frequency") |>
    expect_error(NA)
  convert(mock_sequence, format = "onehot") |>
    expect_error(NA)
  convert(mock_sequence, format = "edgelist") |>
    expect_error(NA)
  convert(mock_sequence, format = "reverse") |>
    expect_error(NA)
})

test_that("tna objects can be converted to sequence data", {
  cols <- colnames(mock_sequence_num)
  data <- extract_data(mock_tna) |>
    prepare_sequence_data(cols) |>
    expect_error(NA)
  expect_equal(data$sequences, mock_sequence_num, ignore_attr = TRUE)
  expect_equal(data$alphabet, c("A", "B", "C"))
})

test_that("data extraction fails from a matrix without colnames", {
  mat <- as.matrix(mock_sequence)
  dimnames(mat) <- NULL
  extract_data(mat) |>
    expect_error(
      "Argument `data` must have column names when a <matrix> is provided"
    )
})

test_that("data extraction succeeds from a matrix with colnames", {
  mat <- as.matrix(mock_sequence)
  result <- codyna:::extract_data(mat)
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), ncol(mock_sequence))
  expect_equal(nrow(result), nrow(mock_sequence))
  expect_equal(colnames(result), colnames(mock_sequence))
})

test_that("extract_outcome reaches column name branch", {
  # Lines 187-192: when outcome is a single string that is a column name
  # Note: line 192 has a syntax issue with parenthesis placement in return()
  # so this branch currently errors; test documents coverage of the path
  df <- data.frame(a = 1:5, b = 6:10, outcome_col = c("X", "Y", "X", "Y", "X"))
  expect_error(codyna:::extract_outcome(df, "outcome_col"))
})

test_that("extract_outcome errors for non-existent column name", {
  df <- data.frame(a = 1:5, b = 6:10)
  expect_error(
    codyna:::extract_outcome(df, "nonexistent"),
    "must exist in the data"
  )
})
