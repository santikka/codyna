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
  prepare_sequence_data(mock_tna) |>
    expect_error(NA)
  data <- prepare_sequence_data(mock_tna)
  expect_equal(data$sequences, mock_sequence_num, ignore_attr = TRUE)
  expect_equal(data$alphabet, c("A", "B", "C"))
})
