test_that("complexity measures can be calculated", {
  complexity(mock_ts) |>
    expect_error(NA)
  complexity(mock_ts, measures = c("autocorrelation", "variance")) |>
    expect_error(NA)
})

test_that("alternative window alignments can be used", {
  complexity(mock_ts, measures = "fluctuation", align = "left") |>
    expect_error(NA)
  complexity(mock_ts, measures = "fluctuation", align = "right") |>
    expect_error(NA)
})
