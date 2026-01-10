test_that("EWS can be computed with rolling window", {
  detect_warnings(mock_ts, method = "rolling") |>
    expect_error(NA)
})

test_that("EWS can be computed with expanding window", {
  detect_warnings(mock_ts, method = "expanding") |>
    expect_error(NA)
})

test_that("detrending methods can be applied", {
  methods <- c(
    "gaussian",
    "loess",
    "linear",
    "first-diff"
  )
  for (method in methods) {
    detect_warnings(mock_ts, detrend = method) |>
      expect_error(NA)
  }
})

test_that("EWS can be applied without demeaning", {
  detect_warnings(mock_ts, method = "rolling", demean = FALSE) |>
    expect_error(NA)
  detect_warnings(mock_ts, method = "expanding", demean = FALSE) |>
    expect_error(NA)
})

test_that("classification accounts for the number of metrics", {
  ews <- detect_warnings(mock_ts, "expanding", metrics = "sd")
  cls <- attr(ews, "classification")
  expect_true(length(levels(cls$state)) == 3)
  ews <- detect_warnings(mock_ts, "expanding", metrics = c("sd", "ar1"))
  cls <- attr(ews, "classification")
  expect_true(length(levels(cls$state)) == 4)
})

test_that("EWS results with rolling window can be plotted", {
  ews <- detect_warnings(mock_ts, method = "rolling")
  plot(ews) |>
    expect_error(NA)
})

test_that("EWS results with expanding window can be plotted", {
  ews <- detect_warnings(mock_ts, method = "expanding")
  plot(ews) |>
    expect_error(NA)
})
