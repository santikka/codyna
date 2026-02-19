# Tests for utilities.R uncovered lines

test_that("roll passes time to fun when time is provided", {
  set.seed(42)
  values <- rnorm(20)
  time_vec <- seq_along(values)
  # A function that takes both values and time arguments
  slope_fun <- function(values, time) {
    if (length(values) < 2) return(NA_real_)
    fit <- lm(values ~ time)
    coef(fit)[2L]
  }
  result <- codyna:::roll(
    slope_fun, values, time_vec, window = 5L, align = "center"
  )
  expect_length(result, 20L)
  # The middle values should not be NA
  expect_false(all(is.na(result)))
  # Edge values should be NA due to center alignment
  expect_true(is.na(result[1]))
  expect_true(is.na(result[20]))
})

test_that("rollmean works with right alignment", {
  values <- as.numeric(1:10)
  result <- codyna:::rollmean(values, window = 3L, align = "right")
  expect_length(result, 10L)
  # First two values should be NA for right alignment with window 3
  expect_true(is.na(result[1]))
  expect_true(is.na(result[2]))
  # Third value should be mean of 1,2,3 = 2
  expect_equal(result[3], 2)
  # Last value should be mean of 8,9,10 = 9
  expect_equal(result[10], 9)
})

test_that("rollmean works with left alignment", {
  values <- as.numeric(1:10)
  result <- codyna:::rollmean(values, window = 3L, align = "left")
  expect_length(result, 10L)
  # Last two values should be NA for left alignment with window 3
  expect_true(is.na(result[9]))
  expect_true(is.na(result[10]))
  # First value should be mean of 1,2,3 = 2
  expect_equal(result[1], 2)
})

test_that("get_cols works with numeric column indices via bang-bang", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  # Use !! to inject a numeric vector so it goes through eval_tidy path
  idx <- c(1L, 3L)
  expr <- rlang::quo(!!idx)
  result <- codyna:::get_cols(expr, df)
  expect_equal(result, c("a", "c"))
})

test_that("get_cols errors for out-of-range numeric indices via bang-bang", {
  df <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  idx <- c(1L, 5L)
  expr <- rlang::quo(!!idx)
  expect_error(
    codyna:::get_cols(expr, df),
    "select columns that don't exist"
  )
})

test_that("get_cols errors for non-character non-numeric type via bang-bang", {
  df <- data.frame(a = 1:3, b = 4:6)
  # Inject a logical vector which is neither character nor numeric
  vals <- c(TRUE, FALSE)
  expr <- rlang::quo(!!vals)
  expect_error(
    codyna:::get_cols(expr, df),
    "character.*vector.*integer.*vector"
  )
})
