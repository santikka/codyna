# Tests for check.R uncovered lines

test_that("check_string returns early when argument is missing", {
  # Line 94: check_string with missing argument returns NULL silently
  f <- function(x) codyna:::check_string(x)
  result <- f()
  expect_null(result)
})
