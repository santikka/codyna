test_that("sequence indices can be computed", {
  expect_error(
    sequence_indices(mock_sequence),
    NA
  )
  expect_error(
    sequence_indices(engagement[1:200, ]),
    NA
  )
})

test_that("favorable states can be specified", {
  expect_error(
    sequence_indices(mock_sequence, favorable = "A"),
    NA
  )
})
