context("Test placeholder.R")

test_that("placeholder doesn't fall over", {
  expect_equal(placeholder(2), 4)
  expect_is(placeholder(4), "numeric")
})
