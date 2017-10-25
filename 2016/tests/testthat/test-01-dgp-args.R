context("dgp_2016 arguments")

test_that("x argument only accepts valid input", {
  expect_error(dgp_2016(NULL, 1, 1))
  expect_error(dgp_2016(list(junk = 10), 1, 1))
  expect_error(dgp_2016(matrix(NA, 0, 0), 1, 1))
})

test_that("parameter argument only accepts valid input", {
  x <- seq_len(5L)
  expect_error(dgp_2016(x, 0, 1))
  expect_error(dgp_2016(x, 78, 1))
  expect_error(dgp_2016(x, c(1, 2), 1))
  expect_error(dgp_2016(x, "not-a-number", 1))
  expect_error(dgp_2016(x, NULL, 1))
})

## random.seed is pretty lenient with seed args, so we are too
test_that("random.seed argument only accepts valid input", {
  x <- seq_len(5L)
  expect_error(dgp_2016(x, 1, "not-a-number"))
  expect_error(dgp_2016(x, 1, list(no_seed_arg = 10)))
  expect_error(dgp_2016(x, 1, 0L))
  expect_error(dgp_2016(x, 1, 101L))
})

