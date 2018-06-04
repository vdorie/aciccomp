context("dgp_2017 arguments")

test_that("parameter argument only accepts valid input", {
  expect_error(dgp_2017(0, 1))
  expect_error(dgp_2017(33, 1))
  expect_error(dgp_2017(c(1, 2), 1))
  expect_error(dgp_2017("not-a-number", 1))
  expect_error(dgp_2017(NULL, 1))
})

## random.seed is pretty lenient with seed args, so we are too
test_that("random.seed argument only accepts valid input", {
  expect_error(dgp_2017(1, "not-a-number"))
  expect_error(dgp_2017(1, list(no_seed_arg = 10)))
  expect_error(dgp_2017(1, 0L))
  expect_error(dgp_2017(1, 251L))
})

