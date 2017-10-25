if (require(testthat, quietly = TRUE)) {
  require(aciccomp2016)
  test_check("aciccomp2016")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}

