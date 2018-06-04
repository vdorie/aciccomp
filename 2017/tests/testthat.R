if (require(testthat, quietly = TRUE)) {
  require(aciccomp2017)
  test_check("aciccomp2017")
} else {
  cat("package 'testthat' not available; cannot run unit tests\n")
}

