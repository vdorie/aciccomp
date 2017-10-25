context("dgp_2016 results")

test_that("dgp_2016 reproduces past simulations", {

  load(system.file("testData.RData", package = "aciccomp2016"))

  numParameterSets <- 20L
  
  data(input_2016, package = "aciccomp2016")
  
  for (i in numParameterSets) {
    sim <- dgp_2016(input_2016, i, 1L)
    expect_equal(sim$y , sims[[i]]$y, tolerance = 1e-8)
    expect_identical(sim$z, sims[[i]]$z)
  }
})
