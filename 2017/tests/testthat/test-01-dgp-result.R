context("dgp_2017 results")

test_that("dgp_2017 reproduces past simulations", {

  load(system.file("testData.RData", package = "aciccomp2017"))

  numParameterSets <- 20L
  
  for (i in numParameterSets) {
    sim <- dgp_2017(i, 1L)
    expect_equal(sim$y , sims[[i]]$y, tolerance = 1e-8)
    expect_identical(sim$z, sims[[i]]$z)
  }
})
