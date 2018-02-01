testthat::context("runMSE with testOM")

testthat::test_that("runMSE works with testOM", {
  OM <- testOM 
  OM@nsim <- 6
  testthat::expect_is(runMSE(OM, silent=TRUE), "MSE")
})

testthat::test_that("runMSE works with testOM - parallel", {
  OM <- testOM 
  OM@nsim <- 96
  DLMtool::setup()
  testthat::expect_is(runMSE(OM, silent=TRUE, parallel = TRUE), "MSE")
})

  
testthat::test_that("runMSE works with testOM - interval = 1", {
  OM <- testOM 
  OM@nsim <- 6
  OM@interval <- 1 
  testthat::expect_is(runMSE(OM, silent=TRUE), "MSE")
})
  


