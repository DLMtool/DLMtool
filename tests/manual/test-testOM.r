
context("runMSE with testOM")

testthat::test_that("runMSE works with testOM for all MP classes", {
  OM <- testOM 
  OM@nsim <- 2
  
  # Output 
  testthat::expect_is(runMSEnomsg(OM, MPs="AvC"), "MSE")
  
  # Size limit 
  testthat::expect_is(runMSEnomsg(OM, MPs=c("matlenlim", "matlenlim2")), "MSE")
  
  # Harvest Slot limit 
  testthat::expect_is(runMSEnomsg(OM, MPs=c("slotlim")), "MSE")
  
  # Effort controls
  testthat::expect_is(runMSEnomsg(OM, MPs=c("curE", "curE75")), "MSE")
  
  # Spatial closures 
  testthat::expect_is(runMSEnomsg(OM, MPs=c("MRnoreal", "MRreal")), "MSE")
  
})

