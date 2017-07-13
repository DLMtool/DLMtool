
context("runMSE with testOM")

test_that("runMSE works with testOM for all MP classes", {
  OM <- testOM 
  OM@nsim <- 2
  
  # Output 
  expect_is(runMSE(OM, MPs="AvC"), "MSE")
  
  # Size limit 
  expect_is(runMSE(OM, MPs=c("matlenlim", "matlenlim2")), "MSE")
  
  # Harvest Slot limit 
  expect_is(runMSE(OM, MPs=c("slotlim")), "MSE")
  
  # Effort controls
  expect_is(runMSE(OM, MPs=c("curE", "curE75")), "MSE")
  
  # Spatial closures 
  expect_is(runMSE(OM, MPs=c("MRnoreal", "MRreal")), "MSE")
  
})

