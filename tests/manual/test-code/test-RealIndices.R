
testthat::context("test Real indices in Data object")



# Simulate some data # 
library(DLMtool)
OM <- DLMtool::testOM
Hist <- runMSE(OM, Hist=TRUE)
MPs <- "AvC"
# Grab indices from one sim 
sim <- sample(1:OM@nsim,1)
Bind <- Hist@TSdata$B[sim,]
VBind <- Hist@TSdata$VB[sim,]
SpBInd<- Hist@TSdata$SSB[sim,]
om <- OM 

testthat::test_that("Works with 3 indices", {
  om@cpars$Data <- new("Data")
  om@cpars$Data@Type <- c("Biomass", "VBiomass", 'SpBiomass')
  om@cpars$Data@RInd <- array(rbind(Bind, VBind, SpBInd), dim=c(1, 3, om@nyears))
  MSE <- runMSE(om, MPs=MPs)
  testthat::expect_is(MSE, "MSE")
})
  
testthat::test_that("Works with 2 indices", {
  om@cpars$Data <- new("Data")
  om@cpars$Data@Type <- c("Biomass", 'SpBiomass')
  om@cpars$Data@RInd <- array(rbind(Bind, SpBInd), dim=c(1, 2, om@nyears))
  MSE <- runMSE(om, MPs=MPs)
  testthat::expect_is(MSE, "MSE")
})
  

testthat::test_that("Works with 1 index", {
  om@cpars$Data <- new("Data")
  om@cpars$Data@Type <- c("VBiomass")
  om@cpars$Data@RInd <- array(VBind, dim=c(1, 1, om@nyears))
  MSE <- runMSE(om, MPs=MPs)
  testthat::expect_is(MSE, "MSE")
})




