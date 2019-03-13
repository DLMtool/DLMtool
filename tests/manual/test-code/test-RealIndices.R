
testthat::context("test Real indices in Data object")

# Simulate some data # 
library(DLMtool)
OM <- DLMtool::testOM
OM@nsim <- 5
Hist <- runMSE(OM, Hist=TRUE, silent=TRUE)
MPs <- "AvC"
# Grab indices frOM one sim 
sim <- sample(1:OM@nsim,1)
Bind <- Hist@TSdata$B[sim,]
VBind <- Hist@TSdata$VB[sim,]
SpBInd<- Hist@TSdata$SSB[sim,]


testthat::test_that("Works with 3 indices", {
  OM@cpars$Data <- new("Data")
  OM@cpars$Data@Type <- c("Biomass", "VBiomass", 'SpBiomass')
  OM@cpars$Data@RInd <- array(rbind(Bind, VBind, SpBInd), dim=c(1, 3, OM@nyears))
  MSE <- runMSE(OM, MPs=MPs, silent = TRUE)
  testthat::expect_is(MSE, "MSE")
})
  
testthat::test_that("Works with 2 indices", {
  OM@cpars$Data <- new("Data")
  OM@cpars$Data@Type <- c("Biomass", 'SpBiomass')
  OM@cpars$Data@RInd <- array(rbind(Bind, SpBInd), dim=c(1, 2, OM@nyears))
  MSE <- runMSE(OM, MPs=MPs, silent = TRUE)
  testthat::expect_is(MSE, "MSE")
})
  

testthat::test_that("Works with 1 index", {
  OM@cpars$Data <- new("Data")
  OM@cpars$Data@Type <- c("VBiomass")
  OM@cpars$Data@RInd <- array(VBind, dim=c(1, 1, OM@nyears))
  MSE <- runMSE(OM, MPs=MPs, silent = TRUE)
  testthat::expect_is(MSE, "MSE")
})


testthat::test_that("Import 3 real indices from Excel", {
  Data <- new("Data")
  Data@Type <- c("Biomass", "VBiomass", 'SpBiomass')
  Data@RInd <- array(rbind(Bind, VBind, SpBInd), dim=c(1, length(Data@Type), OM@nyears))
  Data2csv(Data, file='test.csv', overwrite = TRUE)
  DataIn <- new("Data", "test.csv")
  testthat::expect_equal(DataIn@Type, Data@Type)
  testthat::expect_equal(DataIn@RInd[1,,], Data@RInd[1,,])
})

testthat::test_that("Import 2 real indices from Excel", {
  Data <- new("Data")
  Data@Type <- c("Biomass", 'SpBiomass')
  Data@RInd <- array(rbind(Bind, SpBInd), dim=c(1, length(Data@Type), OM@nyears))
  Data2csv(Data, file='test.csv', overwrite = TRUE)
  DataIn <- new("Data", "test.csv")
  testthat::expect_equal(DataIn@Type, Data@Type)
  testthat::expect_equal(DataIn@RInd[1,,], Data@RInd[1,,])
})

testthat::test_that("Import 1 real index from Excel", {
  Data <- new("Data")
  Data@Type <- c('SpBiomass')
  Data@RInd <- array(SpBInd, dim=c(1, length(Data@Type), OM@nyears))
  Data2csv(Data, file='test.csv', overwrite = TRUE)
  DataIn <- new("Data", "test.csv")
  testthat::expect_equal(DataIn@Type, Data@Type)
  testthat::expect_equal(DataIn@RInd[1,,], Data@RInd[1,,])
})

tt <- file.remove('test.csv')


testthat::test_that("Works with 2 Real indices longer than nyears", {
  om <- OM
  om@nyears <- 30
  om@cpars$Data <- new("Data")
  om@cpars$Data@Type <- c("Biomass", 'SpBiomass')
  om@cpars$Data@RInd <- array(rbind(Bind, SpBInd), dim=c(1, 2, OM@nyears))
  MSE <- runMSE(om, MPs=MPs, silent = TRUE, PPD=TRUE)
  Dataout <- MSE@Misc$Data[[1]]
  testthat::expect_equal(Dataout@RInd[1,1,1:50], om@cpars$Data@RInd[1,1,])
  testthat::expect_equal(Dataout@RInd[1,2,1:50], om@cpars$Data@RInd[1,2,])
})

testthat::test_that("Works with 1 Real index longer than nyears", {
  om <- OM
  om@nyears <- 30
  om@cpars$Data <- new("Data")
  om@cpars$Data@Type <- "Biomass"
  om@cpars$Data@RInd <- array(Bind, dim=c(1, 1, OM@nyears))
  MSE <- runMSE(om, MPs=MPs, silent = TRUE, PPD=TRUE)
  Dataout <- MSE@Misc$Data[[1]]
  testthat::expect_equal(Dataout@RInd[1,1,1:50], om@cpars$Data@RInd[1,1,])
})

testthat::test_that("Works with 1 Real index shorter than nyears", {
  om <- OM
  om@nyears <- 60
  om@cpars$Data <- new("Data")
  om@cpars$Data@Type <- "Biomass"
  om@cpars$Data@RInd <- array(Bind, dim=c(1, 1, OM@nyears))
  MSE <- runMSE(om, MPs=MPs, silent = TRUE, PPD=TRUE)
  Dataout <- MSE@Misc$Data[[1]]
  testthat::expect_equal(Dataout@RInd[1,1,1:50], om@cpars$Data@RInd[1,1,])
})




