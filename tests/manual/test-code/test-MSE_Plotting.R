testthat::context("Test MSE Plotting functions")
windows()
DLMextra()
library(DLMextra)
rm(list=ls())
# setup()
MSEobj <- updateMSE(DLMextra::testMSE)
COSEWICobj <- runCOSEWIC(testOM)
funs <- plotFun(msg=FALSE)

for (ff in funs) {
  testthat::test_that("main plot MSE functions", {
    fun <- get(ff)
    if (grepl("COSEWIC", ff)) {
      testthat::expect_error(fun(COSEWICobj), NA, info=ff)  
    } else {
      if ('Show' %in%names(formals(fun))) {
        testthat::expect_error(fun(MSEobj, Show=FALSE), NA, info=ff) 
      } else {
        testthat::expect_error(fun(MSEobj), NA, info=ff) 
      }
       
    }
    
  })
}

