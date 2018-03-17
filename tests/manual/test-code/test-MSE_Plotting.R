testthat::context("Test MSE Plotting functions")
windows()
DLMextra()
library(DLMextra)
rm(list=ls())
# setup()
MSEobj <- updateMSE(DLMextra::testMSE)

funs <- plotFun(msg=FALSE)
funs <- funs[!funs %in% c('plotFleet', 'plotStock', 
                          "COSEWIC_plot", "DFO_hist",
                          'plotOFL')]

for (ff in funs) {
  testthat::test_that("main plot MSE functions", {
    fun <- get(ff)
    testthat::expect_error(fun(MSEobj), NA, info=ff)
  })
}




