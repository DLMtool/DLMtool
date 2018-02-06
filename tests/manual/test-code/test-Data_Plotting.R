testthat::context("test plot Data")

DLMextra()
library(DLMextra)
Dat <- avail(Data)
Dat <- Dat[!Dat %in% c("SimulatedData", "Simulation_1")]
                    
for (dat in Dat) {
  testthat::test_that(paste("plot works with ", dat), {
    datobj <- get(dat)
    testthat::expect_error(plot(datobj), NA)
  })
}
 
                   
                    

