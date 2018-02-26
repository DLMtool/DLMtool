testthat::context("test plot Data")

DLMextra(TRUE)
library(DLMextra)
Dat <- avail(Data)
Dat <- Dat[!Dat %in% c("SimulatedData", "Simulation_1")]
                    
for (dat in Dat) {
  testthat::test_that(paste("plot works with ", dat), {
    datobj <- get(dat)
    MP <- Can(datobj)
    outputMPs <- MP[which(MPtype(MP)[,2] == "Output")]
    datobj <- TAC(datobj, outputMPs)
    testthat::expect_error(plot(datobj), NA)
  })
}
 
   

