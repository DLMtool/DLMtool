testthat::context("Data_Functions")

DLMextra()
library(DLMextra)
Dat <- avail('Data')
Dat <- Dat[!Dat %in% c("SimulatedData", "Simulation_1")]

output <- avail('Output')
                    
for (dat in Dat) {
  testthat::test_that(paste("runMP, Can, Cant, Needed, Sense, TAC, work with ", dat), {
    datobj <- get(dat)
    cans <- Can(datobj)
    cants <- Cant(datobj)
    avails <- avail("MP")
    testthat::expect_error(out <<- runMP(datobj, silent=TRUE), NA)
    testthat::expect_true(length(out@MPs) == length(cans))
    testthat::expect_true(length(avails) - nrow(cants) == length(cans))
    testthat::expect_true(length(Needed(datobj)) == nrow(cants))
    testthat::expect_error(Sense(out, out@MPs[1]), NA)
    if (dat %in% cans) testthat::expect_error(TAC(datobj), NA)
  })
}
 

  
                    

