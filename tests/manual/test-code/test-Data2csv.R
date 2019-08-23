
testthat::context("Test Data2csv function")
rm(list=ls())

dats<-avail('Data')
file <- "test.csv"
sim <- 1
for (dat in dats) {
  testthat::test_that(paste("Data2csv works with ", dat), {
    Data2csv(get(dat), file, simno = sim,overwrite=T)
    readDat <- new("Data", file)
    for (sl in slotNames('Data')) {
      if (!sl %in% c("TAC", "Sense", "MPrec")) {
        orig <- slot(get(dat), sl)
        if (class(orig) == "integer") orig <- as.numeric(orig)
        read <- slot(readDat, sl)
        testthat::expect_equal(class(orig), class(read))
        if (class(orig) == "character") {
          och <- nchar(orig)
          if (length(och)<1) och <- 0
          testthat::expect_equal(och, nchar(read))
        } else {
          if (class(orig)=="matrix") {
            nonna <- which(!is.na(orig[sim,]))
            testthat::expect_equal(orig[sim,nonna], read[sim,nonna])
          }
          if (class(orig)=="numeric") testthat::expect_equal(orig[sim], read[sim])
        }
      }
    }
  })
}

file.remove(file)






