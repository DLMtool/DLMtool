
testthat::context("Test available MPs can be applied to the Data objects")

library(DLMtool)
DLMextra()
library(DLMextra)
datobjs <- avail('Data')
input <- avail("Input")
output <- avail("Output")

Ntest <- length(datobjs)
canMPs <- list()
for (x in 1:Ntest) {
  dat <- get(datobjs[x])
  testthat::test_that(paste0("Can function works with ", dat@Name), {
    testthat::expect_is(canMPs[[x]] <<- Can(dat), 'character', info=dat@Name)
  }) 
}


for (x in 1:Ntest) {
  dat <- get(datobjs[x])
  testthat::test_that(paste0("Cant function works with ", dat@Name), {
    testthat::expect_is(Cant(dat), 'matrix', info=dat@Name)
  }) 
}


for (x in 1:Ntest) {
  dat <- get(datobjs[x])
  testthat::test_that(paste0("Needed function works with ", dat@Name), {
    testthat::expect_is(Needed(dat), 'character', info=dat@Name)
  }) 
}


for (x in 1:Ntest) {
  dat <- get(datobjs[x])
  testthat::test_that(paste0("Output control works with ", dat@Name), {
    MPs <- canMPs[[x]]
    for (mm in MPs) {
      if (mm %in%  output) {
        info <- paste0(dat@Name, ": ", mm)
        testthat::expect_error(TAC(dat, MPs=mm), NA, info=info)
      } 
    }
  }) 
}

for (x in 1:Ntest) {
  dat <- get(datobjs[x])
  testthat::test_that(paste0("Input control works with ", dat@Name), {
    MPs <- canMPs[[x]]
    for (mm in MPs) {
      if (mm %in%  input) {
        info <- paste0(dat@Name, ": ", mm)
        testthat::expect_error(Input(dat, MPs=mm), NA, info=info)
      } 
    }
  }) 
}





