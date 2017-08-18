
testthat::context("Test available MPs can be applied to the Data objects")

library(DLMtool)
library(DLMdata)
datobjs <- avail('Data')
input <- avail("Input")
output <- avail("Output")

Ntest <- length(datobjs)
canMPs <- list()

testthat::test_that("Can function works with all built-in Data objects", {
  for (x in 1:Ntest) {
    dat <- get(datobjs[x])
    testthat::expect_error(canMPs[[x]] <<- Can(dat), NA, info=dat@Name)
  }
}) 

testthat::test_that("Cant function works with all built-in Data objects", {
  for (x in 1:Ntest) {
    dat <- get(datobjs[x])
    testthat::expect_error(Cant(dat), NA, info=dat@Name)
  }
})
 
testthat::test_that("Needed function works with all built-in Data objects", {
  for (x in 1:Ntest) {
    dat <- get(datobjs[x])
    testthat::expect_error(Cant(dat), NA, info=dat@Name)
  }
})


testthat::test_that("Output control works with all built-in Data objects", {
  for (x in 1:Ntest) {
    # if(interactive()) print(x)
    dat <- get(datobjs[x])
    MPs <- canMPs[[x]]
    for (mm in MPs) {
      # if(interactive()) print(mm)
      if (class(get(mm)) == "Output") {
        info <- paste0(dat@Name, ": ", mm)
        testthat::expect_error(TAC(dat, MPs=mm), NA, info=info)
      }  
    }
  }
})

testthat::test_that("Input control works with all built-in Data objects", {
  for (x in 1:Ntest) {
    # if(interactive()) print(x)
    dat <- get(datobjs[x])
    MPs <- canMPs[[x]]
    for (mm in MPs) {
      if (class(get(mm)) == "Input") {
        # if(interactive()) print(mm)
        info <- paste0(dat@Name, ": ", mm)
        testthat::expect_error(Input(dat, MPs=mm, msg=FALSE, CheckMPs = FALSE), NA, info=info)
      }  
    }
  }
})


