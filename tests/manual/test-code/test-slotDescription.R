testthat::context("slot descriptions")


classes <- c("Stock", "Fleet", "Obs", "Imp", "Data", "OM")

for (cl in classes) {
  Desc <- get(paste0(cl, "Description"))
  slots <- slotNames(cl)
  testthat::test_that(paste("slots are in correct order for ", cl), {
    testthat::expect_true(all(slots == Desc[,1]))
  })
}

cbind(slots, Desc[,1], slots == Desc[,1])

cbind(sort(slots), sort(Desc[,1]), sort(slots) == sort(Desc[,1]))


               
                    
