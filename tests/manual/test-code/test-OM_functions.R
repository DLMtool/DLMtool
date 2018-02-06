testthat::context("test OM_functions")


testthat::test_that("ForceCor works ", {
  OM <- DLMtool::testOM
  OM@nsim <- 6 
  OM <- ForceCor(OM)
  testthat::expect_error(runMSE(OM, silent=TRUE), NA)
})

Stock <- Bluefin_tuna
Fleet <- DecE_NDom
Obs <- Imprecise_Biased
Imp <- Overages
  
types <- c("Stock", "Fleet", "Obs", "Imp")
for (type in types) {
  testthat::test_that(paste("Replace and SubOM work with ", type), {
    OM1 <- DLMtool::testOM
    OM2 <- new("OM", Stock, Fleet, Obs, Imp)
    OMnew <- Replace(OM1, OM2, type)
    
    newobj <- SubOM(OMnew, type)
    nms <- slotNames(newobj)
    nms <- nms[!nms %in% c("Name", "Species", "Region", "Agency", "Latitude", "Longitude", "Source")]
    for (nm in nms) {
      chk <- is.na(slot(OM2, nm))
      slot(OM2, nm)[chk] <- rep(0, length(chk))
      chk <- is.na(slot(newobj, nm))
      slot(newobj, nm)[chk] <- rep(0, length(chk))
      testthat::expect_equal(slot(newobj, nm), slot(OM2, nm))
    }
  })
}
  
  
  
 
                    

