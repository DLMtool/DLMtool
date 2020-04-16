testthat::context("cpars works for all Stock, Fleet, Obs, Imp slots")

rm(list=ls())
library(DLMtool)

# Stock Parameters
OM <- testOM 
OM@nsim <- 5 

slots <- slotNames("Stock")

for (sl in slots) {
  testthat::test_that(paste("cpars works with Stock slot: ", sl), {

    OM2 <- OM 
    sl_val <- suppressWarnings(mean(slot(OM2, sl)))
    if (sl == "Perr") {
      OM2@cpars[[sl]] <- matrix(sl_val[1], OM2@nsim, OM2@nyears+OM2@proyears+OM2@maxage)
    } else {
      OM2@cpars[[sl]] <- rep(sl_val[1], OM2@nsim)  
    }
    
    testthat::expect_error(runMSE(OM2, Hist=TRUE), NA)
  })
}


# Fleet Parameters
OM <- testOM 
OM@nsim <- 5 

slots <- slotNames("Fleet")

for (sl in slots) {
  testthat::test_that(paste("cpars works with Fleet slot: ", sl), {
    OM2 <- OM 
    sl_val <- suppressWarnings(mean(slot(OM2, sl)))
    if(is.na(sl_val)) sl_val <- 0.05
    OM2@cpars[[sl]] <- rep(sl_val[1], OM2@nsim)
    if (sl %in% c("L5", "LFS", "Vmaxlen", "LR5", "LFR", "Rmaxlen")) {
      OM2@isRel <- 'FALSE'
      testthat::expect_error(runMSE(OM2, Hist=TRUE), NA)  
    } else {
      if (sl %in% c("EffYears", "EffUpper", "EffLower")) {
        testthat::expect_error(runMSE(OM2, Hist=TRUE))
      } else {
        testthat::expect_error(runMSE(OM2, Hist=TRUE), NA)  
      }
    }
  
  })
}

testthat::test_that(paste("cpars works with Fleet slots: EffLower, EffUpper, EffYears"), {
  OM2 <- OM
  OM2@cpars$EffYears <- 1:20 
  OM2@cpars$EffLower <- seq(0, 1, length.out=length(OM2@cpars$EffYears))
  OM2@cpars$EffUpper <- OM2@cpars$EffLower + 0.1
  testthat::expect_error(runMSE(OM2, Hist=TRUE), NA)  
})


