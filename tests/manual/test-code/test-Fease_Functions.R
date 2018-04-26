testthat::context("Feasibility Functions")

Data <- avail("Data")
data <- sample(Data, 1)

testthat::test_that(paste("Fease works with ", data), {
  data <- get(data)
  Grid <- expand.grid(TAC=c(TRUE, FALSE), TAE=c(TRUE, FALSE), SL=c(TRUE, FALSE), Spatial=c(TRUE, FALSE), names.only=c(TRUE, FALSE))
  
  for (rr in 1:nrow(Grid)) {
    TAC <- Grid$TAC[rr]
    TAE <- Grid$TAE[rr]
    SL <- Grid$SL[rr]
    Spatial <- Grid$Spatial[rr]
    names.only <- Grid$names.only[rr]
    
    if (all(c(TAC, TAE, SL, Spatial) == FALSE)) {
      testthat::expect_error(Fease(data, TAC, TAE, SL, Spatial, names.only, msg=FALSE), info=Grid[rr,])
    } else {
      testthat::expect_error(Fease(data, TAC, TAE, SL, Spatial, names.only, msg=FALSE), NA, info=Grid[rr,])
      if (names.only) {
        testthat::expect_is(Fease(data, TAC, TAE, SL, Spatial, names.only, msg=FALSE), "character", info=Grid[rr,])
      } else {
        testthat::expect_is(Fease(data, TAC, TAE, SL, Spatial, names.only, msg=FALSE), "data.frame", info=Grid[rr,])
      }
    }
  }
})



